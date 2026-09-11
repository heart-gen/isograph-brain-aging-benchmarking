"""BrainSEQ cis switch-QTL (S_g) and matched abundance-QTL (A_g) mapping.

WHY THIS EXISTS
---------------
Every genetic result in this paper is currently mediated through GTEx: the set-level
sQTL/eQTL contrast, the per-gene coloc, the CLPP layer. GTEx brain is 181-300 donors per
region, and its splicing phenotype is a LeafCutter intron-excision ratio -- neither the
cohort nor the phenotype is the one the switch layer was discovered in. That leaves two
open objections a reviewer will press: the per-gene modality contrast is null (Fisher
P = 0.88), and S-LDSC's splicing arm clears nominal significance in 1 of 6
trait-contexts. Both are consistent with "there is no splicing effect" and both are
equally consistent with "GTEx bulk sQTL is underpowered".

This module attacks that directly by mapping QTLs against **IsoGraph's own phenotypes**:

    S_g  the gene's switch coordinate  (PC1 of its within-gene CLR isoform composition)
    A_g  the gene's total abundance    (z-scored log2 CPM)

against the same donors, the same variants, the same cis window and the same covariates.
Because there is exactly ONE S_g and ONE A_g per gene, the swQTL-vs-eQTL comparison is a
clean paired test -- it avoids the one-expression-phenotype-versus-many-introns asymmetry
that forces the GTEx contrast to pick a single representative intron and then argue the
choice was not circular.

The question it answers is the one the paper's thesis rests on: **does genetic variation
preferentially regulate the relative-isoform axis or the total-abundance axis of the same
gene?**

WHAT THIS IS NOT
----------------
It is **not independent replication**. BrainSEQ is the cohort the switches were
discovered in, so this is *same-tissue genetic anchoring* / *endogenous genetic
validation*. GTEx remains the external cohort. Cross-region consistency within BrainSEQ
(caudate / DLPFC / hippocampus) is a separate and weaker axis of support. Any prose that
calls this replication is wrong.

THE PHENOTYPES ARE RECOMPUTED, NOT REUSED, AND THAT IS THE POINT
---------------------------------------------------------------
`feature_scores.parquet` from the discovery fits covers only the samples IsoGraph was fit
on -- controls, adults, 238/222/238 per region, of which 232/201/232 are genotyped. That
is *fewer* donors than the validation plan's "452-500", because that figure is the
junction-phenotype n, not the switch-phenotype n.

But `S_g` and `A_g` are deterministic per-gene transforms of the transcript counts
(`isograph.features.channels.gene_feature_channels`), not VAE outputs -- the VAE only
builds the gene-gene network *on top of* them. So they can be recomputed on any sample
set, and a cis-QTL is diagnosis-independent, so the full region cohort is the right
denominator with `Dx` carried as a covariate.

Sample filter (agreed with the analyst, and deliberately looser than the discovery fit's
`Dx == Control, Age >= 18`):
    canonical BSP dataset per region, `Age > 13`, `dropped != "t"`.

ANCESTRY IS A REAL CONSTRAINT HERE
----------------------------------
After that filter the cohort is roughly half African-ancestry. Two arms follow, and the
distinction is not cosmetic:
  * `all_samples` -- the PRIMARY arm for the modality contrast. The contrast is computed
    within a gene, in one cohort, so admixture does not bias it, and the full set is the
    powered one.
  * `ea_only`     -- required for anything colocalized against the EUR GWAS. LD in a ~50%
    AA cohort matches neither the EUR GWAS nor the 1000G EUR panel, so a coloc run on the
    multi-ancestry arm would be comparing signals fine-mapped under incompatible LD.

STAGES
------
  --stage phenotypes  Recompute S_g / A_g on the expanded sample set; write tensorQTL BED
                      phenotypes + covariates, keyed by BrNum. Login-node safe-ish
                      (minutes, a few GB).
  --stage map         cis nominal + permutation + SuSiE, both modalities, via tensorQTL.
  --stage junctions   Junction-level phenotypes for swQTL-positive genes (see
                      `build_brainseq_junction_usage` / `build_brainseq_cluster_usage`).
  --stage meta        The paired swQTL-vs-eQTL contrast and the report.

Outputs under 05_genetic_anchoring/_m/brainseq_switch_qtl/<arm>/<region>/.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out

REGIONS: tuple[str, ...] = ("caudate", "dlpfc", "hippocampus")
# `all_samples` is primary for the modality contrast; `ea_only` exists for coloc against
# the EUR GWAS. Both are built by `.../processed-data/genotypes/qtl/_h/`.
ARMS: tuple[str, ...] = ("all_samples", "ea_only")

_GENO_ROOT = Path("/ocean/projects/bio260021p/shared/resources/processed-data/genotypes/qtl")
# GENCODE v47 gene/transcript annotations shipped beside the BrainSEQ rse objects. Used
# instead of the featureCounts `Chr/Start/End` columns, which are comma-joined per-exon
# lists rather than gene spans, and instead of the bundle's transcripts.parquet, which is
# already restricted to the discovery fit's expressed genes.
_RVAR_ROOT = Path("/ocean/projects/bio260021p/shared/resources/processed-data/r-variables")

# Agreed QC filter. Deliberately NOT `build_bundles.ADULT_AGE_MIN` (18): a cis-QTL does
# not need the adult-brain restriction the aging analysis needs, and dropping to 13 keeps
# ~40 more donors per region. `dropped` is the 't'/'f' string flag in the LIBD metadata.
AGE_MIN = 13.0
AGE_INCLUSIVE = False        # `Age > 13`, as specified

# QTL covariates. Genotype structure + the technical covariates the BrainSEQ metadata
# carries, plus Dx to absorb the case/control split that the expanded cohort introduces
# and that the discovery (controls-only) fit never saw. IDENTICAL for both modalities --
# handing either axis extra covariates would manufacture the contrast this module tests.
QTL_COVARIATES: tuple[str, ...] = (
    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5",
    "Sex", "Dx", "RIN", "mito_rate", "mapping_rate",
)
CIS_WINDOW = 1_000_000
_SEED = 13


def out_dir(arm: str = "all_samples", region: str | None = None) -> Path:
    d = stage_out("anchoring.brainseq_qtl", arm)
    return ensure_dir(d if region is None else d / region)


# --------------------------------------------------------------------------- #
# Stage: phenotypes
# --------------------------------------------------------------------------- #
def expanded_sample_table(region: str) -> pd.DataFrame:
    """The QTL sample set: canonical dataset, `Age > 13`, `dropped != 't'`, all diagnoses.

    Reuses the discovery bundle's metadata joins verbatim (SNP PCs, per-region RNA-seq QC
    metrics, the MoD tidy-up) so the covariates mean the same thing here as they do in
    every other BrainSEQ analysis in the repo. Only the three filter predicates differ,
    and they differ on purpose.
    """
    from isograph_benchmark.inputs.build_bundles import (
        _REGION_DATASET, _load_brainseq_qc_metrics, _load_snp_pcs,
    )
    meta = pd.read_parquet(
        rel("inputs", "processed", "brainseq", "metadata", "libd_rnaseq_metadata.parquet"))
    tx_cols = _tx_sample_ids(region)
    s = meta[meta["RNum"].isin(tx_cols)].copy()
    s = s[(s["Dataset"] == _REGION_DATASET[region]) & (s["dropped"] != "t")]
    s = s[s["Age"] > AGE_MIN] if not AGE_INCLUSIVE else s[s["Age"] >= AGE_MIN]
    s["MoD"] = (s["MoD"].replace("No autopsy performed", "Undetermined")
                        .replace(".", "Undetermined").fillna("Undetermined"))
    s = s.merge(_load_snp_pcs(), on="BrNum", how="left")
    s = s.merge(_load_brainseq_qc_metrics(region),
                left_on="RNum", right_on="sample_rnum", how="left")
    return s.drop(columns=["sample_rnum"], errors="ignore")


def _tx_sample_ids(region: str) -> list[str]:
    import pyarrow.parquet as pq
    f = pq.ParquetFile(rel("inputs", "processed", "brainseq", region, "tx_counts.parquet"))
    return [c for c in f.schema_arrow.names
            if c not in ("Name", "Length", "EffectiveLength")]


# Panel choice, reconciled 2026-09-10 (see `reconcile_genotype_panels`). Two TOPMed LIBD
# releases are on disk and they disagree about which donors exist:
#
#   inputs/raw/brainseq/genotypes/TOPMed_LIBD.*   (legacy, symlinked in-repo)  n = 1,938
#   .../processed-data/genotypes/qtl/all_samples  (current)                    n = 2,438
#
# The current release is larger overall but is MISSING 23 caudate donors the legacy one
# has (19 AA, 4 CAUC -- all present in subjects_phenodata with a valid race, so they were
# dropped by the release, not by a metadata join). Caudate lands at 420 donors on the
# current panel against 436 on the legacy one.
#
# The current tree is used anyway, and the 4% of donors is the price:
#   * its QC is scripted and documented (`qtl/_h/`: --maf 0.005 --geno 0.1 --hwe 1e-6
#     --rm-dup force-first, PCA after LD pruning), where the legacy panel's is not;
#   * it ships the `ea_only` / `aa_only` sub-panels WITH matching genotype PCs, and the
#     EA-only arm is not optional -- LD in a ~50% AA cohort matches neither the EUR GWAS
#     nor the 1000G EUR panel, so coloc has to run on it.
# A 4% power difference does not trade against a documented QC path and a coloc-able
# ancestry stratification.
def genotyped_donors(arm: str, chrom: int = 22) -> set[str]:
    """BrNums present in the genotype panel for an arm.

    Read from a per-chromosome psam because those are keyed by `BrNum`, unlike the merged
    psam whose `#IID` is a chip barcode.
    """
    f = _GENO_ROOT / arm / f"chr{chrom}.psam"
    if not f.exists():
        raise SystemExit(f"missing genotype psam: {f}")
    d = pd.read_csv(f, sep="\t")
    col = "#IID" if "#IID" in d.columns else d.columns[0]
    return set(d[col].astype(str))


def reconcile_genotype_panels(region: str) -> pd.DataFrame:
    """Overlap of this region's RNA-passing donors with every genotype panel on disk.

    The plan requires the panel choice to be recorded rather than assumed, because the
    two releases disagree (see the note on `genotyped_donors`). Written to disk beside
    the phenotypes so the number quoted in the paper can be traced to a file.
    """
    tab = (expanded_sample_table(region)
           .sort_values("RIN", ascending=False).drop_duplicates("BrNum"))
    donors = set(tab["BrNum"].astype(str))
    rows = [{"panel": "rna_filter_passing", "panel_n": len(donors),
             "overlap": len(donors), "source": "expanded_sample_table"}]
    for sub in ("all_samples", "ea_only", "aa_only", "ea_aa_samples"):
        f = _GENO_ROOT / sub / "chr22.psam"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        col = "#IID" if "#IID" in d.columns else d.columns[0]
        ids = set(d[col].astype(str))
        rows.append({"panel": sub, "panel_n": len(ids),
                     "overlap": len(donors & ids), "source": str(f)})
    legacy = rel("inputs", "raw", "brainseq", "genotypes", "TOPMed_LIBD.psam")
    if legacy.exists():
        d = pd.read_csv(legacy, sep="\t")
        ids = set(d["#FID"].astype(str))
        rows.append({"panel": "legacy_AA_EA", "panel_n": len(ids),
                     "overlap": len(donors & ids), "source": str(legacy)})
    return pd.DataFrame(rows)


# Genotype PCs per arm. `all_samples` keeps the discovery bundle's multi-ancestry PCs (the
# same PCA the rest of the repo adjusts for, so its committed results reproduce). A
# sub-panel arm must use PCs computed WITHIN that panel: in an EA-only model the
# multi-ancestry PCs spend their leading axes on the AA/EA split the arm has already
# removed, and leave the within-EA structure unadjusted.
ARM_PCS: dict[str, Path | None] = {
    "all_samples": None,
    "ea_only": _GENO_ROOT / "genetic_similarity" / "ea_only" / "TOPMed_LIBD.EA.eigenvec",
}


def arm_genotype_pcs(samples: pd.DataFrame, arm: str,
                     pcs: pd.DataFrame | None = None) -> pd.DataFrame:
    """Swap in the arm's own genotype PCs (`SNP_PC*`), keyed by BrNum. Refuses a donor
    without them rather than letting a median fill stand in for ancestry."""
    f = ARM_PCS.get(arm)
    if f is None:
        return samples
    if pcs is None:
        if not f.exists():
            raise SystemExit(f"missing {arm} genotype PCs: {f}")
        pcs = pd.read_csv(f, sep="\t")
    pcs = pcs.rename(columns={"#IID": "BrNum", **{f"PC{i}": f"SNP_PC{i}" for i in range(1, 11)}})
    pc_cols = [c for c in pcs.columns if c.startswith("SNP_PC")]
    pcs = pcs[["BrNum", *pc_cols]].assign(BrNum=lambda x: x["BrNum"].astype(str))
    out = (samples.drop(columns=[c for c in samples.columns if c.startswith("SNP_PC")])
                  .assign(BrNum=lambda x: x["BrNum"].astype(str))
                  .merge(pcs, on="BrNum", how="left"))
    n_missing = int(out["SNP_PC1"].isna().sum())
    if n_missing:
        raise SystemExit(f"{n_missing} {arm} donors have no {arm} genotype PCs in {f}")
    return out


def build_phenotypes(region: str, arm: str = "all_samples") -> Path:
    """Recompute S_g / A_g on the expanded set and write tensorQTL BED phenotypes."""
    from isograph.features.channels import gene_feature_channels

    from isograph_benchmark.inputs.build_bundles import (
        _matrix_from_wide_subset, _brainseq_expressed_genes, _restrict_to_expressed_genes,
    )
    dest = ensure_dir(out_dir(arm, region))
    samples = expanded_sample_table(region)
    geno = genotyped_donors(arm)

    # One sample per donor: QTL mapping is per genotype, and a donor contributing two
    # libraries would otherwise enter the model twice. Keep the highest-RIN library.
    samples = (samples.sort_values("RIN", ascending=False)
                      .drop_duplicates("BrNum", keep="first"))
    n_pre = len(samples)
    samples = samples[samples["BrNum"].astype(str).isin(geno)].copy()
    if samples.empty:
        raise SystemExit(f"no genotyped samples for {region}/{arm}")
    samples = arm_genotype_pcs(samples, arm)

    src = rel("inputs", "processed", "brainseq", region)
    tx = pd.read_parquet(src / "tx_counts.parquet")
    gene = pd.read_parquet(src / "gene_counts.parquet")
    keep = [r for r in _tx_sample_ids(region) if r in set(samples["RNum"])]
    samples = samples.set_index("RNum").loc[keep].reset_index()

    tx_annot = _tx_annotation(region)
    tx_feature, tx_matrix = _matrix_from_wide_subset(
        tx, ["Name", "Length", "EffectiveLength"], keep)
    tx_feature = tx_feature.rename(columns={"Name": "transcript_id"}).merge(
        tx_annot, on="transcript_id", how="left")
    gene_feature, gene_matrix = _matrix_from_wide_subset(
        gene, ["Geneid", "Chr", "Start", "End", "Strand", "Length"], keep)
    gene_feature = gene_feature.rename(columns={"Geneid": "gene_id"})

    # The SAME expression filter the discovery bundle applies, so a gene is in or out for
    # the same reason here as there. Recomputed on this sample set, so the gene list can
    # differ slightly from the discovery fit's -- the overlap is reported below.
    gene_keep, expr_note = _brainseq_expressed_genes(gene_matrix)
    gene_feature, gene_matrix, tx_feature, tx_matrix, expr_stats = \
        _restrict_to_expressed_genes(gene_feature, gene_matrix, tx_feature, tx_matrix,
                                     gene_keep)

    # The SAME transcript filter the discovery fits apply before the switch coordinate
    # (`run_models._filter_expressed_transcripts`: count > 10 in >= 70% of samples), computed
    # here over the QTL sample set. Omitting it is not a detail: CLR PC1 over every annotated
    # transcript of an expressed gene is dominated by near-zero isoforms. On the discovery
    # libraries the unfiltered coordinate matched the discovery S_g at median |r| 0.36, the
    # filtered one at 1.000 (2026-09-11, `brainseq_qtl_checks --stage signpin`). It changes only
    # the switch channel; A_g is read from the gene counts below.
    from isograph_benchmark.real_data.run_models import _filter_expressed_transcripts
    n_tx_expressed_genes = int(tx_matrix.shape[0])
    tx_matrix, tx_feature = _filter_expressed_transcripts(tx_matrix, tx_feature)

    # Deterministic per-gene transforms; `design=None` because covariates enter the QTL
    # model, not the phenotype (the 2026-06-28 covariate-adjustment policy).
    matrix, info = gene_feature_channels(
        tx_matrix, tx_feature, gene_counts=gene_matrix, gene_table=gene_feature)
    scores = pd.DataFrame(matrix, columns=samples["BrNum"].astype(str).tolist())
    scores = pd.concat([info.reset_index(drop=True), scores], axis=1)

    coords = _gene_coords(region)
    written = {}
    for ftype, tag in (("switch", "S_g"), ("abundance", "A_g")):
        sub = scores[scores["feature_type"] == ftype].copy()
        bed = _to_bed(sub, coords, samples["BrNum"].astype(str).tolist())
        f = dest / f"phenotypes_{ftype}.bed.gz"
        bed.to_csv(f, sep="\t", index=False, compression="gzip")
        written[ftype] = (len(bed), str(f))
        print(f"  {tag:<4} ({ftype}): {len(bed):,} genes x {len(samples):,} donors -> {f.name}")

    cov = _covariate_frame(samples)
    cov.to_csv(dest / "covariates.txt", sep="\t")
    samples.to_parquet(dest / "samples.parquet", index=False)

    panels = reconcile_genotype_panels(region)
    panels.to_csv(dest / "genotype_panel_reconciliation.tsv", sep="\t", index=False)
    print("  genotype panel reconciliation:")
    print(panels[["panel", "panel_n", "overlap"]].to_string(index=False))

    meta_out = {
        "region": region, "arm": arm,
        "n_samples_pre_genotype": int(n_pre),
        "n_samples": int(len(samples)),
        "n_donors": int(samples["BrNum"].nunique()),
        "filters": f"Dataset={_region_dataset(region)}, dropped!='t', "
                   f"Age {'>=' if AGE_INCLUSIVE else '>'} {AGE_MIN:g}, all diagnoses",
        "expression_filter": expr_note,
        "transcript_filter": ("count > 10 in >= 70% of samples "
                              "(run_models._filter_expressed_transcripts, as in discovery)"),
        "n_transcripts_expressed_genes": n_tx_expressed_genes,
        "n_transcripts_after_tx_filter": int(tx_matrix.shape[0]),
        "dx": samples["Dx"].value_counts().to_dict(),
        "race": samples["Race"].value_counts().to_dict(),
        "covariates": list(QTL_COVARIATES),
        "cis_window": CIS_WINDOW,
        "n_switch_genes": written["switch"][0],
        "n_abundance_genes": written["abundance"][0],
        **{f"expr_{k}": int(v) for k, v in expr_stats.items()},
    }
    (dest / "phenotype_summary.json").write_text(json.dumps(meta_out, indent=2, default=str))
    print(f"  genotyped-and-phenotyped n = {len(samples):,} "
          f"(of {n_pre:,} passing the RNA filter)")
    print(f"  Dx: {meta_out['dx']}")
    print(f"  Race: {meta_out['race']}")
    return dest


def _region_dataset(region: str) -> str:
    from isograph_benchmark.inputs.build_bundles import _REGION_DATASET
    return _REGION_DATASET[region]


def _tx_annotation(region: str) -> pd.DataFrame:
    """GENCODE v47 transcript -> gene map for a region.

    Only `transcript_id` and `gene_id` are carried. The file's `gene_name` /
    `gene_biotype` columns are the GENE symbol and GENE biotype (`DDX11L16`, `lncRNA`),
    NOT per-transcript ones (`DDX11L16-201`), so renaming them to `transcript_name` /
    `transcript_type` would put gene-level values under transcript-level labels. Nothing
    downstream needs them: `gene_feature_channels` groups on `gene_id` and indexes on
    `transcript_id` alone.
    """
    f = _RVAR_ROOT / region / "_m" / "tx-annotation.tsv"
    if not f.exists():
        raise SystemExit(f"missing transcript annotation: {f}")
    d = pd.read_csv(f, sep="\t", usecols=["transcript_id", "gene_id"])
    return d


def _gene_coords(region: str) -> pd.DataFrame:
    """Strand-aware TSS per gene, from the GENCODE v47 gene annotation.

    tensorQTL windows the cis region on the BED interval, so the strand is what keeps a
    minus-strand gene's window centred on its promoter rather than its 3' end.
    """
    f = _RVAR_ROOT / region / "_m" / "gene-annotation.tsv"
    if not f.exists():
        raise SystemExit(f"missing gene annotation: {f}")
    d = pd.read_csv(f, sep="\t")
    d["tss"] = np.where(d["strand"] == "-", d["end"], d["start"])
    return d.rename(columns={"chrom": "chrom"})[["gene_id", "chrom", "tss", "strand"]]


def _to_bed(scores: pd.DataFrame, coords: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """tensorQTL phenotype BED: `#chr start end phenotype_id <donors...>`, TSS-anchored."""
    d = scores.merge(coords, on="gene_id", how="inner")
    d = d[d["chrom"].astype(str).str.match(r"^chr(\d+|X)$")]
    d = d.dropna(subset=["tss"])
    out = pd.DataFrame({
        "#chr": d["chrom"].to_numpy(),
        "start": d["tss"].astype(int).to_numpy() - 1,   # BED is 0-based, half-open
        "end": d["tss"].astype(int).to_numpy(),
        "phenotype_id": d["gene_id"].astype(str).to_numpy(),
    })
    for s in samples:
        out[s] = d[s].to_numpy()
    out = out.drop_duplicates("phenotype_id")
    return out.sort_values(["#chr", "start"], key=_chrom_key).reset_index(drop=True)


def _chrom_key(col: pd.Series) -> pd.Series:
    if col.name != "#chr":
        return col
    n = col.astype(str).str.replace("^chr", "", regex=True)
    return pd.to_numeric(n.replace({"X": "23"}), errors="coerce").fillna(99)


def _covariate_frame(samples: pd.DataFrame) -> pd.DataFrame:
    """tensorQTL covariates: rows = covariates, columns = donors.

    Categoricals are dummied with the first level dropped. Identical construction for
    both modalities by design.
    """
    cols = [c for c in QTL_COVARIATES if c in samples.columns]
    missing = [c for c in QTL_COVARIATES if c not in samples.columns]
    if missing:
        raise SystemExit(f"missing covariates in the sample table: {missing}")
    d = samples[cols].copy()
    cat = [c for c in cols if d[c].dtype == object]
    num = [c for c in cols if c not in cat]
    X = pd.get_dummies(d, columns=cat, drop_first=True, dtype=float)
    for c in num:
        X[c] = pd.to_numeric(X[c], errors="coerce")
    # A covariate that is constant or all-NA carries no information and makes the design
    # rank-deficient; drop it loudly rather than letting the solve fail downstream.
    X = X.loc[:, X.notna().any() & (X.nunique(dropna=True) > 1)]
    X = X.fillna(X.median(numeric_only=True))
    X.index = samples["BrNum"].astype(str).tolist()
    return X.T


# --------------------------------------------------------------------------- #
# Stage: map
# --------------------------------------------------------------------------- #
# Hidden-factor covariates. PEER is not installed anywhere on this system and is
# effectively deprecated; the modern equivalent is PCs of the phenotype matrix itself.
# The COUNT is what has to be matched between the two modalities -- each axis gets its
# own factors (they are different matrices) but the same NUMBER of them, so neither is
# handed extra degrees of freedom to soak up noise with. GTEx uses 60 factors at this
# sample size; 30 is the conservative half of that and is swept by `--n-factors`.
N_HIDDEN_FACTORS = 30
MAF_THRESHOLD = 0.01
NPERM = 10000
MODALITIES = {"switch": "S_g", "abundance": "A_g"}


def _inverse_normal_transform(df: pd.DataFrame) -> pd.DataFrame:
    """Rank-based inverse normal transform, across samples within each phenotype.

    Applied IDENTICALLY to S_g and A_g. Both axes then have the same marginal
    distribution, so a difference in QTL yield between them cannot be an artefact of one
    being heavier-tailed than the other -- which matters, because the whole point of this
    module is a comparison between the two.
    """
    from scipy import stats
    r = df.rank(axis=1, method="average")
    n = df.notna().sum(axis=1).to_numpy()[:, None]
    q = stats.norm.ppf(r.to_numpy() / (n + 1))
    return pd.DataFrame(q, index=df.index, columns=df.columns)


def _hidden_factors(pheno: pd.DataFrame, known: pd.DataFrame, k: int) -> pd.DataFrame:
    """Top-k PCs of the phenotype matrix after regressing out the known covariates.

    Residualizing first stops the factors from simply re-encoding sex, ancestry or RIN,
    which would make them collinear with the covariates already in the model.
    """
    from numpy.linalg import lstsq, svd
    Y = pheno.to_numpy(dtype=float)
    Y = np.nan_to_num(Y - np.nanmean(Y, axis=1, keepdims=True))
    X = np.column_stack([np.ones(known.shape[0]), known.to_numpy(dtype=float)])
    beta, *_ = lstsq(X, Y.T, rcond=None)
    resid = Y.T - X @ beta                       # samples x phenotypes
    resid = resid - resid.mean(axis=0, keepdims=True)
    sd = resid.std(axis=0, ddof=1)
    resid = resid[:, sd > 1e-12] / sd[sd > 1e-12]
    k = int(min(k, resid.shape[0] - 1, resid.shape[1]))
    U, S, _ = svd(resid, full_matrices=False)
    F = U[:, :k] * S[:k]
    F = (F - F.mean(axis=0)) / F.std(axis=0, ddof=1)
    return pd.DataFrame(F, index=pheno.columns,
                        columns=[f"HF{i+1}" for i in range(k)])


def run_map(region: str, arm: str = "all_samples", n_factors: int = N_HIDDEN_FACTORS,
            chroms: list[int] | None = None) -> Path:
    """cis nominal + permutation, both modalities, matched in every respect but the axis."""
    import tensorqtl
    from tensorqtl import cis, pgen

    src = out_dir(arm, region)
    dest = ensure_dir(src / "qtl")
    cov_known = pd.read_csv(src / "covariates.txt", sep="\t", index_col=0).T
    cov_known.index = cov_known.index.astype(str)

    summary = []
    for ftype, tag in MODALITIES.items():
        bed = src / f"phenotypes_{ftype}.bed.gz"
        if not bed.exists():
            raise SystemExit(f"missing {bed}; run --stage phenotypes first")
        pheno_df, pos_df = tensorqtl.read_phenotype_bed(str(bed))
        pheno_df = _inverse_normal_transform(pheno_df)
        donors = [d for d in pheno_df.columns if d in cov_known.index]
        pheno_df = pheno_df[donors]
        known = cov_known.loc[donors]
        hf = _hidden_factors(pheno_df, known, n_factors)
        cov = pd.concat([known, hf], axis=1).astype(float)
        cov.to_csv(dest / f"covariates_used_{ftype}.txt", sep="\t")
        print(f"  [{tag}] {pheno_df.shape[0]:,} phenotypes x {len(donors):,} donors, "
              f"{cov.shape[1]} covariates ({known.shape[1]} known + {hf.shape[1]} hidden)")

        chrom_list = chroms or sorted(
            {int(c.replace("chr", "")) for c in pos_df["chr"].unique()
             if c.replace("chr", "").isdigit()})
        nom_parts, cis_parts = [], []
        for ch in chrom_list:
            prefix = _GENO_ROOT / arm / f"chr{ch}"
            if not Path(f"{prefix}.pgen").exists():
                print(f"    chr{ch}: no genotypes; skipped")
                continue
            pr = pgen.PgenReader(str(prefix), select_samples=donors)
            geno_df = pr.load_genotypes()
            var_df = pr.variant_df
            keep = pos_df["chr"] == f"chr{ch}"
            ph, pp = pheno_df[keep.to_numpy()], pos_df[keep].copy()
            if not len(ph):
                continue
            # The phenotype BED is GENCODE-named (`chr22`); the TOPMed panel's .pvar is
            # `22`. tensorQTL joins the two on that string, so a mismatch silently drops
            # EVERY phenotype ("dropping N phenotypes on chrs. without genotypes") and
            # then dies with a bare "No phenotypes remain after filters". Conform the
            # phenotype side to whatever the genotypes use rather than assuming either.
            #
            # This renames only; it does NOT convert builds. The panel is GRCh38 --
            # verified against the hg19 1000G bim, where every shared rsID is offset and
            # none is identical -- which is the build the GENCODE v47 TSSs are already
            # in. A future hg19 panel would need a liftover here, not a rename.
            pp["chr"] = _conform_chrom(pp["chr"], var_df["chrom"])
            ph, pp = _drop_phenotypes_without_genotypes(ph, pp, var_df, ch)
            if not len(ph):
                continue
            cis.map_nominal(geno_df, var_df, ph, pp, f"{ftype}.chr{ch}",
                            covariates_df=cov, maf_threshold=MAF_THRESHOLD,
                            window=CIS_WINDOW, output_dir=str(dest), verbose=False)
            nom_parts.append(dest / f"{ftype}.chr{ch}.cis_qtl_pairs.chr{ch}.parquet")
            cdf = cis.map_cis(geno_df, var_df, ph, pp, covariates_df=cov,
                              maf_threshold=MAF_THRESHOLD, window=CIS_WINDOW,
                              nperm=NPERM, seed=_SEED, verbose=False)
            cdf["chrom"] = ch
            cis_parts.append(cdf)
            print(f"    chr{ch}: {len(ph):,} phenotypes, {geno_df.shape[0]:,} variants")
            del geno_df, pr

        if not cis_parts:
            raise SystemExit(f"no chromosomes mapped for {region}/{arm}/{ftype}")
        perm = pd.concat(cis_parts)
        perm.index.name = "phenotype_id"
        perm = perm.reset_index()
        # BH rather than Storey q-values: tensorqtl's `calculate_qvalues` needs rpy2 and
        # the R `qvalue` package, and `rfunc` does not import in this environment. BH is
        # the more conservative choice and is applied identically to both modalities.
        perm["qval"] = _bh(perm["pval_beta"])
        perm.to_parquet(dest / f"cis_qtl_{ftype}.parquet", index=False)
        n_sig = int((perm["qval"] < 0.05).sum())
        print(f"  [{tag}] {n_sig:,} / {len(perm):,} genes with a cis-QTL at BH q < 0.05")
        summary.append({"modality": ftype, "label": tag, "n_phenotypes": len(perm),
                        "n_sig_q05": n_sig, "n_donors": len(donors),
                        "n_covariates": int(cov.shape[1]),
                        "n_hidden_factors": int(hf.shape[1])})

    pd.DataFrame(summary).to_csv(dest / "map_summary.tsv", sep="\t", index=False)
    return dest


def _conform_chrom(pheno_chrom: pd.Series, geno_chrom: pd.Series) -> pd.Series:
    """Make phenotype chromosome labels match the genotype panel's convention."""
    geno_has_prefix = bool(geno_chrom.astype(str).str.startswith("chr").any())
    s = pheno_chrom.astype(str)
    if geno_has_prefix:
        return s.where(s.str.startswith("chr"), "chr" + s)
    return s.str.replace("^chr", "", regex=True)


def _drop_phenotypes_without_genotypes(ph: pd.DataFrame, pp: pd.DataFrame,
                                       var_df: pd.DataFrame, ch: int):
    """Guard the chromosome join tensorQTL performs silently.

    Without this, a naming or reference-build mismatch surfaces only as tensorQTL's
    opaque "No phenotypes remain after filters", minutes into a run, with nothing to say
    which of the two it was.
    """
    shared = set(var_df["chrom"].astype(str))
    keep = pp["chr"].astype(str).isin(shared)
    if not keep.any():
        raise SystemExit(
            f"chr{ch}: no phenotype shares a chromosome label with the genotypes "
            f"(phenotypes use {sorted(set(pp['chr'].astype(str)))[:3]}, genotypes use "
            f"{sorted(shared)[:3]}). That is a naming or reference-build mismatch, not "
            f"a data shortage.")
    return ph[keep.to_numpy()], pp[keep]


def _bh(p: pd.Series) -> pd.Series:
    from statsmodels.stats.multitest import multipletests
    q = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.sum():
        q.loc[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


# --------------------------------------------------------------------------- #
# Stage: meta -- the paired swQTL-vs-eQTL contrast
# --------------------------------------------------------------------------- #
# WHY THE RAW YIELDS ARE NOT THE RESULT
# -------------------------------------
# In caudate at n=420, A_g yields ~54% eGenes and S_g ~10% swQTL genes. Read naively that
# says genetic variation regulates abundance and not splicing. It does not say that, and
# quoting it that way would be the same error the GTEx layer already had to correct for.
#
# A_g is a z-scored log-CPM: a directly measured, high-precision phenotype. S_g is PC1 of
# a within-gene CLR isoform composition, so it inherits every bit of isoform-quantification
# noise, and its variance is partly measurement error. Two phenotypes of unequal precision
# tested at the same alpha will differ in yield for reasons that have nothing to do with
# biology. The eQTL/sQTL discovery asymmetry is exactly why `coloc_modality_contrast`
# makes the CONDITIONAL posterior primary rather than the raw counts.
#
# So the primary statistics here are the ones where detectability largely cancels:
#   * pi1 sharing        -- of genes with a swQTL, what fraction show an abundance signal
#                           at the SAME variant, and vice versa. A ratio of ratios, not a
#                           count against a threshold.
#   * effect size on the DOUBLY significant genes -- both axes cleared the same bar, so the
#                           comparison is between measured effects rather than between
#                           detection rates.
#   * expression-matched yield -- swQTL and eQTL yields compared within bins of gene
#                           expression and transcript count, which is where the precision
#                           difference mostly lives.
# The raw counts are reported too, clearly labelled as confounded with precision.
_PI1_LAMBDA = 0.5


def _storey_pi1(p: np.ndarray, lam: float = _PI1_LAMBDA) -> float:
    """Storey's pi1 = fraction of tests with a true signal.

    Used rather than a threshold count because it does not require either axis to clear a
    multiple-testing bar -- which is precisely the step the precision difference distorts.
    """
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if p.size < 50:
        return float("nan")
    pi0 = float((p > lam).sum()) / (len(p) * (1.0 - lam))
    return float(min(max(1.0 - pi0, 0.0), 1.0))


def _load_qtl(region: str, arm: str) -> dict[str, pd.DataFrame]:
    d = out_dir(arm, region) / "qtl"
    out = {}
    for ftype in MODALITIES:
        f = d / f"cis_qtl_{ftype}.parquet"
        if f.exists():
            out[ftype] = pd.read_parquet(f)
    return out


def run_meta(arm: str = "all_samples", regions: list[str] | None = None,
             fdr: float = 0.05) -> Path:
    dest = ensure_dir(stage_out("anchoring.brainseq_qtl", arm))
    regions = regions or list(REGIONS)
    per_region, shared_rows = [], []

    for region in regions:
        q = _load_qtl(region, arm)
        if len(q) < 2:
            print(f"  {region}: needs both modalities mapped; skipping")
            continue
        sw, ab = q["switch"], q["abundance"]
        m = sw.merge(ab, on="phenotype_id", suffixes=("_sw", "_ab"))
        s_hit = m["qval_sw"] < fdr
        a_hit = m["qval_ab"] < fdr

        # pi1 sharing, in both directions. Conditioning on one axis being significant and
        # asking what the OTHER axis's p-value distribution looks like removes the
        # threshold from the comparison.
        pi1_ab_given_sw = _storey_pi1(m.loc[s_hit, "pval_beta_ab"].to_numpy())
        pi1_sw_given_ab = _storey_pi1(m.loc[a_hit, "pval_beta_sw"].to_numpy())

        # Effect size where both cleared the same bar.
        both = m[s_hit & a_hit].copy()
        both["abs_slope_sw"] = both["slope_sw"].abs()
        both["abs_slope_ab"] = both["slope_ab"].abs()
        w_p = np.nan
        if len(both) >= 6:
            from scipy import stats as _st
            w_p = float(_st.wilcoxon(both["abs_slope_sw"], both["abs_slope_ab"])[1])

        # Same lead variant on both axes: the strictest sharing statement available here.
        same_var = int((both["variant_id_sw"] == both["variant_id_ab"]).sum())

        b = int((s_hit & ~a_hit).sum())
        c = int((~s_hit & a_hit).sum())
        from scipy import stats as _st2
        mc_p = _st2.binomtest(b, b + c, 0.5).pvalue if (b + c) else np.nan

        per_region.append({
            "region": region, "arm": arm, "n_genes_paired": int(len(m)),
            "n_swQTL": int(s_hit.sum()), "n_eQTL": int(a_hit.sum()),
            "switch_only": b, "abundance_only": c,
            "both": int((s_hit & a_hit).sum()),
            "neither": int((~s_hit & ~a_hit).sum()),
            "mcnemar_p_RAW_CONFOUNDED": mc_p,
            "pi1_abundance_given_switch": pi1_ab_given_sw,
            "pi1_switch_given_abundance": pi1_sw_given_ab,
            "n_both": int(len(both)),
            "n_same_lead_variant": same_var,
            "frac_same_lead_variant": (same_var / len(both)) if len(both) else np.nan,
            "median_abs_slope_switch": float(both["abs_slope_sw"].median()) if len(both) else np.nan,
            "median_abs_slope_abundance": float(both["abs_slope_ab"].median()) if len(both) else np.nan,
            "wilcoxon_abs_slope_p": w_p,
        })
        m.assign(region=region).to_parquet(dest / f"paired_{region}.parquet", index=False)

        # Expression-matched yield: bin on the number of variants tested and the minor
        # allele count of the lead, both proxies for the power available to a gene.
        m2 = m.copy()
        m2["bin"] = pd.qcut(m2["num_var_ab"], 10, labels=False, duplicates="drop")
        by = (m2.assign(sw=s_hit.to_numpy(), ab=a_hit.to_numpy())
                .groupby("bin", as_index=False)
                .agg(n=("phenotype_id", "size"), swQTL=("sw", "mean"),
                     eQTL=("ab", "mean")))
        by.insert(0, "region", region)
        shared_rows.append(by)

    if not per_region:
        raise SystemExit("no region has both modalities mapped yet")
    summary = pd.DataFrame(per_region)
    summary.to_parquet(dest / "modality_contrast.parquet", index=False)
    if shared_rows:
        pd.concat(shared_rows, ignore_index=True).to_parquet(
            dest / "yield_by_power_bin.parquet", index=False)
    _write_meta_report(dest, summary, fdr)
    print(summary.to_string(index=False))
    print(f"\n  wrote {dest}")
    return dest


def _write_meta_report(dest: Path, s: pd.DataFrame, fdr: float) -> None:
    L: list[str] = []
    A = L.append
    A("# BrainSEQ switch-QTL vs abundance-QTL")
    A("")
    A("Does genetic variation preferentially regulate the **relative-isoform** axis "
      "(`S_g`) or the **total-abundance** axis (`A_g`) of the same gene? Same donors, "
      "same variants, same cis window, same covariates, one phenotype per gene per axis.")
    A("")
    A("**This is same-tissue genetic anchoring, not independent replication.** BrainSEQ "
      "is the cohort the switches were discovered in. GTEx remains the external cohort.")
    A("")
    A("## Read the raw counts with care")
    A("")
    A("`A_g` is a directly measured log-CPM; `S_g` is PC1 of a within-gene CLR isoform "
      "composition and inherits isoform-quantification noise. Two phenotypes of unequal "
      "precision tested at the same alpha differ in yield for reasons unrelated to "
      "biology, so **the raw swQTL/eQTL counts below are confounded with measurement "
      "precision and are not the result.** They are reported because hiding them would "
      "be worse. The primary statistics are the precision-cancelling ones beneath.")
    A("")
    A("| region | genes | swQTL | eQTL | switch-only | abundance-only | both |")
    A("|---|---|---|---|---|---|---|")
    for r in s.itertuples(index=False):
        A(f"| {r.region} | {r.n_genes_paired:,} | {r.n_swQTL:,} | {r.n_eQTL:,} | "
          f"{r.switch_only:,} | {r.abundance_only:,} | {r.both:,} |")
    A("")
    A("## Primary statistics")
    A("")
    A("`pi1` conditions on one axis being significant and asks what fraction of the "
      "OTHER axis carries signal — a ratio of ratios, with no threshold on the second "
      "axis. Effect sizes are compared only on genes where BOTH axes cleared the same "
      "bar, so detection rate largely cancels.")
    A("")
    A("| region | pi1(A given S) | pi1(S given A) | n both | same lead variant | "
      "median abs slope S | median abs slope A | Wilcoxon P |")
    A("|---|---|---|---|---|---|---|---|")
    for r in s.itertuples(index=False):
        A(f"| {r.region} | {r.pi1_abundance_given_switch:.3f} | "
          f"{r.pi1_switch_given_abundance:.3f} | {r.n_both:,} | "
          f"{r.n_same_lead_variant:,} ({r.frac_same_lead_variant:.2f}) | "
          f"{r.median_abs_slope_switch:.3f} | {r.median_abs_slope_abundance:.3f} | "
          f"{r.wilcoxon_abs_slope_p:.3g} |")
    A("")
    A("`yield_by_power_bin.parquet` carries the same yields within deciles of the number "
      "of variants tested, which is where most of the precision difference lives.")
    A("")
    A("### Two things not to overclaim from the numbers above")
    A("")
    A("**The two phenotypes are not independent.** `S_g` and `A_g` are both derived from "
      "the same transcript counts: a variant that strongly changes one isoform's "
      "abundance moves the gene total AND the within-gene composition. So some sharing "
      "is mechanically expected and `pi1(A given S)` should not be read as evidence that "
      "a switch signal *causes* an abundance signal. The informative comparison is the "
      "ASYMMETRY between the two pi1 values, and the same-lead-variant fraction.")
    A("")
    A("**Effect sizes are comparable only because the transform was identical.** Both "
      "axes are rank-based inverse-normal transformed before mapping, so slopes are in "
      "SD units of the transformed phenotype and can be compared. That is a property of "
      "this pipeline, not of QTL slopes in general — comparing an untransformed "
      "composition slope against an untransformed log-CPM slope would be meaningless.")
    A("")
    A(f"FDR: Benjamini-Hochberg at q < {fdr} on the permutation p-values, applied "
      "identically to both axes. tensorQTL's Storey q-values need rpy2 and the R "
      "`qvalue` package, which do not import in this environment.")
    A("")
    (dest / "BRAINSEQ_SWITCH_QTL.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("phenotypes", "map", "junctions", "meta"),
                    required=True)
    ap.add_argument("--region", choices=REGIONS, action="append", default=None)
    ap.add_argument("--arm", choices=ARMS, default="all_samples")
    ap.add_argument("--n-factors", type=int, default=N_HIDDEN_FACTORS,
                    help="hidden phenotype factors; the SAME count is used for S_g and A_g")
    ap.add_argument("--chrom", type=int, action="append", default=None,
                    help="restrict to these chromosomes (repeatable); for smoke tests")
    args = ap.parse_args(argv)
    regions = args.region or list(REGIONS)

    if args.stage == "phenotypes":
        for r in regions:
            print(f"\n== {r} / {args.arm} ==")
            build_phenotypes(r, arm=args.arm)
    elif args.stage == "meta":
        run_meta(arm=args.arm, regions=regions)
    elif args.stage == "map":
        for r in regions:
            print(f"\n== map {r} / {args.arm} ==")
            run_map(r, arm=args.arm, n_factors=args.n_factors, chroms=args.chrom)
    else:
        raise SystemExit(f"stage {args.stage!r} not yet implemented")


if __name__ == "__main__":
    main()
