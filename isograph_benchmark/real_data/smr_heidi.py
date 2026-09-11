"""SMR + HEIDI on the signal-level colocalization nominations (Analysis 5).

WHAT THIS ADDS, AND HOW FAR ITS NUMBERS MAY BE READ
---------------------------------------------------
`coloc.susie` asks whether a disease association and a QTL share a causal variant. It gives
no effect direction, and it separates linkage from sharing only through the balance of PP3
and PP4. SMR adds a signed, ratio-type estimate relating the genetically predicted molecular
phenotype to disease (`b_SMR`), and HEIDI asks whether the association pattern across the cis
region is compatible with one shared causal variant rather than distinct variants in LD. That
makes it orthogonal corroboration for loci coloc already nominated -- positioned beneath coloc,
never beside it. Three rules govern every sentence written from this module's output:

  * `b_SMR` is a ratio estimate under SMR's assumptions (one causal variant per probe, no
    horizontal pleiotropy). It does not establish causal direction in the biological sense
    and cannot distinguish causality from pleiotropy.
  * Failing to reject HEIDI is not proof of a shared variant. HEIDI is underpowered at GTEx
    brain sample sizes, so p_HEIDI >= HEIDI_REJECT reads "not rejected", never "shared".
  * A HEIDI rejection does not overrule a strong signal-level colocalization, and a significant
    SMR estimate does not promote a locus coloc did not support. Disagreements are reported as
    disagreements, in `agreement`, and neither method adjudicates the other.

SCOPE
-----
Only the loci the signal-level layer nominated (PP4_sQTL >= 0.8, all-introns arm), only in the
tissues where that sQTL call holds, for every intron phenotype of the gene and for its eQTL.
The colocalizing intron is the pre-specified primary sQTL probe; the gene's other introns are
reported beside it so the SMR evidence is never a maximum chosen over introns after the fact.
This is a scoped follow-up on nominated loci, not a transcriptome-wide scan, and its
multiple-testing family is the probes it actually instrumented.

PRE-SPECIFIED SETTINGS
----------------------
SMR 1.4.2's own defaults, read from its source (src/SMR.cpp) and passed explicitly so a later
release with different defaults cannot change the analysis silently: instrument p < 5e-8, HEIDI
SNPs at p < 1.5654e-3, 3-20 HEIDI SNPs, a 2 Mb cis window, and the allele-frequency QC (0.2 per
SNP, at most 5% failing). Interpretation, fixed before any result existed: SMR significance is
Bonferroni over instrumented probes within (analysis, QTL modality), and p_HEIDI < 0.01 rejects
a single shared causal variant. A relaxed instrument threshold (`--peqtl-smr`) is a sensitivity
arm and writes to its own directory; the ESD/BESD files do not depend on it and are shared.

INPUTS, AND THE BUILD
---------------------
LD reference: the 1000G EUR Phase 3 PLINK panel (hg19) the coloc layer fine-maps on.
GWAS: the per-locus summary statistics the coloc layer used, with the same long-range-LD
exclusions, written as SMR `.ma`. Its frequency column is NA: the per-locus files carry none,
and with a missing GWAS frequency SMR checks the QTL frequency against the reference panel
alone (`freq_check`, src/SMR_data.cpp) rather than against a number invented here.
QTL: GTEx v11 all-pairs nominal statistics. GRCh38 variant ids are bridged to rsIDs and each
ESD position is the panel's hg19 position for that rsID, so GWAS, QTL and LD reference agree on
position. The effect allele is the GTEx ALT allele, which is the allele `af` and `slope` refer
to. Strand-ambiguous and panel-mismatched variants are dropped, as in the coloc layer. Probe
positions are hg19 gene TSSs from MAGMA's NCBI37.3 gene.loc.

BrainSEQ QTL are the plan's second source and are deliberately NOT wired. The cohort is roughly
half African-ancestry while the GWAS and the LD panel are European, so the only defensible input
is an `ea_only` mapping, which has not been run. `--qtl-source brainseq` refuses rather than
produce an SMR estimate on mismatched ancestry.

STAGES
------
  --stage prep   targets from the signal-level nominations, the (tissue, chr) work list, one
                 `.ma` per analysis. Login-node safe.
  --stage besd   per work-list task: GTEx all-pairs -> ESD -> `smr --make-besd`.
  --stage smr    per work-list task: `smr` for every analysis with targets in that cell.
  --stage meta   assemble, apply the multiple-testing family, join coloc per probe, classify,
                 report.

Outputs under 05_genetic_anchoring/_m/smr_heidi/<source>/.
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data import coloc_modality_contrast as cmc
from isograph_benchmark.real_data.coloc_signal_susie import results_dir, signal_root

SMR_BIN = Path("/ocean/projects/bio260021p/shared/opt/SMR/build/Release/smr")
QTL_SOURCES: tuple[str, ...] = ("gtex", "brainseq")
MODALITIES: tuple[str, ...] = ("eQTL", "sQTL")

# SMR 1.4.2 defaults (src/SMR.cpp), passed explicitly on every call.
PEQTL_SMR = 5e-8
PEQTL_HEIDI = 1.5654e-3
HEIDI_MIN_M = 3
HEIDI_MAX_M = 20
CIS_WIND_KB = 2000
DIFF_FREQ = 0.2
DIFF_FREQ_PROP = 0.05

# Interpretation, pre-specified.
SMR_ALPHA = 0.05
HEIDI_REJECT = 0.01
PP4_CALL = cmc.PP4_CALL
P12_PRIMARY = cmc.P12_PRIMARY
MIN_SHARED_SNPS = cmc.MIN_SHARED_SNPS

# File formats, as SMR 1.4.2 parses them: `.ma` in read_gwas_data (8 columns, header starting
# SNP), the ESD in get_esi_info (9 columns, header starting Chr; beta/se/p are columns 7-9),
# the flist in read_probeinfolst (7 columns, header starting Chr).
MA_COLS = ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]
ESD_COLS = ["Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p"]
FLIST_COLS = ["Chr", "ProbeID", "GeneticDistance", "ProbeBp", "Gene", "Orientation", "PathOfEsd"]
SMR_OUT_COLS = ["probeID", "ProbeChr", "Gene", "Probe_bp", "topSNP", "topSNP_chr", "topSNP_bp",
                "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL",
                "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]

AGREEMENT: tuple[str, ...] = (
    "coloc_and_smr_heidi_not_rejected",   # coloc call, SMR significant, HEIDI not rejected
    "coloc_and_smr_heidi_rejected",       # coloc call, SMR significant, HEIDI rejects
    "coloc_and_smr_heidi_untestable",     # coloc call, SMR significant, < HEIDI_MIN_M SNPs
    "coloc_smr_not_significant",          # coloc call, SMR instrumented but not significant
    "coloc_no_instrument",                # coloc call, no cis-QTL reaches the SMR threshold
    "smr_without_coloc",                  # SMR significant where coloc scored below the call
    "smr_coloc_not_scored",               # SMR significant at a probe coloc never scored
    "neither",
)

_AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
_LOCUS_CHR = re.compile(r"_chr(\d+)$")


# --------------------------------------------------------------------------- #
# Layout and the source guard
# --------------------------------------------------------------------------- #
def check_qtl_source(qtl_source: str, arm: str | None = None) -> None:
    """Refuse any QTL source this module cannot use defensibly."""
    if qtl_source not in QTL_SOURCES:
        raise SystemExit(f"unknown --qtl-source {qtl_source!r}; choose from {QTL_SOURCES}")
    if qtl_source != "brainseq":
        return
    if arm != "ea_only":
        raise SystemExit(
            f"BrainSEQ SMR must use the `ea_only` arm, not {arm!r}: the cohort is roughly half "
            "African-ancestry while the GWAS and the 1000G EUR LD reference are European, so "
            "an SMR estimate or a HEIDI test on the mixed-ancestry QTL would rest on mismatched LD.")
    qdir = stage_out("anchoring.brainseq_qtl") / "ea_only"
    raise SystemExit(
        "BrainSEQ SMR is not wired yet. It needs the ea_only QTL mapping under "
        f"{qdir} (`brainseq_switch_qtl --stage map --arm ea_only`), which has not been run.")


def source_root(qtl_source: str = "gtex") -> Path:
    """Targets, GWAS `.ma` and BESD for one QTL source (independent of the SMR threshold)."""
    return stage_out("anchoring.smr") / qtl_source


def run_root(qtl_source: str = "gtex", peqtl_smr: float = PEQTL_SMR) -> Path:
    """SMR runs and meta outputs; a non-default instrument threshold is a sensitivity arm."""
    base = source_root(qtl_source)
    return base if np.isclose(peqtl_smr, PEQTL_SMR) else (
        base / "sensitivity" / f"peqtl_smr_{peqtl_smr:g}")


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))


def locus_chr(locus_id: str) -> int:
    m = _LOCUS_CHR.search(str(locus_id))
    if not m:
        raise SystemExit(f"cannot read a chromosome from LOCUS_ID {locus_id!r}")
    return int(m.group(1))


def probe_gene(probe_id: str) -> str:
    """Bare ENSG of a GTEx probe: `ENSG...v` (eQTL) or `chr:start:end:clu_N_s:ENSG...v` (sQTL)."""
    return str(probe_id).split(":")[-1].split(".")[0]


# --------------------------------------------------------------------------- #
# Stage: prep
# --------------------------------------------------------------------------- #
def build_targets(nom: pd.DataFrame, cells: pd.DataFrame, call: float = PP4_CALL) -> pd.DataFrame:
    """(nominated locus, gene, tissue) cells where the sQTL call holds."""
    keys = ["analysis", "trait", "LOCUS_ID", "gene"]
    c = cells.merge(nom[keys].drop_duplicates(), on=keys, how="inner")
    c = c[c["PP4_sQTL"] >= call].copy()
    c["chr"] = c["LOCUS_ID"].map(locus_chr)
    c = c.rename(columns={"estimator_sQTL": "coloc_estimator_sQTL",
                          "PP4_sQTL": "coloc_PP4_sQTL", "PP4_eQTL": "coloc_PP4_eQTL",
                          "phenotype_id": "coloc_phenotype_id"})
    keep = [*keys, "chr", "symbol", "tissue", "coloc_estimator_sQTL", "coloc_PP4_sQTL",
            "coloc_PP4_eQTL", "coloc_phenotype_id"]
    for k in keep:
        if k not in c.columns:
            c[k] = None
    return c[keep].sort_values(["analysis", "LOCUS_ID", "gene", "tissue"]).reset_index(drop=True)


def ma_from_loci(loci: list[tuple[int, pd.DataFrame]], excl: pd.DataFrame, n: int) -> pd.DataFrame:
    """Per-locus GWAS (rsid, pos, a1, a2, beta, se, p) as SMR `.ma` rows.

    `a1` is the effect allele in these files (stage 22 re-signs z by it). The long-range-LD
    exclusions are the coloc layer's, applied the same way, so both methods see one GWAS.
    """
    parts = []
    for ch, g in loci:
        g = g.copy()
        for r in excl[excl["chr"] == ch].itertuples(index=False):
            g = g[~((g["pos"] >= r.start) & (g["pos"] <= r.stop))]
        parts.append(g)
    if not parts:
        return pd.DataFrame(columns=MA_COLS)
    g = pd.concat(parts, ignore_index=True)
    g = g[np.isfinite(g["beta"]) & np.isfinite(g["se"]) & (g["se"] > 0)]
    g = g.drop_duplicates("rsid")
    return pd.DataFrame({"SNP": g["rsid"].values, "A1": g["a1"].str.upper().values,
                         "A2": g["a2"].str.upper().values, "freq": np.nan,
                         "b": g["beta"].values, "se": g["se"].values, "p": g["p"].values,
                         "n": int(n)})[MA_COLS]


def write_ma(ma: pd.DataFrame, path: Path) -> None:
    ma.to_csv(path, sep="\t", index=False, na_rep="NA")


def run_prep(qtl_source: str = "gtex", arm: str | None = None, call: float = PP4_CALL) -> Path:
    check_qtl_source(qtl_source, arm)
    from isograph_benchmark.real_data.locus_event_audit import load_nominations

    nom, cells = load_nominations("susie", call=call, sqtl_arm="all")
    t = build_targets(nom, cells, call=call)
    if t.empty:
        raise SystemExit("no nominated (locus, gene, tissue) cells; nothing to test")
    dest = ensure_dir(source_root(qtl_source))
    t.to_parquet(dest / "targets.parquet", index=False)
    work = t[["tissue", "chr"]].drop_duplicates().sort_values(["tissue", "chr"])
    work.to_csv(dest / "work_list.tsv", sep="\t", index=False)

    gm = pd.read_csv(cmc.out_dir("switch") / "gwas_meta.tsv", sep="\t").set_index("trait")
    gdir = ensure_dir(dest / "gwas")
    for analysis, sub in t.groupby("analysis"):
        cdir = cmc.coloc_dir(analysis)
        ef = cdir / "exclude_regions.tsv"
        excl = pd.read_csv(ef, sep="\t") if ef.exists() else pd.DataFrame(
            columns=["chr", "start", "stop"])
        loci = [(locus_chr(lid), pd.read_csv(cdir / "susie" / f"{lid}.gwas.tsv", sep="\t"))
                for lid in sorted(set(sub["LOCUS_ID"]))]
        ma = ma_from_loci(loci, excl, n=int(gm.at[sub["trait"].iloc[0], "n_gwas"]))
        write_ma(ma, gdir / f"{analysis}.ma")
        print(f"  {analysis}: {len(ma):,} GWAS SNPs over {len(loci)} loci")
    print(f"  {len(t)} (locus, gene, tissue) targets, {t['gene'].nunique()} genes, "
          f"{len(work)} (tissue, chr) tasks -> {dest}")
    return dest


# --------------------------------------------------------------------------- #
# Stage: besd
# --------------------------------------------------------------------------- #
def esd_rows(q: pd.DataFrame, bridge: pd.DataFrame, bim: pd.DataFrame) -> pd.DataFrame:
    """GTEx nominal statistics as SMR ESD rows, keyed by `phenotype_id`.

    `q`: phenotype_id, variant_id, af, pval_nominal, slope, slope_se. `bridge`: variant_id ->
    rsid. `bim`: the reference panel's bchr, rsid, bp, A1, A2 (hg19).
    """
    cols = ["phenotype_id", *ESD_COLS]
    if q.empty:
        return pd.DataFrame(columns=cols)
    d = (q.merge(bridge[["variant_id", "rsid"]], on="variant_id")
          .merge(bim[["rsid", "bchr", "bp", "A1", "A2"]], on="rsid"))
    if d.empty:
        return pd.DataFrame(columns=cols)
    parts = d["variant_id"].str.split("_", expand=True)
    ref, alt = parts[2].str.upper(), parts[3].str.upper()
    p1, p2 = d["A1"].astype(str).str.upper(), d["A2"].astype(str).str.upper()
    same = ((ref == p1) & (alt == p2)) | ((ref == p2) & (alt == p1))
    ambiguous = pd.Series([(a, b) in _AMBIGUOUS for a, b in zip(ref, alt)], index=d.index)
    # SMR refuses a frequency of exactly 0 or 1, and a zero SE has no ratio estimate.
    ok = same & ~ambiguous & (d["af"] > 0) & (d["af"] < 1) & (d["slope_se"] > 0)
    d, ref, alt = d[ok], ref[ok], alt[ok]
    out = pd.DataFrame({"phenotype_id": d["phenotype_id"].values,
                        "Chr": d["bchr"].astype(int).values, "SNP": d["rsid"].values,
                        "Bp": d["bp"].astype(int).values, "A1": alt.values, "A2": ref.values,
                        "Freq": d["af"].values, "Beta": d["slope"].values,
                        "se": d["slope_se"].values, "p": d["pval_nominal"].values})
    return out.drop_duplicates(["phenotype_id", "SNP"]).reset_index(drop=True)


def gene_tss_hg19(genes: pd.DataFrame) -> pd.DataFrame:
    """(gene, symbol) -> hg19 TSS and strand from MAGMA's NCBI37.3 gene.loc."""
    from isograph_benchmark.real_data.coloc_prep import GENE_LOC_HG19

    loc = pd.read_csv(GENE_LOC_HG19, sep=r"\s+", header=None,
                      names=["entrez", "chr", "start", "stop", "strand", "symbol"],
                      dtype={"chr": str})
    loc = loc.drop_duplicates("symbol")
    loc["tss"] = np.where(loc["strand"] == "-", loc["stop"], loc["start"])
    return (genes[["gene", "symbol"]].drop_duplicates("gene")
            .merge(loc[["symbol", "strand", "tss"]], on="symbol", how="left"))


def flist_rows(esd: pd.DataFrame, genes: pd.DataFrame, esd_dir: Path) -> pd.DataFrame:
    """One `--eqtl-flist` row per probe.

    A probe whose gene has no hg19 TSS falls back to the median ESD position with orientation
    NA. SMR uses the probe position only to centre its 2 Mb cis window, and GTEx's cis
    variants (within 1 Mb of the TSS) already sit inside it.
    """
    g = genes.drop_duplicates("gene").set_index("gene")
    rows = []
    for pid, blk in esd.groupby("phenotype_id", sort=True):
        gene = probe_gene(pid)
        has = gene in g.index
        sym = str(g.at[gene, "symbol"]) if has and pd.notna(g.at[gene, "symbol"]) else gene
        if has and "tss" in g.columns and pd.notna(g.at[gene, "tss"]):
            bp, orient = int(g.at[gene, "tss"]), str(g.at[gene, "strand"])
        else:
            bp, orient = int(blk["Bp"].median()), "NA"
        rows.append({"Chr": int(blk["Chr"].iloc[0]), "ProbeID": pid, "GeneticDistance": 0,
                     "ProbeBp": bp, "Gene": sym, "Orientation": orient,
                     "PathOfEsd": str(esd_dir / f"{_safe(pid)}.esd")})
    return pd.DataFrame(rows, columns=FLIST_COLS)


def read_gtex_qtl(tissue: str, chrom: int, genes: set[str], modality: str) -> pd.DataFrame:
    """Nominal cis statistics for the target genes, by parquet predicate pushdown."""
    import pyarrow.compute as pc
    import pyarrow.dataset as pds

    if modality == "eQTL":
        f = (cmc.GTEX_V11 / "GTEx_Analysis_v11_eQTL_all_associations"
             / f"{tissue}.v11.allpairs.chr{chrom}.parquet")
        col = "gene_id"
    else:
        f = (cmc.GTEX_V11 / "GTEx_Analysis_v11_sQTL_all_associations"
             / f"{tissue}.v11.cis_sqtl.allpairs.chr{chrom}.parquet")
        col = "phenotype_id"
    cols = ["phenotype_id", "variant_id", "af", "pval_nominal", "slope", "slope_se"]
    if not f.exists() or not genes:
        return pd.DataFrame(columns=cols)
    gs = sorted(genes)
    expr = pc.match_substring(pds.field(col), gs[0])
    for g in gs[1:]:
        expr = expr | pc.match_substring(pds.field(col), g)
    t = (pds.dataset(str(f))
         .to_table(columns=[col, "variant_id", "af", "pval_nominal", "slope", "slope_se"],
                   filter=expr)
         .to_pandas().rename(columns={col: "phenotype_id"}))
    # A substring hit is not an identity: keep only probes whose gene IS a target.
    return t[t["phenotype_id"].map(probe_gene).isin(genes)][cols]


def _read_bim(chrom: int) -> pd.DataFrame:
    from isograph_benchmark.real_data.coloc_prep import PANEL_DIR

    return pd.read_csv(PANEL_DIR / f"1000G.EUR.QC.{chrom}.bim", sep=r"\s+", header=None,
                       names=["bchr", "rsid", "cm", "bp", "A1", "A2"],
                       dtype={"A1": str, "A2": str})


def _run(cmd: list[str], log: Path) -> int:
    with open(log, "w") as fh:
        fh.write(" ".join(cmd) + "\n\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode


def _tasks(root: Path, task: int | None) -> pd.DataFrame:
    wf = root / "work_list.tsv"
    if not wf.exists():
        raise SystemExit(f"missing {wf}; run --stage prep first")
    work = pd.read_csv(wf, sep="\t")
    if task is None:
        return work
    if not 0 <= task < len(work):
        raise SystemExit(f"--task {task} out of range (0-{len(work) - 1})")
    return work.iloc[[task]]


def run_besd(qtl_source: str = "gtex", task: int | None = None, arm: str | None = None) -> None:
    check_qtl_source(qtl_source, arm)
    root = source_root(qtl_source)
    targets = pd.read_parquet(root / "targets.parquet")
    for w in _tasks(root, task).itertuples(index=False):
        tissue, ch = str(w.tissue), int(w.chr)
        sub = targets[(targets["tissue"] == tissue) & (targets["chr"] == ch)]
        genes = set(sub["gene"])
        bridge = pd.read_parquet(cmc.BRIDGE_DIR / f"chr{ch}.parquet")
        bim = _read_bim(ch)
        tss = gene_tss_hg19(sub[["gene", "symbol"]])
        for mod in MODALITIES:
            esd = esd_rows(read_gtex_qtl(tissue, ch, genes, mod), bridge, bim)
            if esd.empty:
                print(f"  {tissue} chr{ch} {mod}: no usable QTL rows")
                continue
            prefix = ensure_dir(root / "besd" / tissue) / f"{mod}.chr{ch}"
            esd_dir = ensure_dir(Path(f"{prefix}.esd"))
            for pid, blk in esd.groupby("phenotype_id"):
                blk[ESD_COLS].to_csv(esd_dir / f"{_safe(pid)}.esd", sep="\t", index=False)
            flist = Path(f"{prefix}.flist")
            flist_rows(esd, tss, esd_dir).to_csv(flist, sep="\t", index=False)
            rc = _run([str(SMR_BIN), "--eqtl-flist", str(flist), "--make-besd",
                       "--out", str(prefix)], Path(f"{prefix}.make_besd.log"))
            print(f"  {tissue} chr{ch} {mod}: {esd['phenotype_id'].nunique()} probes, "
                  f"{len(esd):,} rows -> {prefix}.besd (rc={rc})")
            if rc != 0:
                raise SystemExit(f"smr --make-besd failed (rc={rc}); see {prefix}.make_besd.log")


# --------------------------------------------------------------------------- #
# Stage: smr
# --------------------------------------------------------------------------- #
def smr_command(bfile: Path, ma: Path, besd: Path, probes: Path, out: Path,
                peqtl_smr: float = PEQTL_SMR, threads: int = 4) -> list[str]:
    return [str(SMR_BIN), "--bfile", str(bfile), "--gwas-summary", str(ma),
            "--beqtl-summary", str(besd), "--extract-probe", str(probes), "--out", str(out),
            "--peqtl-smr", f"{peqtl_smr:g}", "--peqtl-heidi", f"{PEQTL_HEIDI:g}",
            "--heidi-min-m", str(HEIDI_MIN_M), "--heidi-max-m", str(HEIDI_MAX_M),
            "--cis-wind", str(CIS_WIND_KB), "--diff-freq", f"{DIFF_FREQ:g}",
            "--diff-freq-prop", f"{DIFF_FREQ_PROP:g}", "--thread-num", str(threads)]


def run_smr(qtl_source: str = "gtex", task: int | None = None, peqtl_smr: float = PEQTL_SMR,
            threads: int = 4, arm: str | None = None) -> None:
    check_qtl_source(qtl_source, arm)
    from isograph_benchmark.real_data.coloc_prep import PANEL_DIR

    root, dest = source_root(qtl_source), run_root(qtl_source, peqtl_smr)
    targets = pd.read_parquet(root / "targets.parquet")
    for w in _tasks(root, task).itertuples(index=False):
        tissue, ch = str(w.tissue), int(w.chr)
        cell = targets[(targets["tissue"] == tissue) & (targets["chr"] == ch)]
        for analysis, sub in cell.groupby("analysis"):
            ma = root / "gwas" / f"{analysis}.ma"
            for mod in MODALITIES:
                besd = root / "besd" / tissue / f"{mod}.chr{ch}"
                if not Path(f"{besd}.besd").exists():
                    continue
                epi = pd.read_csv(f"{besd}.epi", sep="\t", header=None,
                                  names=["chr", "probeID", "cm", "bp", "gene", "orientation"])
                probes = epi.loc[epi["probeID"].map(probe_gene).isin(set(sub["gene"])),
                                 "probeID"]
                if probes.empty:
                    continue
                odir = ensure_dir(dest / "smr" / analysis / tissue)
                pf = odir / f"{mod}.chr{ch}.probes"
                probes.to_csv(pf, index=False, header=False)
                out = odir / f"{mod}.chr{ch}"
                cmd = smr_command(PANEL_DIR / f"1000G.EUR.QC.{ch}", ma, besd, pf, out,
                                  peqtl_smr=peqtl_smr, threads=threads)
                rc = _run(cmd, Path(f"{out}.log"))
                # SMR exits non-zero when no probe clears the instrument threshold; that is
                # an outcome (`coloc_no_instrument`), so it is recorded, not raised.
                print(f"  {analysis} {tissue} chr{ch} {mod}: {len(probes)} probes (rc={rc})")


# --------------------------------------------------------------------------- #
# Stage: meta
# --------------------------------------------------------------------------- #
def read_smr(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t")
    missing = [c for c in SMR_OUT_COLS if c not in d.columns]
    if missing:
        raise SystemExit(f"{path} lacks SMR output columns {missing}")
    return d


def smr_threshold(n_instrumented: int, alpha: float = SMR_ALPHA) -> float:
    return alpha / n_instrumented if n_instrumented > 0 else np.nan


def classify(pp4, p_smr, p_heidi, nsnp_heidi, threshold, call: float = PP4_CALL) -> str:
    coloc = pd.notna(pp4) and pp4 >= call
    if pd.isna(p_smr):
        return "coloc_no_instrument" if coloc else "neither"
    sig = pd.notna(threshold) and p_smr <= threshold
    if not coloc:
        # An intron coloc never scored (no GTEx-matched QTL credible set) is not a coloc
        # negative; keep it apart from one that was scored and fell short.
        if not sig:
            return "neither"
        return "smr_without_coloc" if pd.notna(pp4) else "smr_coloc_not_scored"
    if not sig:
        return "coloc_smr_not_significant"
    if pd.isna(p_heidi) or pd.isna(nsnp_heidi) or nsnp_heidi < HEIDI_MIN_M:
        return "coloc_and_smr_heidi_untestable"
    return ("coloc_and_smr_heidi_rejected" if p_heidi < HEIDI_REJECT
            else "coloc_and_smr_heidi_not_rejected")


def probe_coloc(pairs: pd.DataFrame, hierarchy: pd.DataFrame, srep: pd.DataFrame) -> pd.DataFrame:
    """The coloc posterior for each SMR probe, at the resolution the probe has.

    SMR tests one phenotype per probe, so the comparison is per intron, not per cell: a
    signal-level PP4 for that intron (GTEx-matched pairs, primary prior, shared-SNP floor)
    where SuSiE scored it, else the abf posterior -- which exists only for GTEx's
    representative intron, so no other intron of an abf cell inherits it.
    """
    k = ["analysis", "trait", "LOCUS_ID", "tissue", "modality"]
    sp = pairs[(pairs["p12"] == P12_PRIMARY)
               & pairs["cs_matches_gtex"].fillna(False).astype(bool)
               & (pairs["n_shared"] >= MIN_SHARED_SNPS)].copy()
    sp["probe_key"] = np.where(sp["modality"] == "eQTL", sp["gene"].astype(str),
                               sp["phenotype_id"].astype(str))
    s = (sp.groupby([*k, "probe_key"], as_index=False)["PP4"].max()
           .assign(coloc_estimator_probe="susie"))
    ha = hierarchy[hierarchy["estimator"] == "abf"]
    e = ha[ha["modality"] == "eQTL"].assign(probe_key=lambda x: x["gene"].astype(str))
    q = ha[ha["modality"] == "sQTL"].merge(
        srep[["tissue", "gene", "phenotype_id"]].rename(columns={"phenotype_id": "probe_key"}),
        on=["tissue", "gene"], how="inner")
    a = pd.concat([e, q], ignore_index=True)
    a = a[[*k, "probe_key", "PP4"]].assign(coloc_estimator_probe="abf")
    out = pd.concat([s, a], ignore_index=True).rename(columns={"PP4": "coloc_PP4_probe"})
    return out.drop_duplicates([*k, "probe_key"], keep="first").reset_index(drop=True)


_PROBE_FILE = re.compile(r"^(eQTL|sQTL)\.chr(\d+)\.probes$")


def run_meta(qtl_source: str = "gtex", peqtl_smr: float = PEQTL_SMR,
             arm: str | None = None) -> Path:
    check_qtl_source(qtl_source, arm)
    root, dest = source_root(qtl_source), run_root(qtl_source, peqtl_smr)
    targets = pd.read_parquet(root / "targets.parquet")

    tested, results, runs = [], [], []
    for pf in sorted((dest / "smr").rglob("*.probes")):
        m = _PROBE_FILE.match(pf.name)
        if not m:
            continue
        analysis, tissue, mod = pf.parent.parent.name, pf.parent.name, m.group(1)
        pr = pd.read_csv(pf, header=None, names=["probeID"]).assign(
            analysis=analysis, tissue=tissue, modality=mod)
        tested.append(pr)
        sf = Path(str(pf)[:-len(".probes")] + ".smr")
        runs.append({"analysis": analysis, "tissue": tissue, "modality": mod,
                     "chr": int(m.group(2)), "n_probes": len(pr), "smr_file": sf.exists()})
        if sf.exists():
            s = read_smr(sf)
            if len(s):  # a header-only .smr is a run where no probe was instrumented
                results.append(s.assign(analysis=analysis, tissue=tissue, modality=mod))
    if not tested:
        raise SystemExit(f"no SMR runs under {dest / 'smr'}; run --stage smr first")

    d = pd.concat(tested, ignore_index=True)
    stat = ["topSNP", "A1", "A2", "b_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
            "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
    if results:
        r = pd.concat(results, ignore_index=True)[["analysis", "tissue", "modality", "probeID",
                                                   *stat]]
        d = d.merge(r, on=["analysis", "tissue", "modality", "probeID"], how="left")
    else:
        for c in stat:
            d[c] = np.nan
    d["gene"] = d["probeID"].map(probe_gene)
    d = d.merge(targets[["analysis", "trait", "LOCUS_ID", "gene", "symbol", "tissue",
                         "coloc_phenotype_id"]],
                on=["analysis", "gene", "tissue"], how="inner")
    d["probe_key"] = np.where(d["modality"] == "eQTL", d["gene"], d["probeID"])

    sig = results_dir(signal_root("switch"), "all")
    pairs = pd.read_parquet(sig / "signal_pairs.parquet")
    hier = pd.read_parquet(sig / "cells_hierarchy.parquet")
    srep = pd.read_parquet(cmc.out_dir("switch") / "sqtl_representative.parquet")
    d = d.merge(probe_coloc(pairs, hier, srep),
                on=["analysis", "trait", "LOCUS_ID", "tissue", "modality", "probe_key"],
                how="left")
    d["primary_probe"] = (d["modality"] == "eQTL") | (
        (d["probeID"] == d["coloc_phenotype_id"])
        | (d["coloc_phenotype_id"].isna() & (d["coloc_estimator_probe"] == "abf")))

    # Multiple-testing family: instrumented (tissue, probe) within (analysis, modality).
    fam = (d[d["p_SMR"].notna()].drop_duplicates(["analysis", "modality", "tissue", "probeID"])
           .groupby(["analysis", "modality"]).size().rename("n_instrumented").reset_index())
    d = d.merge(fam, on=["analysis", "modality"], how="left")
    d["n_instrumented"] = d["n_instrumented"].fillna(0).astype(int)
    d["smr_threshold"] = d["n_instrumented"].map(smr_threshold)
    d["agreement"] = [classify(r.coloc_PP4_probe, r.p_SMR, r.p_HEIDI, r.nsnp_HEIDI,
                               r.smr_threshold)
                      for r in d.itertuples(index=False)]

    d.to_parquet(dest / "smr_results.parquet", index=False)
    pd.DataFrame(runs).to_csv(dest / "smr_runs.tsv", sep="\t", index=False)
    _write_report(dest, d, qtl_source, peqtl_smr)
    print(f"  {len(d):,} probe rows; agreement: {d['agreement'].value_counts().to_dict()}")
    print(f"  wrote {dest}")
    return dest


def _fmt(v, nd=3) -> str:
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    if isinstance(v, (float, np.floating)):
        return f"{v:.{nd}g}" if (abs(v) < 1e-3 and v != 0) else f"{v:.{nd}f}"
    return str(v)


def _write_report(dest: Path, d: pd.DataFrame, qtl_source: str, peqtl_smr: float) -> None:
    L: list[str] = []
    A = L.append
    A("# SMR + HEIDI on the signal-level colocalization nominations")
    A("")
    A(f"QTL source: `{qtl_source}`. Instrument threshold `--peqtl-smr {peqtl_smr:g}`"
      + ("" if np.isclose(peqtl_smr, PEQTL_SMR) else " -- **SENSITIVITY ARM**") + ". "
      f"HEIDI SNPs at p < {PEQTL_HEIDI:g}, {HEIDI_MIN_M}-{HEIDI_MAX_M} SNPs, "
      f"{CIS_WIND_KB} kb cis window. Significance: Bonferroni at {SMR_ALPHA} over instrumented "
      f"probes per (analysis, modality). HEIDI rejects a single shared variant at p < "
      f"{HEIDI_REJECT}.")
    A("")
    A("## How these numbers may be read")
    A("")
    A("- `b_SMR` is a signed ratio estimate relating genetically predicted phenotype to disease "
      "under a single-causal-variant, no-pleiotropy model. It does **not** establish causal "
      "direction and cannot distinguish causality from horizontal pleiotropy.")
    A("- An sQTL probe is a LeafCutter intron-excision ratio, which is compositional within its "
      "cluster: introns sharing a splice site trade usage, so sibling probes carry opposite "
      "`b_SMR` signs by construction. A sign is read relative to its cluster, never alone.")
    A("- Failing to reject HEIDI is **not** evidence of a shared variant; HEIDI is underpowered at "
      "GTEx brain sample sizes. It reads \"not rejected\".")
    A("- A HEIDI rejection does **not** overrule a strong signal-level colocalization, and SMR "
      "significance does not promote a locus coloc did not support. Disagreements stay "
      "disagreements.")
    A("")
    A("## Agreement with coloc, primary probes")
    A("")
    prim = d[d["primary_probe"]]
    A("| analysis | modality | probes | instrumented | threshold | "
      + " | ".join(f"`{a}`" for a in AGREEMENT) + " |")
    A("|---|---|---|---|---|" + "---|" * len(AGREEMENT))
    for (an, mod), sub in prim.groupby(["analysis", "modality"]):
        cnt = sub["agreement"].value_counts()
        A(f"| {an} | {mod} | {len(sub)} | {int(sub['p_SMR'].notna().sum())} | "
          f"{_fmt(sub['smr_threshold'].iloc[0])} | "
          + " | ".join(str(int(cnt.get(a, 0))) for a in AGREEMENT) + " |")
    A("")
    A("## Primary probes")
    A("")
    A("| gene | trait | tissue | modality | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | "
      "agreement |")
    A("|---|---|---|---|---|---|---|---|---|")
    for r in prim.sort_values(["trait", "symbol", "tissue", "modality"]).itertuples(index=False):
        A(f"| {r.symbol} | {r.trait} | {r.tissue} | {r.modality} | "
          f"{_fmt(r.coloc_PP4_probe)} ({r.coloc_estimator_probe if pd.notna(r.coloc_estimator_probe) else '—'}) | "
          f"{_fmt(r.b_SMR)} ({_fmt(r.se_SMR)}) | {_fmt(r.p_SMR)} | "
          f"{_fmt(r.p_HEIDI)} ({_fmt(r.nsnp_HEIDI, 0)}) | `{r.agreement}` |")
    A("")
    A("Non-primary sQTL probes (the gene's other introns) are in `smr_results.parquet` with "
      "`primary_probe = False`; they are reported so the SMR evidence is not a maximum over "
      "introns, and they are not summarized here.")
    A("")
    (dest / "SMR_HEIDI.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("prep", "besd", "smr", "meta"), required=True)
    ap.add_argument("--qtl-source", choices=QTL_SOURCES, default="gtex")
    ap.add_argument("--arm", default=None, help="BrainSEQ genotype arm (must be ea_only)")
    ap.add_argument("--task", type=int, default=None,
                    help="0-based row of work_list.tsv (besd/smr); default: every row")
    ap.add_argument("--peqtl-smr", type=float, default=PEQTL_SMR,
                    help="instrument threshold; anything but 5e-8 is a sensitivity arm")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args(argv)
    if args.stage == "prep":
        run_prep(args.qtl_source, arm=args.arm)
    elif args.stage == "besd":
        run_besd(args.qtl_source, task=args.task, arm=args.arm)
    elif args.stage == "smr":
        run_smr(args.qtl_source, task=args.task, peqtl_smr=args.peqtl_smr,
                threads=args.threads, arm=args.arm)
    else:
        run_meta(args.qtl_source, peqtl_smr=args.peqtl_smr, arm=args.arm)


if __name__ == "__main__":
    main()
