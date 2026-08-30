"""S-LDSC step 1: build sQTL / eQTL SNP annotations for IsoGraph switch genes (hg19).

The "step beyond MAGMA": stratified LD-score regression partitions SCZ SNP-heritability
into functional annotations, robust to LD and gene-size confounds (unlike MAGMA's
size-confounded competitive gene-set test). To make it on-thesis for the isoform-switch
layer, the annotations echo the paper's sQTL/eQTL specificity contrast at the
heritability level:

  sqtl_switch  — GTEx brain sQTL SNPs of IsoGraph switch-module genes
  eqtl_switch  — GTEx brain eQTL SNPs of the SAME genes
  cis_switch   — union background (any brain QTL SNP of these genes)

Fit jointly on top of the baselineLD model, the per-SNP heritability (tau) of
sqtl_switch beyond eqtl_switch + cis_switch + baseline is splicing-specific SCZ
heritability concentrated in the switch layer. SNP overlap between annotations is
handled by the joint regression (Finucane-style), so the sQTL/eQTL contrast is not an
artifact of shared cis SNPs.

Genes: IsoGraph phenotype-associated co-switch module genes (bare Ensembl), matched to
GTEx signif_pairs (eQTL: phenotype_id = gene; sQTL: group_id = gene). GTEx variants are
hg38; the S-LDSC panel is hg19, so variant ids are lifted via the GTEx v8 lookup's
`variant_id_b37` field (no chain file needed). Autosomes only.

Writes 05_genetic_anchoring/_m/ldsc/<analysis>/beds/{sqtl_switch,eqtl_switch,cis_switch}_hg19.bed
plus annot_snps.tsv (variant -> hg19 chr/pos, kind) for provenance.
"""
from __future__ import annotations

import argparse
import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data.coloc_prep import load_switch_genes
from isograph_benchmark.real_data.qtl_anchoring import _GTEX_TISSUE, _bare
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir
from isograph_benchmark.real_data.switch_bundles import get_bundle

V8_LOOKUP = Path("/ocean/projects/bio250020p/shared/resources/public-data/gtex-v8/"
                 "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze."
                 "lookup_table.txt.gz")
DEFAULT_XQTL_DIR = rel("inputs", "raw", "gtex_v11", "xqtl")
_GTEX_BRAIN = sorted(set(_GTEX_TISSUE.values()))
_QVAL = 0.05
# signif_pairs gene-id column differs by kind (eQTL: phenotype_id; sQTL: group_id).
_KIND = {"sQTL": ("sQTLs", "group_id"), "eQTL": ("eQTLs", "phenotype_id")}


def collect_qtl_variants(xqtl_dir: Path, tissues: list[str], genes: set[str]) -> pd.DataFrame:
    """Unique GTEx brain QTL variant ids (b38) per kind, for the given genes."""
    rows = []
    for kind, (suffix, gene_col) in _KIND.items():
        seen: set[str] = set()
        for tissue in tissues:
            path = xqtl_dir / f"{tissue}.v11.{suffix}.signif_pairs.parquet"
            if not path.exists():
                continue
            d = pd.read_parquet(path, columns=[gene_col, "variant_id"])
            d = d[_bare(d[gene_col]).isin(genes)]
            seen.update(d["variant_id"].unique())
        rows.append(pd.DataFrame({"variant_id": sorted(seen), "kind": kind}))
        print(f"  {kind}: {len(seen):,} unique variant ids across {len(tissues)} tissues")
    return pd.concat(rows, ignore_index=True)


def lift_to_hg19(variant_ids: set[str]) -> pd.DataFrame:
    """Map b38 variant_id -> hg19 (chr, pos) via the GTEx v8 lookup's variant_id_b37.

    Streams the (large) lookup once, keeping only the needed variant ids. v8/v11 share
    the b38 variant_id string for common SNPs; v11-only ids drop out (fine for S-LDSC,
    which only annotates reference-panel SNPs anyway).
    """
    tmp = ensure_dir(stage_out("anchoring.ldsc", "_tmp")) / "_want_variants.txt"
    tmp.write_text("\n".join(sorted(variant_ids)) + "\n")
    # cols: 1 variant_id(b38) ... 8 variant_id_b37 (chr_pos_ref_alt_b37)
    awk = (r'NR==FNR{w[$1]=1;next} ($1 in w) && $8!="" {print $1"\t"$8}')
    zcat = subprocess.Popen(["zcat", str(V8_LOOKUP)], stdout=subprocess.PIPE)
    res = subprocess.run(["awk", "-F", "\t", awk, str(tmp), "-"],
                         stdin=zcat.stdout, capture_output=True, text=True)
    zcat.stdout.close(); zcat.wait()
    tmp.unlink(missing_ok=True)
    if res.returncode != 0:
        raise RuntimeError(f"v8 lookup awk failed: {res.stderr[:400]}")
    d = pd.read_csv(StringIO(res.stdout), sep="\t", header=None,
                    names=["variant_id", "b37"])
    parts = d["b37"].str.split("_", expand=True)
    d["chr"] = parts[0]
    d["pos"] = pd.to_numeric(parts[1], errors="coerce")
    d = d[d["chr"].isin([str(i) for i in range(1, 23)]) & d["pos"].notna()]
    d["pos"] = d["pos"].astype(int)
    return d[["variant_id", "chr", "pos"]]


def _write_bed(df: pd.DataFrame, path: Path, label: str) -> None:
    if df.empty:
        print(f"  {label}: 0 SNPs — skipped")
        return
    bed = df[["chr", "pos"]].drop_duplicates().copy()
    bed["start"] = bed["pos"] - 1        # 0-based BED
    bed["end"] = bed["pos"]
    bed["chr"] = pd.Categorical(bed["chr"], [str(i) for i in range(1, 23)], ordered=True)
    bed = bed.sort_values(["chr", "start"])
    bed["chr"] = "chr" + bed["chr"].astype(str)
    bed[["chr", "start", "end"]].to_csv(path, sep="\t", index=False, header=False)
    print(f"  {label}: {len(bed):,} SNP intervals -> {path.name}")


def _collect_switch_genes(pairs: list[tuple[str, str | None]], variant: str,
                          fdr: float, min_recurrence: int = 1) -> set[str]:
    """Phenotype-associated switch-module genes across (analysis, region) pairs.

    With ``min_recurrence > 1`` (pooled bundles), keep only genes whose switch-module
    membership reproduces in at least that many analyses — a size-controlled "core"
    switch layer rather than the broad union, so the S-LDSC annotation stays a specific
    switch-layer signal rather than the generic cis regions of every gene.
    """
    from collections import Counter
    counts: Counter = Counter()
    for analysis, region in pairs:
        iso_dir = _artifact_dir(analysis, region, variant)
        g = set(load_switch_genes(iso_dir, fdr)["gene"])
        counts.update(g)
        if len(pairs) > 1:
            print(f"  {analysis}/{region or '-'}: {len(g)} switch genes")
    genes = {g for g, n in counts.items() if n >= min_recurrence}
    if min_recurrence > 1:
        print(f"  recurrence >= {min_recurrence}: {len(genes)} of {len(counts)} union genes kept")
    return genes


def _build_annotation(genes: set[str], out_dir: Path, xqtl_dir: Path,
                      tissues: list[str]) -> None:
    qv = collect_qtl_variants(xqtl_dir, tissues, genes)
    hg19 = lift_to_hg19(set(qv["variant_id"]))
    qv = qv.merge(hg19, on="variant_id", how="inner")
    print(f"Variants lifted to hg19 (autosomes): {qv['variant_id'].nunique():,}")

    bed_dir = ensure_dir(out_dir / "beds")
    _write_bed(qv[qv["kind"] == "sQTL"], bed_dir / "sqtl_switch_hg19.bed", "sqtl_switch")
    _write_bed(qv[qv["kind"] == "eQTL"], bed_dir / "eqtl_switch_hg19.bed", "eqtl_switch")
    _write_bed(qv, bed_dir / "cis_switch_hg19.bed", "cis_switch (union)")
    qv.to_csv(out_dir / "annot_snps.tsv", sep="\t", index=False)
    print(f"Wrote annotations to {bed_dir}")


def run(analysis: str, region: str | None, variant: str, fdr: float,
        xqtl_dir: Path, tissues: list[str]) -> None:
    genes = _collect_switch_genes([(analysis, region)], variant, fdr)
    print(f"IsoGraph switch genes (phenotype-associated modules): {len(genes)}")
    out_dir = ensure_dir(stage_out("anchoring.ldsc", analysis + (f"_{region}" if region else "")))
    _build_annotation(genes, out_dir, xqtl_dir, tissues)


def run_bundle(bundle: str, variant: str, fdr: float, xqtl_dir: Path,
               tissues: list[str], min_recurrence: int = 1) -> None:
    pairs = get_bundle(bundle)
    print(f"Pooling switch genes across {len(pairs)} '{bundle}' analyses "
          f"(min_recurrence={min_recurrence}):")
    genes = _collect_switch_genes(pairs, variant, fdr, min_recurrence)
    print(f"'{bundle}' core switch genes: {len(genes)}")
    out_dir = ensure_dir(stage_out("anchoring.ldsc", bundle))
    _build_annotation(genes, out_dir, xqtl_dir, tissues)


def main() -> None:
    p = argparse.ArgumentParser(description="S-LDSC step 1: sQTL/eQTL SNP annotations for switch genes.")
    p.add_argument("--analysis", default="brainseq-sczd")
    p.add_argument("--region", default=None)
    p.add_argument("--bundle", default=None,
                   help="pool switch genes across a named bundle (e.g. 'aging'); "
                        "overrides --analysis/--region.")
    p.add_argument("--variant", default="standard")
    p.add_argument("--fdr", type=float, default=_QVAL)
    p.add_argument("--min-recurrence", type=int, default=3,
                   help="for --bundle: keep genes whose switch membership recurs in "
                        ">= this many analyses (size-controlled core switch layer).")
    p.add_argument("--xqtl-dir", default=str(DEFAULT_XQTL_DIR))
    p.add_argument("--tissue", action="append", help="GTEx tissue(s); default all 13 brain.")
    args = p.parse_args()
    tissues = args.tissue or _GTEX_BRAIN
    if args.bundle:
        run_bundle(args.bundle, args.variant, args.fdr, Path(args.xqtl_dir),
                   tissues, args.min_recurrence)
    else:
        run(args.analysis, args.region, args.variant, args.fdr, Path(args.xqtl_dir), tissues)


if __name__ == "__main__":
    main()
