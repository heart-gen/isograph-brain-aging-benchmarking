"""Ingest GTEx v11 STAR junction counts into per-region within-gene junction usage.

The raw GTEx junction matrix (``inputs/raw/gtex_v11/counts/..._junctions.gct.gz``,
523,817 junctions x 19,788 samples, all tissues) provides split-read counts per
junction, already annotated to a gene via the GCT ``Description`` column. This module
turns it into a per-brain-region ``junction_usage.parquet`` that is the GTEx analogue
of BrainSEQ's PSI events for :mod:`validate_switch_splicing`:

    usage_{j,s} = reads_{j,s} / sum_{j' in gene(j)} reads_{j',s}

i.e. a within-gene junction-usage fraction. This is derived from split reads and does
NOT use the RSEM transcript quantifier that feeds IsoGraph, so testing usage ~ age is a
genuine orthogonal check of IsoGraph's switch calls.

To keep the output small and relevant, junctions are restricted to genes that IsoGraph
switch-scored in that region (the validation universe). One region per invocation; the
full gzip is streamed once with ``usecols`` limited to the region's bundle samples, so
memory stays modest. Heavy enough to run on SLURM.

Output: ``inputs/processed/gtex_v11/<region>/junction_usage.parquet`` with columns
``gene`` (unversioned), ``junction`` (chr:start-end:strand), and one column per sample.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel

_GCT = rel("inputs", "raw", "gtex_v11", "counts",
           "GTEx_Analysis_2025-08-22_v11_STARv2.7.11b_junctions.gct.gz")

GTEX_REGIONS = [
    "amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
    "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
    "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
    "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra",
]


def _strip_ver(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".").str[0]


def _region_samples(region: str) -> list[str]:
    s = pd.read_parquet(rel("inputs", "bundles", "gtex_v11_brain", region, "samples.parquet"),
                        columns=["sample_id"])
    return list(s["sample_id"].astype(str))


def _switch_scored_genes(region: str) -> set[str]:
    fs = pd.read_parquet(rel("real_data", "gtex", region, "_m", "isograph_vae",
                             "feature_scores.parquet"), columns=["gene_id", "feature_type"])
    return set(_strip_ver(fs.loc[fs["feature_type"] == "switch", "gene_id"]))


def _gct_header_samples() -> list[str]:
    # Row 3 of the GCT is the header: Name, Description, <samples...>.
    hdr = pd.read_csv(_GCT, sep="\t", skiprows=2, nrows=0)
    return list(hdr.columns)


def build_region(region: str, min_total_reads: int = 10, min_samples_expressed: int = 30) -> str:
    samples = _region_samples(region)
    header = _gct_header_samples()
    present = [s for s in samples if s in header]
    missing = len(samples) - len(present)
    genes = _switch_scored_genes(region)
    print(f"[{region}] {len(present)}/{len(samples)} bundle samples in GCT "
          f"({missing} missing); restricting to {len(genes)} switch-scored genes ...", flush=True)

    usecols = ["Name", "Description"] + present
    df = pd.read_csv(_GCT, sep="\t", skiprows=2, usecols=usecols)
    df["gene"] = _strip_ver(df["Description"])
    df = df[df["gene"].isin(genes)].copy()
    print(f"[{region}] read {len(df)} junctions in switch-scored genes; computing usage ...",
          flush=True)

    counts = df[present].to_numpy(dtype=float)
    gene_key = df["gene"].to_numpy()
    # Within-gene, per-sample total across the gene's junctions.
    totals = (pd.DataFrame(counts, index=gene_key).groupby(level=0).transform("sum").to_numpy())
    with np.errstate(invalid="ignore", divide="ignore"):
        usage = np.where(totals > 0, counts / totals, np.nan)

    out = pd.DataFrame(usage, columns=present)
    out.insert(0, "junction", df["Name"].to_numpy())
    out.insert(0, "gene", gene_key)

    # Drop junctions that are essentially never expressed / never quantified.
    expressed = (df[present].to_numpy(dtype=float) >= 1).sum(axis=1) >= min_samples_expressed
    total_reads = df[present].to_numpy(dtype=float).sum(axis=1) >= min_total_reads
    keep = expressed & total_reads
    out = out.loc[keep].reset_index(drop=True)

    dest = ensure_dir(rel("inputs", "processed", "gtex_v11", region)) / "junction_usage.parquet"
    out.to_parquet(dest, index=False, compression="zstd")
    print(f"[{region}] wrote {len(out)} junctions x {len(present)} samples "
          f"({out['gene'].nunique()} genes) -> {dest}", flush=True)
    return str(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--region", action="append", dest="regions",
                        help="GTEx brain region(s); repeatable. Default: all 13.")
    parser.add_argument("--min-total-reads", type=int, default=10)
    parser.add_argument("--min-samples-expressed", type=int, default=30)
    args = parser.parse_args()
    for region in (args.regions or GTEX_REGIONS):
        build_region(region, args.min_total_reads, args.min_samples_expressed)


if __name__ == "__main__":
    main()
