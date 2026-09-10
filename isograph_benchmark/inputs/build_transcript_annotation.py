"""Regenerate the BrainSEQ transcript annotation that `build_bundles` depends on.

WHY THIS EXISTS
---------------
`build_bundles.build_brainseq_bundle` reads
`inputs/raw/brainseq/annotations/transcript-annotation.tsv` to attach `gene_id`,
`transcript_name` and `transcript_type` to the RSEM transcript matrix. That file was
never tracked and is no longer on disk, so a clean checkout cannot rebuild the BrainSEQ
bundles at all -- and those bundles feed every IsoGraph fit in the paper. The data was
not lost, only the input that produced it: the committed bundles still carry the columns.

This module reconstructs the file from the GENCODE v47 primary-assembly GTF, which is the
annotation the bundles were built against. It is deliberately a separate, verifiable step
rather than a patch to `build_bundles`:

  * `transcript_name` must be the TRANSCRIPT name (`DRD2-201`), not the gene symbol.
    `run_models._drd2_gene_id` resolves the DRD2 sanity check with
    `transcript_name.str.startswith("DRD2-")`, so substituting gene symbols -- which is
    what the per-region `tx-annotation.tsv` in the shared r-variables tree carries -- would
    not error. It would silently degrade to "WARNING: DRD2 check could not resolve
    gene_id" and the check would stop checking anything.
  * `--verify` compares the regenerated table against a committed bundle's
    `transcripts.parquet` and fails on ANY mismatch. Because this file sits upstream of
    every published fit, "it looks right" is not good enough: the regeneration has to be
    demonstrated to reproduce the annotation the bundles were actually built with, before
    anything downstream is re-run.

Usage:
    python -m isograph_benchmark.inputs.build_transcript_annotation --verify
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel

# The annotation the BrainSEQ bundles were built against. GENCODE v47 primary assembly is
# also what the GTF cache, the switch-consequence layer and the phASER feature BED use --
# the GTF's 387,944 transcripts match the phASER `gene_ae` row count per sample exactly.
DEFAULT_GTF = Path(
    "/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/gtf/"
    "gencode.v47.primary_assembly.annotation.gtf")

# Exactly the columns `build_bundles` selects, plus the gene-level pair, which costs
# nothing and saves the next consumer a second parse.
COLUMNS = ["transcript_id", "gene_id", "transcript_name", "transcript_type",
           "gene_name", "gene_type"]

_ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf(gtf: Path = DEFAULT_GTF) -> pd.DataFrame:
    """One row per transcript, from the GTF's `transcript` feature lines."""
    if not gtf.exists():
        raise SystemExit(f"missing GENCODE GTF: {gtf}")
    rows = []
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "transcript":
                continue
            a = dict(_ATTR.findall(f[8]))
            rows.append([a.get(c) for c in COLUMNS])
    d = pd.DataFrame(rows, columns=COLUMNS)
    missing = [c for c in COLUMNS if d[c].isna().all()]
    if missing:
        raise SystemExit(f"GTF carried no values for {missing}; wrong annotation file?")
    return d


def verify(annot: pd.DataFrame, bundle: Path) -> None:
    """Fail unless the regenerated table reproduces a committed bundle exactly.

    The bundle is the ground truth here: it was written by the original file, so agreement
    proves the regeneration recovers that file's content rather than merely something
    plausible from the same release.
    """
    if not bundle.exists():
        raise SystemExit(f"missing bundle for verification: {bundle}")
    b = pd.read_parquet(bundle)
    need = [c for c in ("transcript_id", "gene_id", "transcript_name", "transcript_type")
            if c in b.columns]
    m = b[need].merge(annot, on="transcript_id", how="left", suffixes=("_bundle", "_new"))

    unresolved = int(m["transcript_name_new"].isna().sum())
    if unresolved:
        raise SystemExit(
            f"{unresolved:,} of {len(b):,} bundle transcripts are absent from the "
            f"regenerated annotation — wrong GENCODE release?")
    bad = []
    for c in ("gene_id", "transcript_name", "transcript_type"):
        if f"{c}_bundle" not in m.columns:
            continue
        diff = m[m[f"{c}_bundle"].astype(str) != m[f"{c}_new"].astype(str)]
        if len(diff):
            bad.append((c, len(diff), diff.head(3)))
    if bad:
        for c, n, ex in bad:
            print(f"  MISMATCH {c}: {n:,} rows differ")
            print(ex.to_string(index=False))
        raise SystemExit("regenerated annotation does NOT reproduce the committed bundle")

    drd2 = int(annot["transcript_name"].str.startswith("DRD2-", na=False).sum())
    if drd2 == 0:
        raise SystemExit("no DRD2- transcript names; the DRD2 sanity check would silently "
                         "stop checking anything")
    print(f"  VERIFIED against {bundle}")
    print(f"    {len(b):,} bundle transcripts, all resolved, "
          f"0 mismatches on {', '.join(c for c in need if c != 'transcript_id')}")
    print(f"    DRD2- transcript names present: {drd2}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--verify", action="store_true",
                    help="check the result reproduces a committed bundle before writing")
    ap.add_argument("--bundle", type=Path,
                    default=rel("inputs", "bundles", "brainseq_v1", "caudate",
                                "transcripts.parquet"))
    args = ap.parse_args(argv)

    out = args.out or rel("inputs", "raw", "brainseq", "annotations",
                          "transcript-annotation.tsv")
    print(f"  parsing {args.gtf}")
    annot = parse_gtf(args.gtf)
    print(f"  {len(annot):,} transcripts, {annot['gene_id'].nunique():,} genes")
    if args.verify:
        verify(annot, args.bundle)
    ensure_dir(out.parent)
    annot.to_csv(out, sep="\t", index=False)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
