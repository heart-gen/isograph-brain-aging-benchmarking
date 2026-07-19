"""Harmonize the IFGC (Ferrari et al. 2014) FTD discovery sumstats for the pipeline.

The raw IFGC file (`SUM-STATS-Discovery-full.txt`) is space-delimited with columns
`CHR BP SNP PVALUE A1 OR SE T`, where SNP is `chr:pos:ref:alt` (b37) with NO rsID and
the effect is an odds ratio. Every downstream tool here is rsID-keyed (MAGMA matches the
g1000_eur panel; S-LDSC/coloc match HapMap3 / the LD panel), so this maps each SNP to the
g1000_eur panel rsID by (chr, pos) with allele-set agreement, converts OR->beta=log(OR),
and writes a clean, tab-delimited harmonized file:

  SNP  CHR  BP  A1  A2  BETA  SE  PVALUE

(A1 = effect allele from the raw file; A2 = the other allele parsed from the SNP id.)
Only SNPs present in the panel are kept (the ones MAGMA / the LD panel can use). N is not
in the sumstats and is supplied downstream as a fixed value (IFGC discovery sample size).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("/ocean/projects/bio260021p/shared/resources/gwas/ftd/SUM-STATS-Discovery-full.txt")
BIM = Path("/ocean/projects/bio250020p/shared/opt/magma-v1.10/g1000_eur.bim")
OUT = Path("/ocean/projects/bio260021p/shared/resources/gwas/ftd/ftd_ifgc2014_rsid.tsv.gz")


def _load_panel() -> pd.DataFrame:
    bim = pd.read_csv(BIM, sep="\t", header=None, usecols=[0, 1, 3, 4, 5],
                      names=["chr", "rsid", "pos", "pA1", "pA2"],
                      dtype={"chr": "int32", "pos": "int64", "rsid": "string",
                             "pA1": "string", "pA2": "string"})
    bim["akey"] = np.where(bim["pA1"] < bim["pA2"], bim["pA1"] + bim["pA2"],
                           bim["pA2"] + bim["pA1"])
    return bim[["chr", "pos", "rsid", "akey"]]


def run(raw: Path, out: Path) -> None:
    ftd = pd.read_csv(raw, sep=r"\s+", engine="c",
                      dtype={"CHR": "int32", "BP": "int64", "SNP": "string",
                             "PVALUE": "float64", "A1": "string", "OR": "float64",
                             "SE": "float64"})
    print(f"raw FTD rows: {len(ftd):,}")
    parts = ftd["SNP"].str.split(":", expand=True)  # chr:pos:ref:alt
    ref, alt = parts[2].str.upper(), parts[3].str.upper()
    a1 = ftd["A1"].str.upper()
    # A2 = the allele of {ref, alt} that is not the effect allele A1
    a2 = np.where(a1 == alt, ref, alt)
    ftd = ftd.assign(A1=a1, A2=pd.array(a2, dtype="string"))
    ftd["akey"] = np.where(ftd["A1"] < ftd["A2"], ftd["A1"] + ftd["A2"],
                           ftd["A2"] + ftd["A1"])
    ftd["BETA"] = np.log(ftd["OR"])

    panel = _load_panel()
    merged = ftd.merge(panel, left_on=["CHR", "BP", "akey"],
                       right_on=["chr", "pos", "akey"], how="inner")
    merged = merged.drop_duplicates("rsid")
    print(f"mapped to panel rsID (allele-matched): {len(merged):,} "
          f"({100*len(merged)/len(ftd):.1f}% of raw)")

    out_df = merged[["rsid", "CHR", "BP", "A1", "A2", "BETA", "SE", "PVALUE"]].rename(
        columns={"rsid": "SNP"})
    out_df = out_df[out_df["PVALUE"].between(0, 1) & out_df["SE"].gt(0)]
    out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(out_df):,} harmonized SNPs -> {out}")


def main() -> None:
    p = argparse.ArgumentParser(description="Harmonize IFGC FTD sumstats to panel rsIDs.")
    p.add_argument("--raw", type=Path, default=RAW)
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args()
    run(args.raw, args.out)


if __name__ == "__main__":
    main()
