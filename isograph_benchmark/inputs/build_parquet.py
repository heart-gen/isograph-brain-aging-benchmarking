from __future__ import annotations

import gzip
import re
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel


def _safe_region(name: str) -> str:
    clean = re.sub(r"^Brain - ", "", name)
    clean = re.sub(r"[^A-Za-z0-9]+", "_", clean).strip("_").lower()
    return clean


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_parquet(path, index=False, compression="zstd")


def convert_brainseq_region(region: str) -> None:
    raw = rel("inputs", "raw", "brainseq")
    out = rel("inputs", "processed", "brainseq", region)
    write_parquet(pd.read_csv(raw / "counts" / region / "tx-counts.tsv", sep="\t"), out / "tx_counts.parquet")
    write_parquet(pd.read_csv(raw / "counts" / region / "gene-counts.tsv", sep="\t"), out / "gene_counts.parquet")
    write_parquet(pd.read_csv(raw / "counts" / region / "psi-events.tsv.gz", sep="\t"), out / "psi_events.parquet")


def convert_brainseq_metadata() -> None:
    raw = rel("inputs", "raw", "brainseq")
    out = rel("inputs", "processed", "brainseq", "metadata")
    meta = pd.read_csv(raw / "metadata" / "libd_rnaseq_metadata.tab", sep="\t")
    write_parquet(meta, out / "libd_rnaseq_metadata.parquet")
    for region_dir in (raw / "metadata").iterdir():
        if not region_dir.is_dir():
            continue
        frames = [pd.read_csv(path, sep="\t") for path in region_dir.glob("*.tsv")]
        if frames:
            write_parquet(pd.concat(frames, ignore_index=True), out / f"{region_dir.name}_rnaseq_metrics.parquet")


def brain_sample_ids() -> dict[str, list[str]]:
    attrs = pd.read_csv(
        rel("inputs", "raw", "gtex_v11", "metadata", "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"),
        sep="\t",
        usecols=["SAMPID", "SMTS", "SMTSD"],
    )
    brain = attrs.loc[attrs["SMTSD"].astype(str).str.startswith("Brain -")].copy()
    return {
        tissue: group["SAMPID"].tolist()
        for tissue, group in brain.groupby("SMTSD", sort=True)
    }


def convert_gtex_transcript_tpm() -> None:
    ids_by_tissue = brain_sample_ids()
    raw_path = rel(
        "inputs", "raw", "gtex_v11", "counts",
        "GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz",
    )
    header = pd.read_csv(raw_path, sep="\t", nrows=0).columns.tolist()
    for tissue, sample_ids in ids_by_tissue.items():
        cols = ["transcript_id", "gene_id", *[s for s in sample_ids if s in header]]
        df = pd.read_csv(raw_path, sep="\t", usecols=cols)
        write_parquet(df, rel("inputs", "processed", "gtex_v11", _safe_region(tissue), "transcript_tpm.parquet"))


def convert_gtex_gct(file_name: str, output_name: str) -> None:
    ids_by_tissue = brain_sample_ids()
    raw_path = rel("inputs", "raw", "gtex_v11", "counts", file_name)
    header = pd.read_csv(raw_path, sep="\t", skiprows=2, nrows=0).columns.tolist()
    for tissue, sample_ids in ids_by_tissue.items():
        cols = ["Name", "Description", *[s for s in sample_ids if s in header]]
        df = pd.read_csv(raw_path, sep="\t", skiprows=2, usecols=cols)
        write_parquet(df, rel("inputs", "processed", "gtex_v11", _safe_region(tissue), f"{output_name}.parquet"))


def _age_midpoint(age_bin: str) -> float:
    lo, hi = age_bin.split("-")
    return (int(lo) + int(hi)) / 2


def convert_gtex_sample_metadata() -> None:
    attrs = pd.read_csv(
        rel("inputs", "raw", "gtex_v11", "metadata", "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"),
        sep="\t",
        dtype={"SMGTC": str},
    )
    subj = pd.read_csv(
        rel("inputs", "raw", "gtex_v11", "metadata", "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"),
        sep="\t",
        usecols=["SUBJID", "SEX", "AGE"],
    )
    subj["AGE"] = subj["AGE"].apply(_age_midpoint)
    attrs["SUBJID"] = attrs["SAMPID"].str.split("-").str[:2].str.join("-")
    attrs = attrs.merge(subj, on="SUBJID", how="left")
    brain = attrs.loc[attrs["SMTSD"].astype(str).str.startswith("Brain -")].copy()
    for tissue, group in brain.groupby("SMTSD", sort=True):
        write_parquet(group, rel("inputs", "processed", "gtex_v11", _safe_region(tissue), "sample_attributes.parquet"))


def main() -> None:
    for region in ["caudate", "hippocampus", "dlpfc"]:
        convert_brainseq_region(region)
    convert_brainseq_metadata()
    convert_gtex_sample_metadata()
    convert_gtex_transcript_tpm()
    convert_gtex_gct("GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.gct.gz", "gene_reads")
    convert_gtex_gct("GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz", "gene_tpm")


if __name__ == "__main__":
    main()
