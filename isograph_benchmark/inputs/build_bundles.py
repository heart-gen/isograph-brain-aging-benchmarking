from __future__ import annotations

import numpy as np
import pandas as pd

from isograph.io.artifacts import DatasetBundle, build_feature_spec, build_matrix_spec, save_dataset_bundle
from isograph.validation import DatasetManifest
from isograph_benchmark.paths import ensure_dir, rel

ADULT_AGE_MIN = 18.0
SNP_PC_COLS = [f"SNP_PC{i}" for i in range(1, 11)]

# Map region name → per-region QC metrics parquet
_REGION_METRICS = {
    "caudate": "caudate_rnaseq_metrics.parquet",
    "hippocampus": "hippocampus_rnaseq_metrics.parquet",
    "dlpfc": "dlpfc_rnaseq_metrics.parquet",
}

# Restrict each region to its BSP2 (RiboZeroGold) dataset only.
# DLPFC metadata contains both BrainSeq_Phase1 (PolyA, BSP1) and
# BrainSeq_Phase2_DLPFC (RiboZeroGold, BSP2) with overlapping RNums;
# keep only BSP2 so library prep matches hippocampus and caudate.
_REGION_DATASET = {
    "caudate": "BrainSeq_Phase3_Caudate",
    "hippocampus": "BrainSeq_Phase2_HIPPO",
    "dlpfc": "BrainSeq_Phase2_DLPFC",
}


def _matrix_from_wide(df: pd.DataFrame, id_cols: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    feature = df[id_cols].copy()
    matrix = df.drop(columns=id_cols).to_numpy(dtype=float)
    return feature, matrix


def _matrix_from_wide_subset(
    df: pd.DataFrame, id_cols: list[str], keep_samples: list[str]
) -> tuple[pd.DataFrame, np.ndarray]:
    feature = df[id_cols].copy()
    matrix = df[keep_samples].to_numpy(dtype=float)
    return feature, matrix


def _load_snp_pcs() -> pd.DataFrame:
    path = rel("inputs", "processed", "brainseq", "genetic_similarity", "_m", "TOPMed_LIBD.eigenvec")
    if not path.exists():
        raise FileNotFoundError(
            f"SNP PC file not found: {path}\n"
            "Run: bash inputs/processed/brainseq/genetic_similarity/_h/compute_snp_pcs.sh"
        )
    pcs = pd.read_csv(path, sep="\t")
    # plink2 eigenvec header: '#IID  PC1  PC2  ...' (no FID when psam has only IID)
    pcs = pcs.rename(columns={"#IID": "BrNum"})
    pc_rename = {f"PC{i}": f"SNP_PC{i}" for i in range(1, 11)}
    pcs = pcs.rename(columns=pc_rename)
    return pcs[["BrNum"] + SNP_PC_COLS]


def _load_brainseq_qc_metrics(region: str) -> pd.DataFrame:
    name = _REGION_METRICS[region]
    path = rel("inputs", "processed", "brainseq", "metadata", name)
    metrics = pd.read_parquet(path, columns=["sample_rnum", "mapping_rate", "mito_mapped", "total_reads", "r_rna_rate"])
    metrics["mito_rate"] = metrics["mito_mapped"] / metrics["total_reads"].clip(lower=1)
    return metrics[["sample_rnum", "mapping_rate", "mito_rate", "r_rna_rate"]]


def build_brainseq_bundle(
    region: str,
    allowed_diagnoses: list[str] | None = None,
    suite_name: str | None = None,
    out_bundle_name: str | None = None,
) -> None:
    if allowed_diagnoses is None:
        allowed_diagnoses = ["Control"]
    if suite_name is None:
        suite_name = "brainseq_v1"
    if out_bundle_name is None:
        out_bundle_name = region

    src = rel("inputs", "processed", "brainseq", region)
    tx = pd.read_parquet(src / "tx_counts.parquet")
    gene = pd.read_parquet(src / "gene_counts.parquet")

    all_sample_ids = tx.columns[3:].tolist()  # cols: Name, Length, EffectiveLength, <RNums...>

    meta = pd.read_parquet(
        rel("inputs", "processed", "brainseq", "metadata", "libd_rnaseq_metadata.parquet")
    )
    pcs = _load_snp_pcs()
    metrics = _load_brainseq_qc_metrics(region)

    dataset_name = _REGION_DATASET[region]
    samples = meta.loc[meta["RNum"].isin(all_sample_ids)].copy()
    samples = samples[
        (samples["Dataset"] == dataset_name)
        & (samples["Dx"].isin(allowed_diagnoses))
        & (samples["dropped"] == "f")
        & (samples["Age"] >= ADULT_AGE_MIN)
    ].copy()

    samples["MoD"] = (
        samples["MoD"]
        .replace("No autopsy performed", "Undetermined")
        .replace(".", "Undetermined")
        .fillna("Undetermined")
    )

    samples = samples.merge(pcs, on="BrNum", how="left")
    samples = samples.merge(metrics, left_on="RNum", right_on="sample_rnum", how="left")
    samples = samples.drop(columns=["sample_rnum"], errors="ignore")

    filtered_ids = [r for r in all_sample_ids if r in set(samples["RNum"])]
    samples = samples.set_index("RNum").loc[filtered_ids].reset_index().rename(columns={"RNum": "sample_id"})

    tx_annot = pd.read_csv(
        rel("inputs", "raw", "brainseq", "annotations", "transcript-annotation.tsv"), sep="\t"
    )
    tx_feature, tx_matrix = _matrix_from_wide_subset(tx, ["Name", "Length", "EffectiveLength"], filtered_ids)
    tx_feature = tx_feature.rename(columns={"Name": "transcript_id"}).merge(
        tx_annot[["transcript_id", "gene_id", "transcript_name", "transcript_type"]],
        on="transcript_id",
        how="left",
    )

    gene_feature, gene_matrix = _matrix_from_wide_subset(gene, ["Geneid", "Chr", "Start", "End", "Strand", "Length"], filtered_ids)
    gene_feature = gene_feature.rename(columns={"Geneid": "gene_id"})

    dx_label = "+".join(sorted(allowed_diagnoses))
    manifest = DatasetManifest(
        dataset_name=f"brainseq_{out_bundle_name}_v1",
        suite_name=suite_name,
        description=(
            f"BrainSEQ {region} IsoGraph bundle — Dx={dx_label}, adults (Age≥{ADULT_AGE_MIN}), "
            "QC-passed, with SNP PCs"
        ),
        sample_table="samples.parquet",
        feature_tables=[
            build_feature_spec("gene", "genes.parquet", gene_feature),
            build_feature_spec("transcript", "transcripts.parquet", tx_feature),
        ],
        matrices=[
            build_matrix_spec("gene_counts", "gene_counts.npz", gene_matrix),
            build_matrix_spec("transcript_counts", "transcript_counts.npz", tx_matrix),
        ],
        provenance={
            "region": region,
            "source": "BrainSEQ",
            "filters": f"Dx={dx_label}, dropped=f, Age>={ADULT_AGE_MIN}",
            "snp_pcs": "TOPMed-imputed, computed via inputs/_h/compute_snp_pcs.sh",
        },
    )
    bundle = DatasetBundle(
        manifest=manifest,
        sample_table=samples,
        feature_tables={"gene": gene_feature, "transcript": tx_feature},
        matrices={"gene_counts": gene_matrix, "transcript_counts": tx_matrix},
        truth_tables={},
    )
    out_path = ensure_dir(rel("inputs", "bundles", suite_name, out_bundle_name))
    save_dataset_bundle(bundle, out_path)
    n_snp = samples[SNP_PC_COLS[0]].notna().sum()
    print(f"  {out_bundle_name}: {len(samples)} samples ({n_snp} with SNP PCs) → {out_path}")


def build_gtex_bundle(region_dir_name: str) -> None:
    src = rel("inputs", "processed", "gtex_v11", region_dir_name)
    tx = pd.read_parquet(src / "transcript_tpm.parquet")
    samples = pd.read_parquet(src / "sample_attributes.parquet").rename(columns={"SAMPID": "sample_id"})
    sample_ids = [col for col in tx.columns if col not in {"transcript_id", "gene_id"}]
    samples = samples.set_index("sample_id").loc[sample_ids].reset_index()
    tx_feature, tx_matrix = _matrix_from_wide(tx, ["transcript_id", "gene_id"])
    gene_feature = (
        tx_feature[["gene_id"]]
        .dropna()
        .drop_duplicates()
        .sort_values("gene_id")
        .reset_index(drop=True)
    )
    manifest = DatasetManifest(
        dataset_name=f"gtex_v11_{region_dir_name}_v1",
        suite_name="gtex_v11_brain",
        description=f"GTEx v11 {region_dir_name} transcript TPM IsoGraph bundle",
        sample_table="samples.parquet",
        feature_tables=[
            build_feature_spec("gene", "genes.parquet", gene_feature),
            build_feature_spec("transcript", "transcripts.parquet", tx_feature),
        ],
        matrices=[
            build_matrix_spec("transcript_counts", "transcript_counts.npz", tx_matrix),
        ],
        provenance={
            "region": region_dir_name,
            "source": "GTEx v11",
            "assay_note": "transcript TPM is stored as transcript_counts for IsoGraph compatibility",
        },
    )
    bundle = DatasetBundle(
        manifest=manifest,
        sample_table=samples,
        feature_tables={"gene": gene_feature, "transcript": tx_feature},
        matrices={"transcript_counts": tx_matrix},
        truth_tables={},
    )
    save_dataset_bundle(bundle, ensure_dir(rel("inputs", "bundles", "gtex_v11_brain", region_dir_name)))


def main() -> None:
    print("Building BrainSEQ bundles (controls only)...")
    for region in ["caudate", "hippocampus", "dlpfc"]:
        build_brainseq_bundle(region)
    print("\nBuilding BrainSEQ SCZD+Control bundle (caudate)...")
    build_brainseq_bundle(
        "caudate",
        allowed_diagnoses=["Control", "SCZD"],
        suite_name="brainseq_sczd",
        out_bundle_name="caudate",
    )
    print("\nBuilding GTEx v11 bundles...")
    gtex_root = rel("inputs", "processed", "gtex_v11")
    if gtex_root.exists():
        for region_dir in sorted(path for path in gtex_root.iterdir() if path.is_dir()):
            build_gtex_bundle(region_dir.name)
            print(f"  {region_dir.name}: done")


if __name__ == "__main__":
    main()
