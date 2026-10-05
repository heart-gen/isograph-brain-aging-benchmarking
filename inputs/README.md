# Input data processing

This directory holds the real-data inputs for the IsoGraph brain-aging analyses and
documents the path from raw RNA-seq sources to the self-contained IsoGraph dataset
bundles that every downstream analysis consumes. The processing *code* lives in the
`isograph_benchmark.inputs` package (`copy_raw.py`, `build_parquet.py`,
`build_bundles.py`, plus the helpers `build_transcript_annotation.py` and
`build_gtex_junction_usage.py`); this `inputs/` tree holds its *data products* and a few
static reference resources. The whole pipeline
is reproducible with one script:

```bash
cd <repo-root>
bash inputs/_h/build_data_pipeline.sh           # full run
bash inputs/_h/build_data_pipeline.sh --skip-raw --skip-pcs   # from the Zenodo deposit
```

Raw source files are intentionally git-ignored and are planned for Zenodo archiving
(DOI pending; they are not yet in `zenodo/MANIFEST.tsv`, which lists only heavy
module-discovery outputs). Steps 1–3 below are fully reproducible from that deposit; the controlled
-access genotype step (SNP PCs) ships its precomputed output in the deposit.

## Directory layout

| Path | Contents | Tracked? |
|---|---|---|
| `inputs/raw/` | Raw BrainSEQ TSVs and GTEx GCT/TXT count matrices + metadata (not present in a fresh clone) | git-ignored (Zenodo) |
| `inputs/processed/` | Per-region Parquet count/PSI matrices and sample metadata (`brainseq/`, `gtex_v11/`); also holds the SNP-PC step under `brainseq/genetic_similarity/` | tracked, except `gtex_v11/*/junction_usage.parquet` |
| `inputs/bundles/` | Self-contained IsoGraph dataset bundles (`brainseq_v1/`, `brainseq_sczd/`, `gtex_v11_brain/`, one dir per region) | tracked |
| `inputs/go_annotations/` | GO annotations (`go-basic.obo`, `goa_human.gaf.gz`, `Homo_sapiens.gene_info.gz`) for module enrichment | tracked |
| `inputs/rbp_motifs/` | ATtRACT RBP motif database and PWMs (`ATtRACT_db.txt`, `pwm.txt`) used by the RBP scans | tracked (downloaded `*.zip` ignored) |
| `inputs/_m/longread_aged_dlpfc/` | Public long-read (ONT, Bambu) DLPFC transcript counts + `provenance.json`, used for orthogonal switch confirmation | tracked (`inputs/_m/*.log` ignored) |
| `inputs/tin/` | Derived TIN caches (regenerable) | git-ignored |
| `inputs/_h/` | Build/acquisition scripts (below) | tracked |

Scripts in `inputs/_h/`:

| Script | Purpose |
|---|---|
| `build_data_pipeline.sh` | Steps 0–3 below |
| `rebuild_gtex_counts.sh` | SLURM rebuild of the GTEx bundles on the count scale |
| `build_gtex_junction_usage.sh` | SLURM array (one task per region): GTEx STAR junction counts → `inputs/processed/gtex_v11/<region>/junction_usage.parquet` (GTEx analogue of BrainSEQ PSI; git-ignored, regenerable) |
| `download_neuronal_clip.sh` | SLURM download of public neuronal CLIP inputs into `inputs/raw/neuronal_clip/` (config: `configs/neuronal_clip.yaml`) |
| `EGA_ACCESS_REQUEST.md` | Notes for the optional controlled-access EGA eCLIP datasets |

## Data sources

- **BrainSEQ** — LIBD postmortem brain RNA-seq (Salmon transcript counts, gene counts,
  PSI splice events) for caudate (Phase 3), hippocampus and DLPFC (Phase 2), plus a
  master sample manifest and per-region QC metrics.
- **GTEx v11** — 13 brain regions, RSEM transcript expected counts and TPM, RNASeQC
  gene reads and TPM, and sample/subject annotations. Exact subject age is resolved
  from the v8 phenotype file where available and falls back to v11 age-band midpoints.

## Pipeline

The pipeline runs in four steps, orchestrated by `inputs/_h/build_data_pipeline.sh`:

**Step 0 — Copy raw data** (`isograph_benchmark.inputs.copy_raw`). Copies raw files from
HPC storage into `inputs/raw/` and writes a SHA-256 manifest
(`reports/raw_copy_manifest.parquet`). Source paths are read from `configs/data_sources.yaml`. Skipped with `--skip-raw` when `inputs/raw/` is
already populated from the Zenodo deposit.

**Step 1 — Convert to Parquet** (`isograph_benchmark.inputs.build_parquet`). Reads
`inputs/raw/` and writes compressed Parquet to `inputs/processed/`:

- BrainSEQ per region: `tx_counts.parquet`, `gene_counts.parquet`, `psi_events.parquet`;
  plus `libd_rnaseq_metadata.parquet` and per-region QC metrics under `metadata/`.
- GTEx per region (13 regions): `transcript_reads.parquet` (RSEM expected counts),
  `transcript_tpm.parquet`, `gene_reads.parquet`, `gene_tpm.parquet`, and
  `sample_attributes.parquet`.

**Step 2 — BrainSEQ SNP principal components**
(`inputs/processed/brainseq/genetic_similarity/_h/compute_snp_pcs.sh`, plink2). LD-prunes
each autosome (MAF ≥ 0.05, geno ≤ 0.05, HWE p > 1e-6, r² < 0.2), merges chr1–chr22, and
computes 10 PCs, written to
`inputs/processed/brainseq/genetic_similarity/_m/TOPMed_LIBD.eigenvec` (and
`.eigenval`). Genotypes are
controlled access; skip with `--skip-pcs` (samples then receive NaN PC covariates).
Runtime ≈ 20–30 min on 16 cores.

**Step 3 — Build IsoGraph bundles** (`isograph_benchmark.inputs.build_bundles`). Reads
`inputs/processed/` (and the SNP PCs) and writes one bundle per region under
`inputs/bundles/`. Expression filtering is applied during this step:

- BrainSEQ: gene **CPM ≥ 1** in ≥ max(10, 10% of samples), on Salmon gene counts.
- GTEx: gene **CPM ≥ 1** in ≥ max(10, 10% of samples), on RNASeQC gene read counts.

BrainSEQ bundle builds additionally require the transcript annotation described under
[Dependencies and gotchas](#dependencies-and-gotchas).

Both suites use the same CPM ≥ 1 filter. GTEx bundles are built from RSEM expected
counts (count scale), so the same library-size-normalized CPM threshold applies as for
BrainSEQ — this replaced an earlier TPM ≥ 0.1 filter (see [GTEx counts rebuild](#gtex-counts-rebuild)).

## Bundle format

Each bundle is a self-contained IsoGraph dataset directory:

| File | Contents |
|---|---|
| `manifest.json` | Schema, dimensions, filters, and provenance (source, expression filter, feature counts before/after filtering, SNP-PC note) |
| `samples.parquet` | Sample metadata: `sample_id` + all covariates and trait columns |
| `genes.parquet` | Gene feature table (`gene_id`, genomic coordinates) |
| `transcripts.parquet` | Transcript feature table (`transcript_id` → `gene_id`) |
| `gene_counts.npz` | Gene count matrix `[genes × samples]` |
| `transcript_counts.npz` | Transcript count matrix `[transcripts × samples]` — providing this activates IsoGraph's abundance channel alongside the switch channel |

Example provenance (BrainSEQ caudate): 78,932 → **20,365 genes** and 384,354 →
**232,176 transcripts** after the CPM ≥ 1 filter, 238 Control adult samples.

## Bundle inventory

| Suite | Region(s) | Cohort / filters | Trait |
|---|---|---|---|
| `brainseq_v1/` | `caudate` (Phase 3), `hippocampus`, `dlpfc` (Phase 2) | Dx = Control, not dropped, Age ≥ 18 | Age |
| `brainseq_sczd/` | `caudate` | Dx = Control + SCZD, not dropped, Age ≥ 18 (disease associations) | Dx |
| `gtex_v11_brain/` | 13 brain regions (`amygdala` … `substantia_nigra`) | GTEx v11 brain | AGE (exact; v8 preferred) |

Covariates: BrainSEQ bundles carry Sex, MoD, RIN, mapping rate, mito rate, and
SNP PCs (`SNP_PC1`–`SNP_PC10`; samples without genotypes get NaN); GTEx bundles carry SEX, SMRIN, SMTSISCH, and SMMAPRT.

## GTEx counts rebuild

`inputs/_h/rebuild_gtex_counts.sh` (SLURM) regenerates the 13 GTEx bundles from RSEM
transcript **expected counts** (count scale) with a **gene CPM ≥ 1** filter — the same
threshold as BrainSEQ — replacing the earlier TPM-based bundles that used a TPM ≥ 0.1
filter. Moving to the count scale plus the CPM filter resolves a VAE divergence seen on
TPM input and an out-of-memory failure from the looser TPM ≥ 0.1 filter (~2× the
features). It runs Step 1's `convert_gtex_transcript_reads()` then rebuilds each region
with `build_gtex_bundle()`, which applies the CPM filter via the shared
`_brainseq_expressed_genes` helper.

## Dependencies and gotchas

- **Transcript annotation.** `build_bundles` reads
  `inputs/raw/brainseq/annotations/transcript-annotation.tsv`, which is not tracked.
  Regenerate it from the GENCODE v47 primary-assembly GTF with
  `python -m isograph_benchmark.inputs.build_transcript_annotation --verify`; `--verify`
  fails on any mismatch with a committed bundle's `transcripts.parquet`.
- **Stale references in scripts.** The header comments of `build_data_pipeline.sh` still
  describe the GTEx filter as TPM ≥ 0.1 and list `transcript_tpm` as the GTEx transcript
  input; the code (`build_bundles.build_gtex_bundle`) uses RSEM expected counts with
  CPM ≥ 1. `compute_snp_pcs.sh`'s own header and the bundle manifests'
  `snp_pcs` provenance string point at `inputs/_h/compute_snp_pcs.sh`; the script lives under
  `inputs/processed/brainseq/genetic_similarity/_h/`.
- **Junction usage** (`junction_usage.parquet`) is not part of Steps 0–3, is git-ignored,
  and needs the GTEx junction GCT in `inputs/raw/`; run `build_gtex_junction_usage.sh`
  after the GTEx bundles exist (it restricts to bundle samples).

## Data availability

The data sources are referenced at build time through `configs/data_sources.yaml`, which
points at machine-specific paths (HPC and local); that file is environment-specific and is not the
authoritative record of where the data come from. The authoritative pointers are:

| Resource | In this repository? | Where to obtain |
|---|---|---|
| Processed per-region matrices (`inputs/processed/`) | **Yes** — tracked Parquet | Included here; regenerable from raw via Steps 0–1 |
| IsoGraph dataset bundles (`inputs/bundles/`) | **Yes** — tracked Parquet/`.npz` | Included here; regenerable via Step 3 |
| Raw BrainSEQ + GTEx count matrices and metadata (`inputs/raw/`) | No — git-ignored | Zenodo deposit, DOI `<pending>` (archived at publication) |
| GTEx v11 expression and sample annotations | No | GTEx Portal `<https://gtexportal.org>`; protected-access data via AnVIL/dbGaP accession `<pending>` |
| BrainSEQ (LIBD) RNA-seq and phenotypes | No | LIBD / BrainSEQ consortium; controlled-access via dbGaP accession `<pending>` |
| BrainSEQ TOPMed-imputed genotypes (for SNP PCs) | No — controlled access | dbGaP accession `<pending>`; precomputed PCs included in the Zenodo deposit |
| IsoGraph software | No — separate repository | `<https://github.com/heart-gen/IsoGraph>`, PyPI package `isograph`, software DOI `<pending>` |

Steps 1–3 of the pipeline are fully reproducible from the Zenodo deposit once it is
released; the controlled-access genotype step (SNP PCs) ships its precomputed output so
no protected data are required to rebuild the bundles. Replace each `<pending>`
placeholder with the final DOI/accession before manuscript submission.

## Reproducibility notes

- Raw files are git-ignored; the Step 0 manifest records a SHA-256 checksum per file so
  the Zenodo deposit can be verified bit-for-bit.
- `inputs/processed/` and `inputs/bundles/` are tracked so analyses are reproducible
  without re-downloading raw data; large matrices are stored as Parquet/`.npz`.
- The pipeline is import-guarded (`python3 -c "import isograph_benchmark"`) and must be
  run from the repo root with the project environment active
  (`/ocean/projects/bio260021p/shared/opt/envs/isograph`).
- Requirements: Python ≥ 3.11 (project env), R ≥ 4.3 (downstream), plink2 ≥ 2.00 (Step 2).
