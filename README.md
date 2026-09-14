# IsoGraph — brain aging benchmarking

Benchmarking and brain-aging analysis for **IsoGraph**, an isoform-switch network
method. The repository is laid out as the manuscript's argument: each numbered stage
answers one question and feeds the next.

> **North star.** IsoGraph is a **complementary isoform-switch network method**, not a
> globally superior one. The de-confounded gene-level test shows abundance dominates;
> IsoGraph's defensible value is a small, specific **DTU-without-DGE** layer that is
> structurally invisible to any DGE/WGCNA pipeline. Every claim in this repo is written
> to that bound.

## Stages

| Stage | Question | Headline output |
|---|---|---|
| [`01_synthetic_benchmark/`](01_synthetic_benchmark/README.md) | Does the method recover switch modules where truth is known? | Wins where switching dominates, loses when abundance-dominated or degraded: module recovery favours IsoGraph in **7 of 15 scenarios** — **Fig 1** |
| [`02_module_discovery/`](02_module_discovery/README.md) | What modules exist in postmortem human brain? | IsoGraph fits + 3 matched WGCNA baselines; the cohort × region artifact store |
| [`03_module_trust/`](03_module_trust/README.md) | Are they reproducible, or is a fine partition noise? | 250/266 chance-trusted; 23 aging replications — **Fig 2** |
| [`04_module_characterization/`](04_module_characterization/README.md) | What are they, and do they add signal beyond abundance? | 6 of 8 disease modules GO-invisible, 4 with near-complete switch coverage; axes near-orthogonal; **the scope bound** |
| [`05_genetic_anchoring/`](05_genetic_anchoring/README.md) | Are they genetically real? Where does disease risk land? | Supporting, set level: splicing-specificity contrast 1.111, IsoGraph-only, **not supported per gene**; 42 signal-level coloc nominations (SNCA the worked example) — **Fig 3, Fig 4** |
| [`06_switch_mechanism/`](06_switch_mechanism/README.md) | Are the switches real, and what do they do? | Productive UTR/CDS remodeling, not decay |
| [`07_rbp_regulation/`](07_rbp_regulation/README.md) | What trans factors could drive the co-switching? | 714 module × RBP hits on raw counts, 310 under the opportunity-adjusted GLM — but only 43 in common (candidates, not binding data) |
| [`manuscript/`](manuscript/README.md) | — | Every display item + its builder |

[`ANALYSIS_MAP.md`](ANALYSIS_MAP.md) is the one-row-per-analysis index:
analysis → CLI → wrapper → outputs → display item.

## Supporting directories

- `inputs/` — `raw/` (gitignored local copies, Zenodo-bound), `processed/` parquet,
  `bundles/` IsoGraph dataset bundles, plus GO annotations, RBP motifs and TIN caches.
- `isograph_benchmark/` — the Python package: every analysis is a committed,
  parametrized CLI module here, invoked by a stage wrapper.
- `configs/` — single source of truth for grids, covariates and traits.
- `env/`, `scripts/`, `tests/`, `zenodo/`, `develop/` (local scratch, gitignored).

## Conventions

Each stage follows the group layout: `_h/` holds the SLURM wrappers, numbered `01..N`
in **dependency order**; `_m/` holds its outputs; `README.md` explains the stage.

**Paths are never hardcoded.** `isograph_benchmark/paths.py` holds `OUTPUT_DIRS`, the
single definition of where each stage writes; code addresses stages by logical bucket
(`stage_out`, `region_store`, `region_artifact_dir`). Moving a stage is a one-line edit.

Heavy regenerable artifacts (`feature_scores`, `feature_reconstruction`,
`high_vs_low_table`, `edges`, coloc/LDSC per-locus intermediates) are gitignored and
distributed via Zenodo — see `zenodo/`. Lean result tables are tracked in git-LFS.

## Running

Heavy compute is SLURM-only (account `bio260021p`; memory is `--cpus-per-task` × 2000M,
so do not pass `--mem`). Login nodes are for reads, small aggregation and plotting.
Submit from the repo root, or set `ISOGRAPH_BENCHMARK_ROOT`:

```bash
sbatch 05_genetic_anchoring/_h/01.qtl_anchoring.sh
```

Interpreters:

| Use | Path |
|---|---|
| Python | `/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python` |
| Figure R | `/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript` |
| WGCNA + GWAS R | `/ocean/projects/bio250020p/shared/opt/env/R_env` |

Terminology: synthetic *nonlinear* settings are interactions in feature space; real-data
*spline aging* analyses are spline models of age against module eigengenes.
