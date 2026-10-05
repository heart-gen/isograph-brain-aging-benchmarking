# IsoGraph brain analysis

Analysis and benchmarking repository for **IsoGraph resolves coordinated transcript
choice and gene abundance in the aging human brain**, by Alexis Bennett, Elisa Kain
Johnson, and Kynon J. M. Benjamin.

- [Manuscript source](https://github.com/heart-gen/isograph-brain-manuscript)
- [Read the manuscript](https://heart-gen.github.io/isograph-brain-manuscript/) · [PDF](https://heart-gen.github.io/isograph-brain-manuscript/manuscript.pdf)
- [IsoGraph software](https://github.com/heart-gen/IsoGraph) · [software archive cited in the manuscript](https://doi.org/10.5281/zenodo.21707653)

Release **v1.0.0** captures the initial manuscript submission. The companion
manuscript metadata and framing are aligned to manuscript commit
[`14a74b1fa6faf32c9b2c008ab467236f70a30136`](https://github.com/heart-gen/isograph-brain-manuscript/tree/14a74b1fa6faf32c9b2c008ab467236f70a30136).
See [release notes](RELEASE_NOTES.md) and [archiving instructions](zenodo/README.md).
Preprint and analysis-archive DOIs are pending.

IsoGraph represents relative transcript usage and gene abundance as paired network
channels. The analyses evaluate synthetic recovery, brain aging associations,
module stability and transfer, held-out differential transcript usage (DTU),
transcript structure, and genetic and orthogonal support. Recovery depends on the
signal regime: abundance and degradation settings can favor abundance-based methods.
Disease analyses prioritize candidate transcript events, including alternative
terminal-exon usage at **PRDM2** associated with amyotrophic lateral sclerosis.

## Analyses

| Stage | Contents |
|---|---|
| [`01_synthetic_benchmark/`](01_synthetic_benchmark/README.md) | Synthetic module recovery, interpretation, robustness, genetics and scale benchmarks |
| [`02_module_discovery/`](02_module_discovery/README.md) | BrainSEQ and GTEx module fits and matched co-expression baselines |
| [`03_module_characterization/`](03_module_characterization/README.md) | Usage and abundance contributions, enrichment and module interpretation |
| [`04_module_trust/`](04_module_trust/README.md) | Donor-subset stability, cross-cohort transfer and functional preservation |
| [`05_genetic_anchoring/`](05_genetic_anchoring/README.md) | QTL anchoring, allelic evidence and disease colocalization |
| [`06_switch_mechanism/`](06_switch_mechanism/README.md) | Transcript architecture, junction and independent long-read support |
| [`07_rbp_regulation/`](07_rbp_regulation/README.md) | Candidate RBP regulation and binding-support analyses |
| [`08_integration/`](08_integration/README.md) | Integrated gene/event evidence and held-out DTU analyses |
| [`manuscript/`](manuscript/README.md) | Figure and table builders, results and supplementary exports |

Final display numbering and legends are maintained in the companion manuscript.
Stage READMEs also describe historical analyses and intermediate outputs.

## Supporting directories

- `inputs/` — `raw/` (gitignored local copies, Zenodo-bound), `processed/` parquet,
  `bundles/` IsoGraph dataset bundles, plus GO annotations, RBP motifs and TIN caches.
- `isograph_benchmark/` — the Python package: every analysis is a committed,
  parametrized CLI module here, invoked by a stage wrapper.
- `configs/` — single source of truth for grids, covariates and traits.
- `00_scripts/` — including `slurm_dag.sh`, the dependency-graph submitter behind every `run_stage.sh`.
- `env/`, `tests/`, `zenodo/`, `develop/` (local scratch, gitignored).

## Conventions

Each stage follows the group layout: `_h/` holds the SLURM wrappers, `_m/` holds its
outputs, `README.md` explains the stage.

**Names encode run order.** Stages are ordered so that every stage reads only earlier
stages. Within a stage, a wrapper's leading number is its **tier** and the letter tells
siblings apart (`01a`, `01b`, `02a`, …): a step reads only tiers before its own, so every
step of a tier can run in parallel once the earlier tiers have finished. The few same-tier
steps that must wait on a sibling — because both write one shared file — say so in the
stage's `run_stage.sh`. A later-tier mode of an analysis (a report, an aggregate, a meta
step) is its own numbered wrapper rather than a flag on an earlier one. Helper scripts a
wrapper calls (`stability_wgcna.R`, `ldsc_wrapper.py`, …) are unnumbered. Retired wrappers
live unnumbered in `_h/retired/`; their outputs survive only at the git tag named in the
stage README.

**The run order is committed, not described.** `<stage>/_h/run_stage.sh` submits that
stage as a SLURM dependency graph, and `00_scripts/run_pipeline.sh` chains the stages.

**Paths are never hardcoded.** `isograph_benchmark/paths.py` holds `OUTPUT_DIRS`, the
single definition of where each stage writes; code addresses stages by logical bucket
(`stage_out`, `region_store`, `region_artifact_dir`). Moving a stage is a one-line edit.

Heavy regenerable artifacts (`feature_scores`, `feature_reconstruction`,
`high_vs_low_table`, `edges`, coloc/LDSC per-locus intermediates) are gitignored and
intended for a separate Zenodo deposit — see `zenodo/`. Lean result tables are tracked in git-LFS.

## Running

Heavy compute is SLURM-only (account `bio260021p`; memory is `--cpus-per-task` × 2000M,
so do not pass `--mem`). Login nodes are for reads, small aggregation and plotting.
Run from the repo root, or set `ISOGRAPH_BENCHMARK_ROOT`:

```bash
bash 00_scripts/run_pipeline.sh --dry-run                      # every stage's plan, nothing submitted
bash 00_scripts/run_pipeline.sh --stages 05,06,07,08           # submit stages, each waiting on the last
bash 05_genetic_anchoring/_h/run_stage.sh --dry-run # one stage's plan
bash 05_genetic_anchoring/_h/run_stage.sh --from 06 # resume a stage at a tier
sbatch 05_genetic_anchoring/_h/01a.qtl_anchoring.sh # a single step
```

Two steps need a login node because compute nodes have no outbound network
(`06_switch_mechanism/_h/01f`, `07_rbp_regulation/_h/03b`). The runners hold whatever
waits on them and print the command; rerun with `--login-done <id>` once it has run.
Options for every runner are documented in `00_scripts/slurm_dag.sh`.

Inputs are built by `inputs/_h/build_data_pipeline.sh`; the synthetic benchmark in
`01_synthetic_benchmark/` is archived and not part of `00_scripts/run_pipeline.sh`.

Interpreters:

| Use | Path |
|---|---|
| Python | `/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python` |
| Figure R | `/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript` |
| WGCNA + GWAS R | `/ocean/projects/bio250020p/shared/opt/env/R_env` |

Terminology: synthetic *nonlinear* settings are interactions in feature space; real-data
*spline aging* analyses are spline models of age against module eigengenes.

## Citation and licensing

[CITATION.cff](CITATION.cff) provides the repository version, author order, ORCIDs,
and companion manuscript reference. Cite this repository and the manuscript;
add the version-specific archive DOI to citations once the archive is deposited.
The IsoGraph package DOI above identifies the separate software package.

Original analysis code is licensed under [Apache-2.0](LICENSE). Original derived
data, tables, figures and documentation are licensed under
[CC-BY-4.0](LICENSE-DATA.md). Third-party resources and controlled-access inputs
retain their original terms; these licenses do not grant rights to redistribute them.

## Obtaining result artifacts

Clone with Git LFS installed, then retrieve the tracked artifacts:

```bash
git clone https://github.com/heart-gen/isograph-brain-aging-benchmarking.git
cd isograph-brain-aging-benchmarking
git checkout v1.0.0
git lfs pull
```

Environment snapshots are in [`env/`](env/). Source-data access and preparation
are described in [`inputs/README.md`](inputs/README.md). A Git source archive may
contain LFS pointers rather than artifact contents; see the archive instructions
for verifying a complete deposit. Large ignored intermediates require the separate
data deposit or regeneration from authorized source inputs.
