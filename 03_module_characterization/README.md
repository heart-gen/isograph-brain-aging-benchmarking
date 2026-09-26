# 03 — Module characterization

**Question:** what *are* these modules, and do they carry information abundance
pipelines miss?

This stage feeds the honest scope bound for the whole paper: on per-module rates
IsoGraph is **not** globally superior to WGCNA. Phenotype signal lives in the switch
features (both switch-fed methods win) and GO enrichment is abundance-dominated. The one
clean method effect is `isograph` > `wgcna_multiplex` on identical features. Read every
result here against the three matched baselines from stage 02. The three-baseline rate
comparison that states the bound also reads the cross-cohort replication tables, so it runs
in stage 04 (`04_module_trust/_h/03h.baseline_comparison.sh`).

Stage order: this stage was `04_module_characterization` until 2026-09-15. It now precedes
module trust because the trust funnel reads its interpretation, enrichment and
composition-unique genes.

## Order

Run the whole stage with `bash 03_module_characterization/_h/run_stage.sh` (add `--dry-run`
to print the plan). The leading number of a wrapper is its tier; steps in one tier run in
parallel.

| Step | Wrapper | Waits on | Produces |
|---|---|---|---|
| 01a–01b | `interpret_modules_{brainseq,gtex}` | stage 02 | Switch drivers, transcript polarity, high-vs-low contrasts, GENCODE v47 structural annotation (run with `--force`) |
| 01c–01d | `module_enrichment_{brainseq,gtex}` | stage 02 | Per-module GO:BP + network metrics joined to phenotype FDR, all methods |
| 01e–01f | `incremental_association_{brainseq,gtex}` | stage 02 | De-confounded test: does the switch channel add signal *beyond* gene-level abundance? |
| 01g–01h | `celltype_composition_{brainseq,gtex}` | stage 02 | Cell-type composition (MuSiC deconvolution; GTEx via `gtex_music_deconv.R`) and the composition-adjusted incremental arm |
| 02a | `go_invisible_gate` | 01a, 01c | The GO-invisible disease-module gate (BrainSEQ SCZD) |
| 02b | `characterize_composition_unique` | 01c, 01e | Composition-unique gene sets |
| 02c | `composition_meta` | 01g, 01h | With-vs-without composition adjustment rollup across BrainSEQ + GTEx (`celltype_composition meta`) |
| 03a | `abundance_structure` | 01e, 01f, 02b | Abundance/switch orthogonality for the BrainSEQ SCZD store, the pooled incremental summary and the example gene (figSeparation B-C) |
| 03b | `abundance_structure_all` | 03a | Abundance/switch orthogonality for **all 17** stores (array; `--orthogonality-only`) |
| 03c | `abundance_structure_rollup` | 03a, 03b | *Local, not in the DAG.* Pools the 17 tables into `_m/axis_orthogonality_{all,summary}.parquet` + `AXIS_ORTHOGONALITY.md` and rebuilds figSeparation (panel A faceted by analysis) |

Retired (`_h/retired/`; outputs only at tag `legacy_expression_filter`): `project_tiers_pilot`,
`tiers_fanout` (tier projection, cited nowhere).

## Covariate policy

Residualization is **discovery-only**: confounds are regressed out when learning the
network, and trait inference happens downstream on raw `feature_scores`
(`residualize_composition=False`). Do not re-residualize at inference — that double-corrects.

## Key results

- **6 of 8** phenotype-significant SCZD switch modules (FDR ≤ 0.10) are **GO-invisible**. Four
  carry real anticorrelated transcript pairs in nearly every member (M026 23/23, M022 29/29,
  M023 27/27, M010 65/83) with functional consequence at or above background; M025 (11/24) is
  partial and M020 (4/30) weak. The two GO-visible disease modules switch comparably (M012
  40/58, M011 59/65), so GO-invisibility reflects GO's gene-level bias, not low module quality.
  *(Corrected 2026-09-12 from a stale "4/4"; see `GO_INVISIBLE_GATE_SUMMARY.md`.)*
- Abundance and switch axes are near-orthogonal in every analysis, not only the SCZD store
  (per-analysis median |r| and fraction |r| < 0.1 in `_m/AXIS_ORTHOGONALITY.md`); 34 SCZD / 43
  aging-caudate composition-unique genes carry phenotype signal total abundance misses.

**CLIs:** `isograph_benchmark/real_data/{interpret_modules,module_enrichment,go_enrichment,go_invisible_gate,incremental_association,abundance_structure_separation,characterize_composition_unique,celltype_composition}.py`.

## Display items

S-real-3 `figGoInvisible`, S-real-4 `figSeparation`, **Fig 5** `figCompositionRobustness`;
supplementary table S6 (GO-invisible modules), S13 (composition). S-real-1 `figBaselineRates` and
tables S1/S2 come from the baseline comparison in stage 04.
