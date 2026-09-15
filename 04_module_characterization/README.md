# 04 — Module characterization

**Question:** what *are* these modules, and do they carry information abundance
pipelines miss?

This stage produces the honest scope bound for the whole paper: on per-module rates
IsoGraph is **not** globally superior to WGCNA. Phenotype signal lives in the switch
features (both switch-fed methods win) and GO enrichment is abundance-dominated. The one
clean method effect is `isograph` > `wgcna_multiplex` on identical features. Read every
result here against the three matched baselines from stage 02.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–02 | `interpret_modules_{brainseq,gtex}` | Switch drivers, transcript polarity, high-vs-low contrasts, GENCODE v47 structural annotation |
| 03–04 | `module_enrichment_{brainseq,gtex}` | Per-module GO:BP + network metrics joined to phenotype FDR, all methods |
| 05 | `go_invisible_gate` | The GO-invisible disease-module gate (BrainSEQ SCZD) |
| 06–07 | `incremental_association_{brainseq,gtex}` | De-confounded test: does the switch channel add signal *beyond* gene-level abundance? |
| 08–09 | `abundance_structure`, `characterize_composition_unique` | Abundance/switch orthogonality; composition-unique gene sets |
| 10–11 | `celltype_composition_brainseq`, `gtex_composition` | Cell-type composition (MuSiC deconvolution) |
| 12–13 | `project_tiers_pilot`, `tiers_fanout` | Tier projection — **retired 2026-09-14**: cited nowhere; not re-run on the switching filter, outputs only at tag `legacy_expression_filter` |
| 14 | `baseline_comparison` | Three-baseline per-module rate comparison (the scope bound) |

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
- Abundance and switch axes are near-orthogonal (median |r| ≈ 0.13); 34 SCZD / 43
  aging-caudate composition-unique genes carry phenotype signal total abundance misses.

**CLIs:** `isograph_benchmark/real_data/{interpret_modules,module_enrichment,go_enrichment,go_invisible_gate,incremental_association,abundance_structure_separation,characterize_composition_unique,celltype_composition,project_tiers,baseline_comparison}.py`.

## Display items

S-real-1 `figBaselineRates`, S-real-3 `figGoInvisible`, S-real-4 `figSeparation`;
supplementary tables S1/S2 (baseline rates), S6 (GO-invisible modules).
