# Module trust: display tables

Written by `isograph_benchmark/real_data/module_trust_tables.py`. Presentation only -- every column is copied from the ledger named beside it; regenerate rather than editing.

| file | source ledger |
|---|---|
| `split_half_modules.csv` | `module_stability__<cohort>__<region>__<method>.parquet` |
| `split_half_pairs.csv` | `within_cohort__<cohort>__<region>__<method>.parquet` |
| `projection_modules.csv` | `eigengene_projection/eigengene_projection_all.parquet` |
| `projection_summary.csv` | `eigengene_projection/eigengene_projection_summary.parquet` |
| `functional_preservation.csv` | `functional_preservation__<method>__<model>__stats.json` |
| `crosscohort_permutation.csv` | `replication_permutation__*__stats.json` |
| `driver_structure.csv` | `module_complementarity__<cohort>__<region>__<method>.parquet` |
| `resolution_sensitivity.csv` | `stability_summary.parquet` |
| `projection_sign_scale.csv` | `eigengene_projection/eigengene_projection_all.parquet + eigengene_projection/switch_axis_alignment__<region>.parquet` |

## The four claims

| claim | IsoGraph | WGCNA |
|---|---|---|
| modules_chance_trusted | 93/118 (0.79) | 62/72 (0.86) |
| split_half_age_sign_concordance | 22/23 (0.96) | 14/14 (1.00) |
| aging_axis_transfers__brainseq_to_gtex | 33/38 (0.87) | 28/52 (0.54) |
| aging_axis_transfers__gtex_to_brainseq | 44/55 (0.80) | 3/10 (0.30) |

## Functional preservation of matched pairs (`linear` arm)

| method | measure | n | matched mean | size-matched null | p_emp |
|---|---|---|---|---|---|
| isograph | `celltype_r` | 38 | 0.2081 | 0.0396 | 0.00599 |
| isograph | `go_jaccard` | 24 | 0.0971 | 0.0296 | 0.000999 |
| isograph | `structure_r` | 0 | untested (no computable pair) | n/a | n/a |
| wgcna | `celltype_r` | 0 | untested (no computable pair) | n/a | n/a |
| wgcna | `go_jaccard` | 52 | 0.1356 | 0.1393 | 1 |
| wgcna | `structure_r` | 0 | untested (no computable pair) | n/a | n/a |

## Matched-pair count: the two arms disagree

- **isograph**, `pearson` vs the matching null: 1/38, p_emp = 0.912 *(quoted in the text)*
- **isograph**, `spline_f` vs the matching null: 0/38, p_emp = 1 *(stage's pre-registered primary)*
- **wgcna**, `pearson` vs the matching null: 2/52, p_emp = 0.151 *(quoted in the text)*
- **wgcna**, `spline_f` vs the matching null: 5/52, p_emp = 0.161 *(stage's pre-registered primary)*

## Why the projected-age statistic is null-standardised

Raw sign agreement is counted among the modules significant in both cohorts on raw correlations; the standardised column counts every age-testable module. A region that agrees unanimously in one direction and a region that disagrees unanimously are the same phenomenon: the target cohort's own age-correlated structure enters every projection.

| method | direction | region | raw (both-sig) | standardised (all) |
|---|---|---|---|---|
| isograph | brainseq_to_gtex | caudate | 8/8 (1.00) | 11/12 (0.92) |
| isograph | brainseq_to_gtex | dlpfc_ba9 | 8/8 (1.00) | 9/10 (0.90) |
| isograph | brainseq_to_gtex | hippocampus | 0/13 (0.00) | 13/16 (0.81) |
| isograph | gtex_to_brainseq | caudate | 1/1 (1.00) | 20/22 (0.91) |
| isograph | gtex_to_brainseq | dlpfc_ba9 | 2/2 (1.00) | 11/12 (0.92) |
| isograph | gtex_to_brainseq | hippocampus | 0/14 (0.00) | 13/21 (0.62) |
| wgcna | brainseq_to_gtex | caudate | 1/1 (1.00) | 14/23 (0.61) |
| wgcna | brainseq_to_gtex | dlpfc_ba9 | 3/3 (1.00) | 12/24 (0.50) |
| wgcna | brainseq_to_gtex | hippocampus | 1/4 (0.25) | 2/5 (0.40) |
| wgcna | gtex_to_brainseq | caudate | 1/1 (1.00) | 1/6 (0.17) |
| wgcna | gtex_to_brainseq | dlpfc_ba9 | 1/1 (1.00) | 1/3 (0.33) |
| wgcna | gtex_to_brainseq | hippocampus | 0/1 (0.00) | 1/1 (1.00) |

Report both. The pre-registered rule in `REPLICATION_PERMUTATION.md` forbids the word *replication* for this arm; the honest wording is "matched modules with concordant age effects".
