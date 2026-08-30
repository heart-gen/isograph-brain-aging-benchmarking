# Synthetic sample-table refresh

2290 unique datasets in the run grid. Only `samples.parquet` is rewritten; the regenerated count matrices, feature tables and truth tables are fingerprinted against the cached copies and then discarded, so a dataset that fails to reproduce is reported rather than overwritten.

| status | n | meaning |
|---|---|---|
| `mismatch` | 1204 | **a shared artifact did not reproduce; left untouched** |
| `unchanged` | 560 | sample table already matched the generator |
| `refreshed` | 526 | sample table rewritten; all other artifacts reproduced exactly |

Datasets whose cached sample table was missing at least one of `RIN`, `neuron_frac`, `batch`, `library_size` before this pass: **1730**. Those are the datasets on which `isograph_vae_residual` was silently reducing to `isograph_vae`.

## Datasets that did not reproduce

These block the refresh. The cached expression data is not recoverable from the committed config, which invalidates the assumption behind every cached benchmark result -- escalate rather than re-running.

| dataset_id | scenario | seed | n mismatched | artifacts |
|---|---|---|---|---|
| `026327dfc26eb89c` | feature_space_interactions | 13 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `79ae2de5a7a13829` | feature_space_interactions | 13 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `2a4f315f3fdb58cc` | feature_space_interactions | 14 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `a21933aa6240a746` | feature_space_interactions | 14 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `2154e3ceaed45bba` | feature_space_interactions | 15 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `e211bc27e132a2c8` | feature_space_interactions | 15 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `234a471b8b399378` | feature_space_interactions | 16 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `9e0284a6ffaca782` | feature_space_interactions | 16 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `378676446dcea9e0` | feature_space_interactions | 17 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `ed940ebf582bc267` | feature_space_interactions | 17 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `3d2cf9fc3205f557` | feature_space_interactions | 18 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `df0ec80275beb60a` | feature_space_interactions | 18 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `692c9e3064a825d6` | feature_space_interactions | 19 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `f30d2b5ca53ba363` | feature_space_interactions | 19 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `3230ae30a69a65a4` | feature_space_interactions | 20 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `245392f0ff993d69` | feature_space_interactions | 21 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `a2e1fbdda5fdeb6c` | feature_space_interactions | 21 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `171b121bf76b4c65` | feature_space_interactions | 22 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `61ff512700b9d679` | feature_space_interactions | 22 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `304e8b5b13492da7` | feature_space_interactions | 23 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `ae521e4a6a3585c1` | feature_space_interactions | 23 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `3602531e6146e5b2` | feature_space_interactions | 24 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `b84b44da94fa93d7` | feature_space_interactions | 24 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `3a9c8648989f0756` | feature_space_interactions | 25 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `b733f61b7650f8f3` | feature_space_interactions | 25 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `61b60a2a9a8b2b52` | feature_space_interactions | 26 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `fdca75eb2138e87c` | feature_space_interactions | 26 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `75f106276bcbc312` | feature_space_interactions | 27 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `e52ce80a95562f4f` | feature_space_interactions | 27 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `8438dcfc5cb62eb0` | feature_space_interactions | 28 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `0a485b041b7ec6e9` | feature_space_interactions | 29 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `9626df97ac389642` | feature_space_interactions | 29 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `11771b1677ab2a84` | feature_space_interactions | 30 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `bb12b3fcd42a1641` | feature_space_interactions | 30 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `190b1565801e9b45` | feature_space_interactions | 31 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `a62cd26a73e3f5aa` | feature_space_interactions | 31 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `211424736b8af09e` | feature_space_interactions | 32 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `ae0decf98486bbbb` | feature_space_interactions | 32 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `3ca228687055d16b` | idealized_switching | 13 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `07a0243b7798762b` | idealized_switching | 14 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `87f3643ede55ee9b` | idealized_switching | 14 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `5c0e73098f9a4e91` | idealized_switching | 15 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `050851351bae7001` | idealized_switching | 16 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `8eea73e0b82c186e` | idealized_switching | 16 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `5faa8f892e3d41aa` | idealized_switching | 17 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `4752d26d52bfb706` | idealized_switching | 18 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `e134d9b220db2d8b` | idealized_switching | 18 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `5f7b3f59a23fedfa` | idealized_switching | 19 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `00903940206ecb93` | idealized_switching | 20 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |
| `c7e4c64547b3943b` | idealized_switching | 20 | 6 | `f:transcript,f:truth_module,m:gene_counts,m:psi,m:transcript_counts,t:truth_modules.parquet` |

**Generator additions:** 1810 dataset(s) predate one or more artifacts the generator now emits (for example `truth_switch_event`). Nothing cached differs; those artifacts are simply absent. This tool writes only `samples.parquet` and does not backfill them -- that is a separate decision, since downstream metrics may or may not require them.
