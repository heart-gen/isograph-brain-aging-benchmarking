# Biology gate — GO-invisible switch modules (brainseq-sczd)

Phenotype-significant IsoGraph switch modules (pheno_fdr <= 0.1): **4** — 0 GO-enriched, **4 GO-invisible** (M026, M020, M010, M023).

Reproduce: `python -m isograph_benchmark.real_data.go_invisible_gate --analysis brainseq-sczd`

## Per-module switch coherence

| module | go_invisible | n_genes | pheno_fdr | genes_with_real_switch | max_switch_strength | n_sig_switch_tx | top_switch_genes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M026 | True | 23 | 0.0019 | 21 | 1.098 | 93 | ABCB6, JAM3, XPO1, DDX3X, TEX2, DBNDD2, ENPP2, IFNAR2 |
| M020 | True | 30 | 0.0172 | 31 | 1.392 | 148 | SNF8, UBE2Z, SFXN1, STXBP5, ZDHHC17, PAM16, SPTBN1, CANX |
| M010 | True | 83 | 0.0272 | 73 | 1.251 | 459 | SEPTIN8, FAM107B, SPART, HAGLR, ENSG00000291176, SOX2-OT, PHLDB1, SH3TC2-DT |
| M023 | True | 27 | 0.0471 | 28 | 1.309 | 121 | DHX36, ARNT2, SLC25A12, SLC25A46, FAM234B, PBX1, CAPZA2, DGKH |

Pooled background functional-consequence fractions: cds_changed 0.84, coding_status_change 0.671, biotype_switch 0.736, utr_changed 0.612.

## Reading

- `genes_with_real_switch` near `n_genes` = members carry genuine DTU (both up- and down-correlated transcripts), not abundance shifts.
- Functional-consequence fractions at/above background and indistinguishable from GO-visible disease modules => GO-invisible modules are not lower quality; GO-invisibility reflects GO's gene-level/abundance bias, blind to DTU.
- `top_switch_genes` heterogeneous within a module (shared switch axis, not a shared GO process) => frame as a complementary DTU-without-DGE layer, NOT pathway discovery WGCNA misses.
