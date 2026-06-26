# Biology gate — GO-invisible switch modules (brainseq-sczd)

Phenotype-significant IsoGraph switch modules (pheno_fdr <= 0.1): **8** — 2 GO-enriched, **6 GO-invisible** (M010, M017, M020, M023, M021, M026).

Reproduce: `python -m isograph_benchmark.real_data.go_invisible_gate --analysis brainseq-sczd`

## Per-module switch coherence

| module | go_invisible | n_genes | pheno_fdr | genes_with_real_switch | max_switch_strength | n_sig_switch_tx | top_switch_genes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M015 | False | 36 | 0.0 | 27 | 1.207 | 129 | IL6, OSMR-DT, CSF3, SBNO2, GPRC5A, IRF1, FSTL3, SECTM1 |
| M013 | False | 50 | 0.0386 | 44 | 1.412 | 350 | BTK, CCR5AS, PTPN6, ENSG00000250604, PCED1B-AS1, IRF8, SAMSN1, LGALS9 |
| M010 | True | 90 | 0.0386 | 73 | 1.251 | 459 | SEPTIN8, FAM107B, SPART, HAGLR, ENSG00000291176, SOX2-OT, PHLDB1, SH3TC2-DT |
| M017 | True | 33 | 0.0386 | 10 | 1.01 | 17 | ETV1, ZNF805, RLIM, NUP58, PHLDA1, GNA13, ZBTB41, C4orf3 |
| M020 | True | 31 | 0.0386 | 31 | 1.392 | 148 | SNF8, UBE2Z, SFXN1, STXBP5, ZDHHC17, PAM16, SPTBN1, CANX |
| M023 | True | 28 | 0.0386 | 28 | 1.309 | 121 | DHX36, ARNT2, SLC25A12, SLC25A46, FAM234B, PBX1, CAPZA2, DGKH |
| M021 | True | 31 | 0.0413 | 11 | 1.319 | 59 | CLASP2, GNAL, LNX1, PRKCB, FRMPD4, ZNF584, KCNAB1, UBC |
| M026 | True | 21 | 0.0999 | 21 | 1.098 | 93 | ABCB6, JAM3, XPO1, DDX3X, TEX2, DBNDD2, ENPP2, IFNAR2 |

Pooled background functional-consequence fractions: cds_changed 0.84, coding_status_change 0.671, biotype_switch 0.736, utr_changed 0.612.

## Reading

- `genes_with_real_switch` near `n_genes` = members carry genuine DTU (both up- and down-correlated transcripts), not abundance shifts.
- Functional-consequence fractions at/above background and indistinguishable from GO-visible disease modules => GO-invisible modules are not lower quality; GO-invisibility reflects GO's gene-level/abundance bias, blind to DTU.
- `top_switch_genes` heterogeneous within a module (shared switch axis, not a shared GO process) => frame as a complementary DTU-without-DGE layer, NOT pathway discovery WGCNA misses.
