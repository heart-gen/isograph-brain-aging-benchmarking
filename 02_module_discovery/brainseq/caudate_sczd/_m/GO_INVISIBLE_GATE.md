# Biology gate — GO-invisible switch modules (brainseq-sczd)

Phenotype-significant IsoGraph switch modules (pheno_fdr <= 0.1): **8** — 2 GO-enriched, **6 GO-invisible** (M026, M020, M010, M022, M023, M025).

Reproduce: `python -m isograph_benchmark.real_data.go_invisible_gate --analysis brainseq-sczd`

## Per-module switch coherence

| module | go_invisible | n_genes | pheno_fdr | genes_with_real_switch | max_switch_strength | n_sig_switch_tx | top_switch_genes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M012 | False | 58 | 0.0001 | 40 | 1.234 | 225 | OSMR-DT, SOCS3, LUCAT1, IL6, SBNO2, IRF1, CSF3, GPRC5A |
| M011 | False | 65 | 0.0448 | 59 | 1.55 | 450 | BTK, CCR5AS, PTPN6, ENSG00000250604, PCED1B-AS1, SAMSN1, IRF8, LNCAROD |
| M026 | True | 23 | 0.0019 | 23 | 1.662 | 185 | UNC119, DDX27, HEXA, ZNG1A, UBE2B, PAPOLA, NBN, CEP290 |
| M020 | True | 30 | 0.0172 | 4 | 1.021 | 7 | RAB14, ZBTB41, RLIM |
| M010 | True | 83 | 0.0272 | 65 | 1.312 | 509 | SH3TC2-DT, SEPTIN8, DMBT1, SPART, LINC02882, ENSG00000291176, FAM107B, ENSG00000286856 |
| M022 | True | 29 | 0.033 | 29 | 1.549 | 200 | BBIP1, TERF1, TTC21B, ACER3, DDX5, PLCG1, RAF1, CEP57 |
| M023 | True | 27 | 0.0471 | 27 | 1.615 | 158 | ARNT2, PBX1, BCR, SLC25A12, SEM1, ST3GAL2, DGKH, SLC25A46 |
| M025 | True | 24 | 0.071 | 11 | 1.719 | 47 | ZC3H7B, MOB3B, CCDC6, ATP5PD, NDUFA13, ANKRD33B, HS3ST5, HSP90AA1 |

Pooled background functional-consequence fractions: cds_changed 0.844, coding_status_change 0.677, biotype_switch 0.753, utr_changed 0.601.

## Reading

- `genes_with_real_switch` near `n_genes` = members carry genuine DTU (both up- and down-correlated transcripts), not abundance shifts.
- Functional-consequence fractions at/above background and indistinguishable from GO-visible disease modules => GO-invisible modules are not lower quality; GO-invisibility reflects GO's gene-level/abundance bias, blind to DTU.
- `top_switch_genes` heterogeneous within a module (shared switch axis, not a shared GO process) => frame as a complementary DTU-without-DGE layer, NOT pathway discovery WGCNA misses.
