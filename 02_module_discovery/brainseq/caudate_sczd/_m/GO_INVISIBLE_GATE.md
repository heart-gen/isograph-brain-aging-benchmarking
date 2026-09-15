# Biology gate — GO-invisible switch modules (brainseq-sczd)

Phenotype-significant IsoGraph switch modules (pheno_fdr <= 0.1): **10** — 1 GO-enriched, **9 GO-invisible** (M010, M017, M015, M013, M031, M030, M012, M018, M014).

Reproduce: `python -m isograph_benchmark.real_data.go_invisible_gate --analysis brainseq-sczd`

## Per-module switch coherence

| module | go_invisible | n_genes | pheno_fdr | genes_with_real_switch | max_switch_strength | n_sig_switch_tx | top_switch_genes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M020 | False | 43 | 0.0065 | 33 | 1.72 | 144 | NDUFS8, PPME1, AP2M1, ZNHIT1, NDUFB9, CNOT2, BLOC1S1, FH |
| M010 | True | 145 | 0.0002 | 126 | 1.226 | 277 | LINC03000, ENSG00000293389, PALS2, MYEF2, NDRG4, FRG1HP, PASK, ABL2 |
| M017 | True | 59 | 0.0017 | 57 | 1.558 | 228 | CWC27, DCAF6, GDAP2, BAZ2B, COQ9, PPP6R1, ATXN2, MRPS33 |
| M015 | True | 66 | 0.0066 | 65 | 1.646 | 239 | BCCIP, SYT7, RPL15, H2BC6, SC5D, GSTM2, GTF3C1, SLC38A1 |
| M013 | True | 105 | 0.0113 | 103 | 1.718 | 485 | ATP5IF1, ERCC6L2, FLCN, H2BC26, ACTR1B, TM2D3, KCTD13, CDK10 |
| M031 | True | 21 | 0.0113 | 21 | 1.389 | 75 | PPP1R12C, DNAJC8, PDAP1, USP4, SLC6A8, FAM98C, OSBPL1A, C1orf174 |
| M030 | True | 22 | 0.0364 | 22 | 1.423 | 60 | LAMP1, RNF13, NRXN2, SYNGR2, TUBA1A, PCCB, CLDND1, LPAR1 |
| M012 | True | 113 | 0.0712 | 98 | 1.748 | 423 | NAPA, REEP1, LRP11, ANKLE2, UBE3C, CABP1, MACROD2, JPH4 |
| M018 | True | 53 | 0.0712 | 50 | 1.779 | 326 | ENSA, NEDD8, CYTH2, NDUFV1, MSANTD1, AVL9, UQCRC1, PRMT7 |
| M014 | True | 69 | 0.082 | 69 | 1.652 | 378 | CISD1, NDUFA3, EIF3L, UBB, MEIS2, DDX18, CDIPT, RAN |

Pooled background functional-consequence fractions: cds_changed 0.929, coding_status_change 0.761, biotype_switch 0.847, utr_changed 0.606.

## Reading

- `genes_with_real_switch` near `n_genes` = members carry genuine DTU (both up- and down-correlated transcripts), not abundance shifts.
- Functional-consequence fractions at/above background and indistinguishable from GO-visible disease modules => GO-invisible modules are not lower quality; GO-invisibility reflects GO's gene-level/abundance bias, blind to DTU.
- `top_switch_genes` heterogeneous within a module (shared switch axis, not a shared GO process) => frame as a complementary DTU-without-DGE layer, NOT pathway discovery WGCNA misses.
