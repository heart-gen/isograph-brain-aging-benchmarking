# Candidate RBP regulons among IsoGraph co-switch modules (intronic scope)

Motifs are scanned on **intronic splice-site flanks** (pre-mRNA sense; up to 100 nt into each intron), the binding niche for splicing-regulatory RBPs invisible to the mature-transcript scan.

For each module and RBP, whether RBP binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool (hypergeometric, BH across module x RBP tests). A significant RBP marks the module as a candidate regulon for it.

- switch genes tested: **11859** over 14 regions
- module x RBP tests: **33431**; candidate regulons at q<0.05: **684** (GO-invisible 265)

| region | module | RBP | module size | switched | enrichment | q | GO-inv |
|--------|--------|-----|-------------|----------|------------|---|--------|
| frontal_cortex_ba9 | M001 | RBMS3 | 677 | 250 | 1.79 | 5.09e-21 | no |
| frontal_cortex_ba9 | M001 | ZFP36L2 | 677 | 358 | 1.47 | 5.29e-17 | no |
| frontal_cortex_ba9 | M001 | ZCRB1 | 677 | 170 | 1.93 | 1.40e-15 | no |
| frontal_cortex_ba9 | M001 | AKAP1 | 677 | 326 | 1.47 | 6.79e-15 | no |
| frontal_cortex_ba9 | M008 | PPIE | 273 | 132 | 1.99 | 6.79e-15 | no |
| frontal_cortex_ba9 | M007 | RALY | 406 | 229 | 1.56 | 3.32e-14 | no |
| frontal_cortex_ba9 | M001 | PABPC5 | 677 | 243 | 1.60 | 8.31e-14 | no |
| frontal_cortex_ba9 | M001 | HNRNPCL1 | 677 | 351 | 1.41 | 1.02e-13 | no |
| frontal_cortex_ba9 | M008 | ELAVL4 | 273 | 118 | 2.03 | 2.28e-13 | no |
| frontal_cortex_ba9 | M001 | SNRPB2 | 677 | 362 | 1.37 | 3.30e-12 | no |
| frontal_cortex_ba9 | M001 | RBM41 | 677 | 172 | 1.74 | 8.81e-12 | no |
| frontal_cortex_ba9 | M003 | TIA1 | 514 | 98 | 2.22 | 1.13e-11 | yes |
| frontal_cortex_ba9 | M001 | RALY | 677 | 338 | 1.38 | 1.23e-11 | no |
| frontal_cortex_ba9 | M007 | ZFP36L2 | 406 | 220 | 1.50 | 2.92e-11 | no |
| frontal_cortex_ba9 | M001 | CPEB2 | 677 | 357 | 1.35 | 3.64e-11 | no |
| frontal_cortex_ba9 | M003 | PTBP1 | 514 | 77 | 2.44 | 6.84e-11 | yes |
| frontal_cortex_ba9 | M003 | HNRNPA1 | 514 | 82 | 2.36 | 6.84e-11 | yes |
| frontal_cortex_ba9 | M003 | SRSF1 | 514 | 83 | 2.32 | 9.90e-11 | yes |
| frontal_cortex_ba9 | M003 | NOVA1 | 514 | 77 | 2.40 | 1.40e-10 | yes |
| frontal_cortex_ba9 | M001 | MATR3 | 677 | 145 | 1.79 | 1.68e-10 | no |
| frontal_cortex_ba9 | M001 | PABPC4 | 677 | 278 | 1.43 | 2.31e-10 | no |
| frontal_cortex_ba9 | M003 | PCBP2 | 514 | 83 | 2.28 | 2.53e-10 | yes |
| frontal_cortex_ba9 | M008 | CELF2 | 273 | 100 | 2.00 | 2.80e-10 | no |
| frontal_cortex_ba9 | M001 | SART3 | 677 | 253 | 1.47 | 4.76e-10 | no |
| frontal_cortex_ba9 | M007 | RBMS3 | 406 | 144 | 1.72 | 6.23e-10 | no |
| frontal_cortex_ba9 | M001 | A1CF | 677 | 237 | 1.48 | 1.41e-09 | no |
| frontal_cortex_ba9 | M003 | ZFP36 | 514 | 85 | 2.17 | 1.52e-09 | yes |
| frontal_cortex_ba9 | M003 | SSB | 514 | 135 | 1.77 | 1.63e-09 | yes |
| frontal_cortex_ba9 | M003 | CELF2 | 514 | 157 | 1.67 | 1.72e-09 | yes |
| frontal_cortex_ba9 | M008 | ELAVL2 | 273 | 102 | 1.91 | 2.47e-09 | no |
