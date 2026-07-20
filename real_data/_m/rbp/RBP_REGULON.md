# Candidate RBP regulons among IsoGraph co-switch modules (mature scope)

Motifs are scanned on the **mature transcript** (exonic + UTR) sequence.

For each module and RBP, whether RBP binding-site switching (motif gained/lost between the switch-pair isoforms) is over-represented among the module's genes vs the region's switch-gene pool (hypergeometric, BH across module x RBP tests). A significant RBP marks the module as a candidate regulon for it.

- switch genes tested: **11859** over 14 regions
- module x RBP tests: **29204**; candidate regulons at q<0.05: **824** (GO-invisible 177)

| region | module | RBP | module size | switched | enrichment | q | GO-inv |
|--------|--------|-----|-------------|----------|------------|---|--------|
| frontal_cortex_ba9 | M008 | PPIE | 273 | 142 | 3.22 | 2.40e-40 | no |
| frontal_cortex_ba9 | M008 | ELAVL4 | 273 | 96 | 2.99 | 4.33e-21 | no |
| frontal_cortex_ba9 | M008 | KHDRBS1 | 273 | 158 | 1.98 | 7.37e-20 | no |
| frontal_cortex_ba9 | M008 | HNRNPD | 273 | 144 | 2.10 | 8.79e-20 | no |
| frontal_cortex_ba9 | M007 | A1CF | 406 | 271 | 1.55 | 1.99e-19 | no |
| frontal_cortex_ba9 | M007 | RBMS3 | 406 | 245 | 1.62 | 9.87e-19 | no |
| frontal_cortex_ba9 | M008 | SYNCRIP | 273 | 144 | 2.02 | 3.29e-18 | no |
| frontal_cortex_ba9 | M008 | RNASEL | 273 | 115 | 2.29 | 4.04e-17 | no |
| frontal_cortex_ba9 | M007 | HNRNPCL1 | 406 | 304 | 1.40 | 2.51e-16 | no |
| frontal_cortex_ba9 | M007 | CPEB2 | 406 | 300 | 1.40 | 3.50e-16 | no |
| frontal_cortex_ba9 | M008 | ELAVL3 | 273 | 163 | 1.74 | 6.01e-15 | no |
| frontal_cortex_ba9 | M007 | RBM41 | 406 | 195 | 1.68 | 1.63e-14 | no |
| frontal_cortex_ba9 | M007 | RBMS1 | 406 | 162 | 1.82 | 5.49e-14 | no |
| frontal_cortex_ba9 | M008 | OAS1 | 273 | 73 | 2.81 | 7.86e-14 | no |
| frontal_cortex_ba9 | M007 | ADAR | 406 | 258 | 1.45 | 1.92e-13 | no |
| frontal_cortex_ba9 | M007 | RALY | 406 | 285 | 1.39 | 2.38e-13 | no |
| frontal_cortex_ba9 | M007 | DDX58 | 406 | 223 | 1.53 | 6.34e-13 | no |
| anterior_cingulate_cortex_ba24 | M023 | PPIE | 77 | 41 | 3.65 | 1.33e-12 | no |
| frontal_cortex_ba9 | M007 | DHX9 | 406 | 237 | 1.48 | 1.56e-12 | no |
| frontal_cortex_ba9 | M008 | IGF2BP3 | 273 | 98 | 2.15 | 4.58e-12 | no |
| frontal_cortex_ba9 | M001 | ZCRB1 | 677 | 256 | 1.51 | 5.40e-12 | no |
| frontal_cortex_ba9 | M008 | TIAL1 | 273 | 52 | 3.29 | 6.84e-12 | no |
| frontal_cortex_ba9 | M008 | U2AF2 | 273 | 175 | 1.55 | 1.03e-11 | no |
| frontal_cortex_ba9 | M007 | ZFP36L2 | 406 | 293 | 1.33 | 1.32e-11 | no |
| frontal_cortex_ba9 | M001 | RBM41 | 677 | 281 | 1.45 | 1.57e-11 | no |
| frontal_cortex_ba9 | M007 | ZNF346 | 406 | 161 | 1.70 | 3.18e-11 | no |
| frontal_cortex_ba9 | M001 | RALY | 677 | 435 | 1.27 | 4.09e-11 | no |
| frontal_cortex_ba9 | M001 | HNRNPCL1 | 677 | 454 | 1.25 | 4.89e-11 | no |
| frontal_cortex_ba9 | M008 | ELAVL2 | 273 | 50 | 3.22 | 4.89e-11 | no |
| frontal_cortex_ba9 | M007 | PABPC4 | 406 | 289 | 1.32 | 1.00e-10 | no |
