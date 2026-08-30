# Clinical-consequence of switched exons — gtex_frontal_cortex_ba9

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **139640** (switched 130586, background 9054; CDS-overlapping 114948). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 6555 | 2.805 | 15.833 | 0.18 | 0.002 | 2258 |
| all_exons | go_invisible | 2942 | 2.647 | 19.312 | 0.14 | 0.002 | 1059 |
| all_exons | go_visible | 3613 | 2.936 | 12.056 | 0.24 | 0.0819 | 1199 |
| cds | all | 5925 | 3.330 | 16.243 | 0.21 | 0.002 | 2141 |
| cds | go_invisible | 2493 | 3.291 | 19.899 | 0.17 | 0.002 | 986 |
| cds | go_visible | 3432 | 3.360 | 12.307 | 0.27 | 0.162 | 1155 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 5685 | 0.760 | 0.936 | 1.05e-151 |
| go_invisible | 2401 | 0.775 | 0.936 | 3.74e-64 |
| go_visible | 3284 | 0.748 | 0.936 | 1.4e-100 |
