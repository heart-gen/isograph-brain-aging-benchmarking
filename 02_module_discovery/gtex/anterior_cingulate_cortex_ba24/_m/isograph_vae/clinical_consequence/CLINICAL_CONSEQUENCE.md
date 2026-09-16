# Clinical-consequence of switched exons — gtex_anterior_cingulate_cortex_ba24

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **63540** (switched 59336, background 4204; CDS-overlapping 53516). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 3020 | 2.224 | 15.814 | 0.14 | 0.002 | 1048 |
| all_exons | go_invisible | 467 | 1.965 | 10.618 | 0.19 | 0.002 | 128 |
| all_exons | go_visible | 2553 | 2.267 | 16.464 | 0.14 | 0.002 | 920 |
| cds | all | 2833 | 2.552 | 16.245 | 0.16 | 0.002 | 1003 |
| cds | go_invisible | 394 | 2.480 | 12.149 | 0.20 | 0.002 | 114 |
| cds | go_visible | 2439 | 2.563 | 16.699 | 0.15 | 0.002 | 889 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2729 | 0.720 | 0.936 | 7.57e-110 |
| go_invisible | 384 | 0.911 | 0.936 | 0.0273 |
| go_visible | 2345 | 0.691 | 0.936 | 6.79e-121 |
