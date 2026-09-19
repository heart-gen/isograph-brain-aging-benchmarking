# Clinical-consequence of switched exons — gtex_frontal_cortex_ba9

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **109241** (switched 102029, background 7212; CDS-overlapping 90889). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 4996 | 2.560 | 14.857 | 0.17 | 0.002 | 1729 |
| all_exons | go_invisible | 678 | 4.363 | 15.503 | 0.28 | 0.002 | 262 |
| all_exons | go_visible | 4318 | 2.288 | 14.742 | 0.16 | 0.002 | 1467 |
| cds | all | 4599 | 2.996 | 15.239 | 0.20 | 0.002 | 1650 |
| cds | go_invisible | 569 | 5.614 | 15.912 | 0.35 | 0.002 | 243 |
| cds | go_visible | 4030 | 2.639 | 15.119 | 0.17 | 0.002 | 1407 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 4409 | 0.721 | 0.936 | 1.44e-172 |
| go_invisible | 547 | 0.765 | 0.936 | 1.62e-20 |
| go_visible | 3862 | 0.714 | 0.936 | 2.76e-158 |
