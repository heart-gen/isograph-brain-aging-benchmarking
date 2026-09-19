# Clinical-consequence of switched exons — brainseq_caudate_sczd

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **36652** (switched 34589, background 2063; CDS-overlapping 29202). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1788 | 2.389 | 27.585 | 0.09 | 0.002 | 508 |
| all_exons | go_invisible | 1437 | 2.321 | 29.112 | 0.08 | 0.002 | 409 |
| all_exons | go_visible | 351 | 2.835 | 12.713 | 0.22 | 0.002 | 99 |
| cds | all | 1565 | 2.934 | 28.532 | 0.10 | 0.002 | 471 |
| cds | go_invisible | 1241 | 2.851 | 30.040 | 0.09 | 0.002 | 381 |
| cds | go_visible | 324 | 3.467 | 13.459 | 0.26 | 0.002 | 90 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 1512 | 0.841 | 0.936 | 4.05e-19 |
| go_invisible | 1201 | 0.839 | 0.936 | 2.45e-16 |
| go_visible | 311 | 0.854 | 0.936 | 0.000102 |
