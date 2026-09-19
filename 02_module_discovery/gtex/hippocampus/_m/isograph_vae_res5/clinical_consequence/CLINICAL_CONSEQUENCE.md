# Clinical-consequence of switched exons — gtex_hippocampus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **31861** (switched 29860, background 2001; CDS-overlapping 26774). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1549 | 1.964 | 7.304 | 0.27 | 0.002 | 546 |
| all_exons | go_invisible | 40 | 1.563 | 14.031 | 0.11 | 0.002 | 10 |
| all_exons | go_visible | 1509 | 1.977 | 7.192 | 0.27 | 0.002 | 536 |
| cds | all | 1453 | 2.265 | 7.558 | 0.30 | 0.002 | 517 |
| cds | go_invisible | 38 | 1.746 | 14.341 | 0.12 | 0.002 | 9 |
| cds | go_visible | 1415 | 2.281 | 7.444 | 0.31 | 0.002 | 508 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 1391 | 0.699 | 0.936 | 2.03e-66 |
| go_invisible | 38 | 0.695 | 0.936 | 0.000146 |
| go_visible | 1353 | 0.699 | 0.936 | 7.7e-64 |
