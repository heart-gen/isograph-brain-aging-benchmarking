# Clinical-consequence of switched exons — gtex_anterior_cingulate_cortex_ba24

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **75899** (switched 71159, background 4740; CDS-overlapping 62081). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 3333 | 2.716 | 14.930 | 0.18 | 0.00599 | 1196 |
| all_exons | go_invisible | 1533 | 2.668 | 15.131 | 0.18 | 0.002 | 601 |
| all_exons | go_visible | 1800 | 2.758 | 14.734 | 0.19 | 0.509 | 595 |
| cds | all | 3006 | 3.240 | 15.383 | 0.21 | 0.016 | 1130 |
| cds | go_invisible | 1381 | 3.162 | 15.549 | 0.20 | 0.002 | 564 |
| cds | go_visible | 1625 | 3.307 | 15.222 | 0.22 | 0.695 | 566 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2871 | 0.688 | 0.936 | 1.26e-147 |
| go_invisible | 1331 | 0.661 | 0.936 | 1.06e-89 |
| go_visible | 1540 | 0.708 | 0.936 | 1.48e-67 |
