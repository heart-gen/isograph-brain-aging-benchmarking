# Clinical-consequence of switched exons — gtex_hippocampus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **59731** (switched 55926, background 3805; CDS-overlapping 46368). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 2709 | 2.544 | 13.920 | 0.18 | 0.002 | 961 |
| all_exons | go_invisible | 1050 | 2.929 | 22.359 | 0.13 | 0.00799 | 399 |
| all_exons | go_visible | 1659 | 2.305 | 8.771 | 0.26 | 0.002 | 562 |
| cds | all | 2266 | 3.232 | 14.278 | 0.23 | 0.002 | 902 |
| cds | go_invisible | 902 | 3.574 | 22.802 | 0.16 | 0.00599 | 381 |
| cds | go_visible | 1364 | 3.005 | 9.025 | 0.33 | 0.002 | 521 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2172 | 0.732 | 0.936 | 1.45e-82 |
| go_invisible | 863 | 0.779 | 0.936 | 8.14e-23 |
| go_visible | 1309 | 0.697 | 0.936 | 2.99e-66 |
