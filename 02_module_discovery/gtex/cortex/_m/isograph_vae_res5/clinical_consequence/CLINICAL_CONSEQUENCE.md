# Clinical-consequence of switched exons — gtex_cortex

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **62986** (switched 58644, background 4342; CDS-overlapping 53707). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 2714 | 2.730 | 20.278 | 0.13 | 0.002 | 952 |
| all_exons | go_invisible | 549 | 4.046 | 26.626 | 0.15 | 0.002 | 241 |
| all_exons | go_visible | 2165 | 2.421 | 16.872 | 0.14 | 0.002 | 711 |
| cds | all | 2550 | 3.125 | 20.554 | 0.15 | 0.002 | 912 |
| cds | go_invisible | 502 | 4.806 | 26.914 | 0.18 | 0.002 | 230 |
| cds | go_visible | 2048 | 2.744 | 17.126 | 0.16 | 0.002 | 682 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2450 | 0.676 | 0.936 | 2.32e-152 |
| go_invisible | 476 | 0.758 | 0.936 | 1.94e-18 |
| go_visible | 1974 | 0.659 | 0.936 | 3.76e-141 |
