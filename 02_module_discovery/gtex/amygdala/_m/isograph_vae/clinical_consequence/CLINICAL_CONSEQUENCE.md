# Clinical-consequence of switched exons — gtex_amygdala

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **1358** (switched 1295, background 63; CDS-overlapping 1164). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 68 | 4.198 | 29.616 | 0.14 | 0.002 | 23 |
| all_exons | go_visible | 68 | 4.198 | 29.616 | 0.14 | 0.004 | 23 |
| cds | all | 66 | 4.850 | 32.180 | 0.15 | 0.004 | 22 |
| cds | go_visible | 66 | 4.850 | 32.180 | 0.15 | 0.002 | 22 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 61 | 0.791 | 0.936 | 0.047 |
| go_invisible | 0 | nan | 0.936 | nan |
| go_visible | 61 | 0.791 | 0.936 | 0.047 |
