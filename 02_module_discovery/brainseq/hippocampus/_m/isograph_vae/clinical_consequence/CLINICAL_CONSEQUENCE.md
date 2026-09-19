# Clinical-consequence of switched exons — brainseq_hippocampus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **1847** (switched 1607, background 240; CDS-overlapping 1640). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 75 | 6.619 | 40.157 | 0.16 | 0.002 | 35 |
| all_exons | go_invisible | 75 | 6.619 | 40.157 | 0.16 | 0.002 | 35 |
| cds | all | 74 | 7.172 | 40.157 | 0.18 | 0.002 | 35 |
| cds | go_invisible | 74 | 7.172 | 40.157 | 0.18 | 0.002 | 35 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 71 | 0.627 | 0.936 | 3.43e-08 |
| go_invisible | 71 | 0.627 | 0.936 | 3.43e-08 |
| go_visible | 0 | nan | 0.936 | nan |
