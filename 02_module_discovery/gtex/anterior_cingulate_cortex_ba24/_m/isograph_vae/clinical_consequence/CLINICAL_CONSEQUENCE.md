# Clinical-consequence of switched exons — gtex_anterior_cingulate_cortex_ba24

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **60904** (switched 56483, background 4421; CDS-overlapping 51252). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 2856 | 2.383 | 20.531 | 0.12 | 0.002 | 1041 |
| all_exons | go_visible | 2856 | 2.383 | 20.531 | 0.12 | 0.002 | 1041 |
| cds | all | 2675 | 2.740 | 20.831 | 0.13 | 0.002 | 997 |
| cds | go_visible | 2675 | 2.740 | 20.831 | 0.13 | 0.002 | 997 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2570 | 0.703 | 0.936 | 1.94e-122 |
| go_invisible | 0 | nan | 0.936 | nan |
| go_visible | 2570 | 0.703 | 0.936 | 1.94e-122 |
