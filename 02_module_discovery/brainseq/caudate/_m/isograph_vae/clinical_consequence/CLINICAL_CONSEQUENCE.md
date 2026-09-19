# Clinical-consequence of switched exons — brainseq_caudate

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **15733** (switched 14665, background 1068; CDS-overlapping 12628). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 743 | 1.837 | 11.814 | 0.16 | 0.002 | 215 |
| all_exons | go_visible | 743 | 1.837 | 11.814 | 0.16 | 0.002 | 215 |
| cds | all | 642 | 2.224 | 11.970 | 0.19 | 0.002 | 204 |
| cds | go_visible | 642 | 2.224 | 11.970 | 0.19 | 0.002 | 204 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 623 | 0.784 | 0.936 | 2.82e-17 |
| go_invisible | 0 | nan | 0.936 | nan |
| go_visible | 623 | 0.784 | 0.936 | 2.82e-17 |
