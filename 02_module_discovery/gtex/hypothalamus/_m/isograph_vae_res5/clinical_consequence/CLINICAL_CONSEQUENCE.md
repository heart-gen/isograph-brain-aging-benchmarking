# Clinical-consequence of switched exons — gtex_hypothalamus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **8768** (switched 7905, background 863; CDS-overlapping 7249). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 373 | 3.155 | 16.795 | 0.19 | 0.00799 | 176 |
| all_exons | go_invisible | 373 | 3.155 | 16.795 | 0.19 | 0.00799 | 176 |
| cds | all | 352 | 3.646 | 16.934 | 0.22 | 0.024 | 173 |
| cds | go_invisible | 352 | 3.646 | 16.934 | 0.22 | 0.028 | 173 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 339 | 0.699 | 0.936 | 6.82e-23 |
| go_invisible | 339 | 0.699 | 0.936 | 6.82e-23 |
| go_visible | 0 | nan | 0.936 | nan |
