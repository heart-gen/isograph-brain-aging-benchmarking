# Clinical-consequence of switched exons — gtex_hypothalamus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **28038** (switched 25890, background 2148; CDS-overlapping 22932). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1243 | 3.195 | 15.310 | 0.21 | 0.002 | 497 |
| all_exons | go_invisible | 1156 | 3.321 | 15.627 | 0.21 | 0.00799 | 472 |
| all_exons | go_visible | 87 | 1.124 | 2.204 | 0.51 | 0.0759 | 25 |
| cds | all | 1134 | 3.783 | 15.860 | 0.24 | 0.02 | 472 |
| cds | go_invisible | 1057 | 3.921 | 16.160 | 0.24 | 0.018 | 452 |
| cds | go_visible | 77 | 1.399 | 2.465 | 0.57 | 0.102 | 20 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 1093 | 0.758 | 0.936 | 8.64e-43 |
| go_invisible | 1019 | 0.738 | 0.936 | 1.43e-49 |
| go_visible | 74 | 1.036 | 0.936 | 0.987 |
