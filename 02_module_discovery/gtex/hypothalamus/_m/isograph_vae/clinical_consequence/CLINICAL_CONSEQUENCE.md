# Clinical-consequence of switched exons — gtex_hypothalamus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **12048** (switched 11192, background 856; CDS-overlapping 9748). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 580 | 2.906 | 28.424 | 0.10 | 0.002 | 230 |
| all_exons | go_invisible | 234 | 3.634 | 35.795 | 0.10 | 0.002 | 109 |
| all_exons | go_visible | 346 | 2.302 | 16.277 | 0.14 | 0.002 | 121 |
| cds | all | 527 | 3.412 | 28.751 | 0.12 | 0.002 | 224 |
| cds | go_invisible | 225 | 4.096 | 36.141 | 0.11 | 0.002 | 107 |
| cds | go_visible | 302 | 2.795 | 16.514 | 0.17 | 0.004 | 117 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 509 | 0.803 | 0.936 | 1.47e-13 |
| go_invisible | 219 | 0.738 | 0.936 | 1.52e-13 |
| go_visible | 290 | 0.857 | 0.936 | 0.000383 |
