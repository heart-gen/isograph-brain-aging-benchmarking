# Clinical-consequence of switched exons — gtex_hippocampus

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **26518** (switched 24753, background 1765; CDS-overlapping 22512). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1274 | 1.962 | 8.483 | 0.23 | 0.002 | 474 |
| all_exons | go_invisible | 43 | 1.912 | 14.803 | 0.13 | 0.0579 | 22 |
| all_exons | go_visible | 1231 | 1.964 | 8.176 | 0.24 | 0.002 | 452 |
| cds | all | 1219 | 2.224 | 8.607 | 0.26 | 0.002 | 453 |
| cds | go_invisible | 36 | 2.447 | 15.143 | 0.16 | 0.044 | 20 |
| cds | go_visible | 1183 | 2.217 | 8.292 | 0.27 | 0.002 | 433 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 1162 | 0.673 | 0.936 | 4.2e-73 |
| go_invisible | 33 | 0.924 | 0.936 | 0.264 |
| go_visible | 1129 | 0.665 | 0.936 | 2.14e-74 |
