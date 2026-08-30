# Clinical-consequence of switched exons — gtex_cortex

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **34222** (switched 31783, background 2439; CDS-overlapping 28416). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1467 | 3.150 | 19.406 | 0.16 | 0.002 | 576 |
| all_exons | go_invisible | 859 | 3.145 | 13.639 | 0.23 | 0.002 | 360 |
| all_exons | go_visible | 608 | 3.156 | 28.239 | 0.11 | 0.002 | 216 |
| cds | all | 1344 | 3.674 | 19.558 | 0.19 | 0.002 | 560 |
| cds | go_invisible | 829 | 3.478 | 13.678 | 0.25 | 0.002 | 355 |
| cds | go_visible | 515 | 4.002 | 28.667 | 0.14 | 0.002 | 205 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 1303 | 0.683 | 0.936 | 5.69e-84 |
| go_invisible | 807 | 0.648 | 0.936 | 4.52e-68 |
| go_visible | 496 | 0.754 | 0.936 | 1.34e-21 |
