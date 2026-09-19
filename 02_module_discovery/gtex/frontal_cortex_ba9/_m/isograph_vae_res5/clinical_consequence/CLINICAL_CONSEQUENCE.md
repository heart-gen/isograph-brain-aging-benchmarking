# Clinical-consequence of switched exons — gtex_frontal_cortex_ba9

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **102536** (switched 96099, background 6437; CDS-overlapping 86331). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 4626 | 2.410 | 15.074 | 0.16 | 0.002 | 1569 |
| all_exons | go_invisible | 1974 | 2.501 | 16.978 | 0.15 | 0.002 | 631 |
| all_exons | go_visible | 2652 | 2.334 | 13.424 | 0.17 | 0.002 | 938 |
| cds | all | 4304 | 2.777 | 15.470 | 0.18 | 0.002 | 1503 |
| cds | go_invisible | 1810 | 2.906 | 17.475 | 0.17 | 0.002 | 608 |
| cds | go_visible | 2494 | 2.670 | 13.742 | 0.19 | 0.002 | 895 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 4136 | 0.703 | 0.936 | 1.55e-180 |
| go_invisible | 1742 | 0.711 | 0.936 | 3.86e-83 |
| go_visible | 2394 | 0.699 | 0.936 | 3.96e-109 |
