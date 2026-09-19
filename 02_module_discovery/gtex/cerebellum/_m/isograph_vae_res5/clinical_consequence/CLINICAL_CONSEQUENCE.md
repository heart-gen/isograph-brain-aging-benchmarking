# Clinical-consequence of switched exons — gtex_cerebellum

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **54204** (switched 51397, background 2807; CDS-overlapping 46592). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 2281 | 2.698 | 11.453 | 0.24 | 0.002 | 690 |
| all_exons | go_invisible | 668 | 2.440 | 14.669 | 0.17 | 0.002 | 177 |
| all_exons | go_visible | 1613 | 2.792 | 10.761 | 0.26 | 0.002 | 513 |
| cds | all | 2137 | 3.090 | 11.647 | 0.27 | 0.002 | 673 |
| cds | go_invisible | 601 | 2.890 | 14.830 | 0.19 | 0.002 | 170 |
| cds | go_visible | 1536 | 3.157 | 10.965 | 0.29 | 0.002 | 503 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2068 | 0.694 | 0.936 | 7.69e-102 |
| go_invisible | 575 | 0.843 | 0.936 | 1.37e-07 |
| go_visible | 1493 | 0.650 | 0.936 | 1.12e-109 |
