# Clinical-consequence of switched exons — gtex_cerebellum

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **21002** (switched 19737, background 1265; CDS-overlapping 17177). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 853 | 3.014 | 16.930 | 0.18 | 0.002 | 276 |
| all_exons | go_invisible | 527 | 2.944 | 20.355 | 0.14 | 0.002 | 183 |
| all_exons | go_visible | 326 | 3.147 | 8.331 | 0.38 | 0.002 | 93 |
| cds | all | 752 | 3.647 | 17.199 | 0.21 | 0.002 | 264 |
| cds | go_invisible | 469 | 3.453 | 20.770 | 0.17 | 0.002 | 175 |
| cds | go_visible | 283 | 4.033 | 8.391 | 0.48 | 0.004 | 89 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 729 | 0.710 | 0.936 | 6.82e-44 |
| go_invisible | 455 | 0.635 | 0.936 | 7.31e-43 |
| go_visible | 274 | 0.804 | 0.936 | 2.1e-07 |
