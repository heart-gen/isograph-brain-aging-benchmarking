# Clinical-consequence of switched exons — gtex_substantia_nigra

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **2468** (switched 2216, background 252; CDS-overlapping 1926). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 90 | 2.368 | 25.856 | 0.09 | 0.002 | 35 |
| all_exons | go_invisible | 40 | 1.610 | 15.839 | 0.10 | 0.012 | 14 |
| all_exons | go_visible | 50 | 2.860 | 31.175 | 0.09 | 0.002 | 21 |
| cds | all | 74 | 2.904 | 26.609 | 0.11 | 0.002 | 33 |
| cds | go_invisible | 34 | 1.827 | 15.839 | 0.12 | 0.026 | 14 |
| cds | go_visible | 40 | 3.721 | 32.587 | 0.11 | 0.002 | 19 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 72 | 0.592 | 0.936 | 1.86e-10 |
| go_invisible | 34 | 0.609 | 0.936 | 2.55e-06 |
| go_visible | 38 | 0.533 | 0.936 | 8e-06 |
