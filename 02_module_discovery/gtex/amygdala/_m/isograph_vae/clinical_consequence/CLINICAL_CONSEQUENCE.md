# Clinical-consequence of switched exons — gtex_amygdala

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **66876** (switched 63172, background 3704; CDS-overlapping 57296). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 3005 | 2.539 | 11.234 | 0.23 | 0.002 | 990 |
| all_exons | go_invisible | 694 | 1.778 | 7.720 | 0.23 | 0.002 | 176 |
| all_exons | go_visible | 2311 | 2.733 | 11.755 | 0.23 | 0.002 | 814 |
| cds | all | 2855 | 2.894 | 11.554 | 0.25 | 0.002 | 965 |
| cds | go_invisible | 610 | 2.180 | 7.898 | 0.28 | 0.002 | 169 |
| cds | go_visible | 2245 | 3.060 | 12.099 | 0.25 | 0.002 | 796 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2735 | 0.705 | 0.936 | 9.09e-119 |
| go_invisible | 591 | 0.928 | 0.936 | 0.0655 |
| go_visible | 2144 | 0.662 | 0.936 | 2.09e-144 |
