# Clinical-consequence of switched exons — brainseq_caudate_sczd

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **5364** (switched 4904, background 460; CDS-overlapping 4080). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 251 | 3.175 | 16.893 | 0.19 | 0.004 | 80 |
| all_exons | go_invisible | 150 | 4.014 | 18.366 | 0.22 | 0.002 | 53 |
| all_exons | go_visible | 101 | 1.342 | 11.536 | 0.12 | 0.949 | 27 |
| cds | all | 203 | 4.146 | 17.178 | 0.24 | 0.002 | 73 |
| cds | go_invisible | 122 | 5.188 | 18.695 | 0.28 | 0.002 | 47 |
| cds | go_visible | 81 | 1.816 | 11.689 | 0.16 | 0.921 | 26 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 194 | 0.845 | 0.936 | 0.000517 |
| go_invisible | 117 | 0.756 | 0.936 | 7.02e-07 |
| go_visible | 77 | 0.965 | 0.936 | 0.768 |
