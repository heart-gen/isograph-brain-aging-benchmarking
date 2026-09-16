# Clinical-consequence of switched exons — brainseq_caudate_sczd

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **9601** (switched 8747, background 854; CDS-overlapping 7900). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 447 | 2.096 | 34.680 | 0.06 | 0.002 | 139 |
| all_exons | go_invisible | 414 | 1.875 | 34.872 | 0.05 | 0.002 | 129 |
| all_exons | go_visible | 33 | 7.692 | 18.450 | 0.42 | 0.617 | 10 |
| cds | all | 409 | 2.485 | 35.110 | 0.07 | 0.002 | 133 |
| cds | go_invisible | 376 | 2.227 | 35.309 | 0.06 | 0.002 | 123 |
| cds | go_visible | 33 | 8.430 | 18.450 | 0.46 | 0.41 | 10 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 397 | 0.832 | 0.936 | 8.59e-07 |
| go_invisible | 364 | 0.831 | 0.936 | 3.11e-07 |
| go_visible | 33 | 0.909 | 0.936 | 0.481 |
