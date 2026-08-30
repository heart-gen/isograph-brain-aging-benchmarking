# Clinical-consequence of switched exons — brainseq_caudate

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **22481** (switched 20785, background 1696; CDS-overlapping 18087). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1067 | 1.520 | 9.550 | 0.16 | 0.002 | 345 |
| all_exons | go_invisible | 1067 | 1.520 | 9.550 | 0.16 | 0.002 | 345 |
| cds | all | 932 | 1.828 | 10.031 | 0.18 | 0.002 | 316 |
| cds | go_invisible | 932 | 1.828 | 10.031 | 0.18 | 0.002 | 316 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 886 | 0.768 | 0.936 | 3.32e-28 |
| go_invisible | 886 | 0.768 | 0.936 | 3.32e-28 |
| go_visible | 0 | nan | 0.936 | nan |
