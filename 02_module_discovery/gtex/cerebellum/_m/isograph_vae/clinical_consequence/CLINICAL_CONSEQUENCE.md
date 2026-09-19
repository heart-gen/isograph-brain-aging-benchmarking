# Clinical-consequence of switched exons — gtex_cerebellum

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **69991** (switched 66465, background 3526; CDS-overlapping 59535). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 3120 | 2.458 | 11.534 | 0.21 | 0.002 | 920 |
| all_exons | go_invisible | 1083 | 2.199 | 13.589 | 0.16 | 0.002 | 288 |
| all_exons | go_visible | 2037 | 2.576 | 10.940 | 0.24 | 0.002 | 632 |
| cds | all | 2911 | 2.834 | 11.748 | 0.24 | 0.002 | 891 |
| cds | go_invisible | 985 | 2.604 | 13.894 | 0.19 | 0.002 | 273 |
| cds | go_visible | 1926 | 2.932 | 11.136 | 0.26 | 0.002 | 618 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2820 | 0.737 | 0.936 | 7.48e-95 |
| go_invisible | 948 | 0.855 | 0.936 | 1.76e-08 |
| go_visible | 1872 | 0.679 | 0.936 | 2.88e-105 |
