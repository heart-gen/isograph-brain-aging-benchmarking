# Clinical-consequence of switched exons — gtex_cortex

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **65941** (switched 60869, background 5072; CDS-overlapping 55338). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 2900 | 2.753 | 18.098 | 0.15 | 0.002 | 1076 |
| all_exons | go_invisible | 873 | 3.411 | 22.242 | 0.15 | 0.002 | 360 |
| all_exons | go_visible | 2027 | 2.477 | 15.301 | 0.16 | 0.002 | 716 |
| cds | all | 2690 | 3.203 | 18.363 | 0.17 | 0.002 | 1027 |
| cds | go_invisible | 792 | 4.066 | 22.488 | 0.18 | 0.002 | 345 |
| cds | go_visible | 1898 | 2.852 | 15.563 | 0.18 | 0.002 | 682 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 2583 | 0.696 | 0.936 | 1.43e-134 |
| go_invisible | 758 | 0.762 | 0.936 | 2.65e-28 |
| go_visible | 1825 | 0.671 | 0.936 | 5.6e-114 |
