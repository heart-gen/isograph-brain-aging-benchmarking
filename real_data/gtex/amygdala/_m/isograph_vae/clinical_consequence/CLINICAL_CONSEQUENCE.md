# Clinical-consequence of switched exons — gtex_amygdala

**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND (constitutive exons of the switching isoforms); `ratio` = switched / background; `p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by GO-invisible module membership.

- exons scored: **32807** (switched 30504, background 2303; CDS-overlapping 26496). `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only.

| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |
|-------|---------|-------|-------------|-------|-------|---|------------|
| all_exons | all | 1349 | 2.673 | 12.844 | 0.21 | 0.002 | 514 |
| all_exons | go_invisible | 568 | 3.620 | 14.596 | 0.25 | 0.002 | 221 |
| all_exons | go_visible | 781 | 1.904 | 11.056 | 0.17 | 0.002 | 293 |
| cds | all | 1166 | 3.264 | 13.236 | 0.25 | 0.002 | 485 |
| cds | go_invisible | 512 | 4.248 | 15.109 | 0.28 | 0.002 | 214 |
| cds | go_visible | 654 | 2.393 | 11.337 | 0.21 | 0.002 | 271 |

**gnomAD LOEUF** of the switch genes vs all genes in the constraint table (Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.

| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |
|---------|----------------|-----------------------|--------------|-------|
| all | 1122 | 0.605 | 0.936 | 4.35e-101 |
| go_invisible | 494 | 0.576 | 0.936 | 3.17e-53 |
| go_visible | 628 | 0.627 | 0.936 | 5.75e-52 |
