# Switched-exon clinical consequence — cross-region rollup

**Primary anchor = gnomAD LOEUF**: are switch genes more loss-of-function constrained than genome-wide (lower median LOEUF; Fisher-combined MWU p). The exon-level ClinVar columns are a direction-neutral secondary readout: how many regions have switched exons with higher (ratio>1) vs lower (ratio<1) P/LP density than constitutive exons at two-sided within-gene permutation p < 0.05, and the median ratio. `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only. Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the CDS scope is the fairer coding-vs-coding contrast.

| scope | stratum | regions | median LOEUF (switch) | LOEUF p | ratio>1 (p<.05) | ratio<1 (p<.05) | median ratio | Fisher p |
|-------|---------|---------|-----------------------|---------|-----------------|-----------------|--------------|----------|
| all_exons | all | 10 | 0.729 | 1.28e-89 | 0 | 10 | 0.16 | 4.51e-17 |
| all_exons | go_invisible | 7 | 0.765 | 2.31e-54 | 0 | 6 | 0.16 | 8.17e-11 |
| all_exons | go_visible | 9 | 0.714 | 1.24e-64 | 0 | 8 | 0.16 | 5.67e-14 |
| cds | all | 10 | 0.729 | 1.28e-89 | 0 | 10 | 0.18 | 5.83e-16 |
| cds | go_invisible | 7 | 0.765 | 2.31e-54 | 0 | 7 | 0.18 | 1.29e-10 |
| cds | go_visible | 9 | 0.714 | 1.24e-64 | 0 | 8 | 0.19 | 4.04e-14 |
