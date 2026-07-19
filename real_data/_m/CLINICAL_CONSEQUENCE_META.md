# Switched-exon clinical consequence — cross-region rollup

**Primary anchor = gnomAD LOEUF**: are switch genes more loss-of-function constrained than genome-wide (lower median LOEUF; Fisher-combined MWU p). The exon-level ClinVar columns are a direction-neutral secondary readout: how many regions have switched exons with higher (ratio>1) vs lower (ratio<1) P/LP density than constitutive exons at two-sided within-gene permutation p < 0.05, and the median ratio. `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only. Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the CDS scope is the fairer coding-vs-coding contrast.

| scope | stratum | regions | median LOEUF (switch) | LOEUF p | ratio>1 (p<.05) | ratio<1 (p<.05) | median ratio | Fisher p |
|-------|---------|---------|-----------------------|---------|-----------------|-----------------|--------------|----------|
| all_exons | all | 10 | 0.721 | 1.09e-93 | 0 | 10 | 0.18 | 2.09e-16 |
| all_exons | go_invisible | 10 | 0.700 | 1.82e-92 | 0 | 10 | 0.15 | 6.81e-16 |
| all_exons | go_visible | 9 | 0.748 | 1.05e-62 | 0 | 6 | 0.17 | 5.62e-10 |
| cds | all | 10 | 0.721 | 1.09e-93 | 0 | 10 | 0.21 | 2.67e-16 |
| cds | go_invisible | 10 | 0.700 | 1.82e-92 | 0 | 10 | 0.17 | 1.03e-15 |
| cds | go_visible | 9 | 0.748 | 1.05e-62 | 0 | 6 | 0.21 | 3.70e-09 |
