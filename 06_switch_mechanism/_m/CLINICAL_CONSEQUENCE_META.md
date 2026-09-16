# Switched-exon clinical consequence — cross-region rollup

**Primary anchor = gnomAD LOEUF**: are switch genes more loss-of-function constrained than genome-wide (lower median LOEUF; Fisher-combined MWU p). The exon-level ClinVar columns are a direction-neutral secondary readout: how many regions have switched exons with higher (ratio>1) vs lower (ratio<1) P/LP density than constitutive exons at two-sided within-gene permutation p < 0.05, and the median ratio. `scope` = all switch-pair exons vs coding (CDS-overlapping) exons only. Alt-spliced exons are usually less constrained, so ratio < 1 is the expected baseline; the CDS scope is the fairer coding-vs-coding contrast.

| scope | stratum | regions | median LOEUF (switch) | LOEUF p | ratio>1 (p<.05) | ratio<1 (p<.05) | median ratio | Fisher p |
|-------|---------|---------|-----------------------|---------|-----------------|-----------------|--------------|----------|
| all_exons | all | 10 | 0.701 | 2.27e-96 | 0 | 10 | 0.16 | 1.48e-16 |
| all_exons | go_invisible | 10 | 0.763 | 9.44e-60 | 0 | 10 | 0.16 | 6.81e-16 |
| all_exons | go_visible | 8 | 0.676 | 4.55e-66 | 0 | 7 | 0.20 | 5.88e-12 |
| cds | all | 10 | 0.701 | 2.27e-96 | 0 | 10 | 0.18 | 3.77e-16 |
| cds | go_invisible | 10 | 0.763 | 9.44e-60 | 0 | 10 | 0.18 | 3.79e-15 |
| cds | go_visible | 8 | 0.676 | 4.55e-66 | 0 | 7 | 0.22 | 4.15e-12 |
