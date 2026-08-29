# Switch coding-consequence enrichment — brainseq_caudate_sczd

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **3440** (GO-invisible 2106, GO-visible 1334)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.412 | 0.459 | 0.90 | 1 |
| all | cds_changed | 0.669 | 0.642 | 1.04 | 0.000999 |
| all | coding_consequence | 0.669 | 0.642 | 1.04 | 0.000999 |
| all | coding_status_change | 0.300 | 0.336 | 0.89 | 1 |
| all | first_exon_changed | 0.956 | 0.967 | 0.99 | 1 |
| all | internal_exon_difference | 0.894 | 0.919 | 0.97 | 1 |
| all | last_exon_changed | 0.928 | 0.951 | 0.98 | 1 |
| all | nmd_switch | 0.131 | 0.144 | 0.91 | 0.995 |
| all | utr_changed | 0.430 | 0.340 | 1.26 | 0.000999 |
| go_invisible | biotype_switch | 0.434 | 0.481 | 0.90 | 1 |
| go_invisible | cds_changed | 0.670 | 0.649 | 1.03 | 0.000999 |
| go_invisible | coding_consequence | 0.670 | 0.649 | 1.03 | 0.000999 |
| go_invisible | coding_status_change | 0.318 | 0.353 | 0.90 | 1 |
| go_invisible | first_exon_changed | 0.958 | 0.968 | 0.99 | 0.994 |
| go_invisible | internal_exon_difference | 0.904 | 0.934 | 0.97 | 1 |
| go_invisible | last_exon_changed | 0.938 | 0.963 | 0.97 | 1 |
| go_invisible | nmd_switch | 0.136 | 0.153 | 0.89 | 0.997 |
| go_invisible | utr_changed | 0.415 | 0.325 | 1.28 | 0.000999 |
| go_visible | biotype_switch | 0.377 | 0.423 | 0.89 | 1 |
| go_visible | cds_changed | 0.667 | 0.630 | 1.06 | 0.000999 |
| go_visible | coding_consequence | 0.667 | 0.630 | 1.06 | 0.000999 |
| go_visible | coding_status_change | 0.271 | 0.310 | 0.88 | 1 |
| go_visible | first_exon_changed | 0.952 | 0.966 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.879 | 0.895 | 0.98 | 0.994 |
| go_visible | last_exon_changed | 0.911 | 0.933 | 0.98 | 1 |
| go_visible | nmd_switch | 0.122 | 0.129 | 0.95 | 0.794 |
| go_visible | utr_changed | 0.454 | 0.365 | 1.24 | 0.000999 |
