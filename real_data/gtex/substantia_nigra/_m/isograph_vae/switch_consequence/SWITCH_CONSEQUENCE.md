# Switch coding-consequence enrichment — gtex_substantia_nigra

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **1310** (GO-invisible 595, GO-visible 715)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.378 | 0.470 | 0.80 | 1 |
| all | cds_changed | 0.731 | 0.708 | 1.03 | 0.002 |
| all | coding_consequence | 0.731 | 0.708 | 1.03 | 0.002 |
| all | coding_status_change | 0.285 | 0.352 | 0.81 | 1 |
| all | first_exon_changed | 0.952 | 0.974 | 0.98 | 1 |
| all | internal_exon_difference | 0.931 | 0.951 | 0.98 | 1 |
| all | last_exon_changed | 0.887 | 0.956 | 0.93 | 1 |
| all | nmd_switch | 0.131 | 0.155 | 0.85 | 0.996 |
| all | utr_changed | 0.458 | 0.363 | 1.26 | 0.000999 |
| go_invisible | biotype_switch | 0.382 | 0.485 | 0.79 | 1 |
| go_invisible | cds_changed | 0.771 | 0.750 | 1.03 | 0.041 |
| go_invisible | coding_consequence | 0.771 | 0.750 | 1.03 | 0.041 |
| go_invisible | coding_status_change | 0.282 | 0.359 | 0.79 | 1 |
| go_invisible | first_exon_changed | 0.945 | 0.977 | 0.97 | 1 |
| go_invisible | internal_exon_difference | 0.906 | 0.949 | 0.95 | 1 |
| go_invisible | last_exon_changed | 0.899 | 0.965 | 0.93 | 1 |
| go_invisible | nmd_switch | 0.124 | 0.161 | 0.77 | 0.997 |
| go_invisible | utr_changed | 0.504 | 0.402 | 1.26 | 0.000999 |
| go_visible | biotype_switch | 0.375 | 0.459 | 0.82 | 1 |
| go_visible | cds_changed | 0.698 | 0.674 | 1.04 | 0.023 |
| go_visible | coding_consequence | 0.698 | 0.674 | 1.04 | 0.023 |
| go_visible | coding_status_change | 0.287 | 0.346 | 0.83 | 1 |
| go_visible | first_exon_changed | 0.958 | 0.971 | 0.99 | 0.979 |
| go_visible | internal_exon_difference | 0.951 | 0.952 | 1.00 | 0.594 |
| go_visible | last_exon_changed | 0.877 | 0.948 | 0.93 | 1 |
| go_visible | nmd_switch | 0.137 | 0.149 | 0.92 | 0.851 |
| go_visible | utr_changed | 0.420 | 0.331 | 1.27 | 0.000999 |
