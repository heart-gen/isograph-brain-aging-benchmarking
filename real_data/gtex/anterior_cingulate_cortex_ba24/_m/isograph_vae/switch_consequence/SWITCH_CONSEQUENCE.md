# Switch coding-consequence enrichment — gtex_anterior_cingulate_cortex_ba24

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **43876** (GO-invisible 21081, GO-visible 22795)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.453 | 0.519 | 0.87 | 1 |
| all | cds_changed | 0.772 | 0.736 | 1.05 | 0.000999 |
| all | coding_consequence | 0.772 | 0.736 | 1.05 | 0.000999 |
| all | coding_status_change | 0.322 | 0.380 | 0.85 | 1 |
| all | first_exon_changed | 0.966 | 0.972 | 0.99 | 1 |
| all | internal_exon_difference | 0.913 | 0.930 | 0.98 | 1 |
| all | last_exon_changed | 0.924 | 0.960 | 0.96 | 1 |
| all | nmd_switch | 0.153 | 0.162 | 0.95 | 1 |
| all | utr_changed | 0.501 | 0.386 | 1.30 | 0.000999 |
| go_invisible | biotype_switch | 0.452 | 0.527 | 0.86 | 1 |
| go_invisible | cds_changed | 0.780 | 0.740 | 1.05 | 0.000999 |
| go_invisible | coding_consequence | 0.780 | 0.740 | 1.05 | 0.000999 |
| go_invisible | coding_status_change | 0.314 | 0.381 | 0.83 | 1 |
| go_invisible | first_exon_changed | 0.964 | 0.973 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.914 | 0.937 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.929 | 0.965 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.163 | 0.174 | 0.94 | 1 |
| go_invisible | utr_changed | 0.515 | 0.386 | 1.33 | 0.000999 |
| go_visible | biotype_switch | 0.453 | 0.511 | 0.89 | 1 |
| go_visible | cds_changed | 0.765 | 0.733 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.765 | 0.733 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.329 | 0.379 | 0.87 | 1 |
| go_visible | first_exon_changed | 0.968 | 0.972 | 1.00 | 1 |
| go_visible | internal_exon_difference | 0.911 | 0.923 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.919 | 0.956 | 0.96 | 1 |
| go_visible | nmd_switch | 0.144 | 0.151 | 0.96 | 1 |
| go_visible | utr_changed | 0.489 | 0.386 | 1.27 | 0.000999 |
