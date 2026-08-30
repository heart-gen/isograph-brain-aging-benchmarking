# Switch coding-consequence enrichment — gtex_cortex

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **20253** (GO-invisible 11891, GO-visible 8362)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.465 | 0.533 | 0.87 | 1 |
| all | cds_changed | 0.784 | 0.748 | 1.05 | 0.000999 |
| all | coding_consequence | 0.784 | 0.748 | 1.05 | 0.000999 |
| all | coding_status_change | 0.323 | 0.382 | 0.85 | 1 |
| all | first_exon_changed | 0.965 | 0.971 | 0.99 | 1 |
| all | internal_exon_difference | 0.919 | 0.938 | 0.98 | 1 |
| all | last_exon_changed | 0.919 | 0.959 | 0.96 | 1 |
| all | nmd_switch | 0.171 | 0.183 | 0.93 | 1 |
| all | utr_changed | 0.515 | 0.395 | 1.30 | 0.000999 |
| go_invisible | biotype_switch | 0.497 | 0.573 | 0.87 | 1 |
| go_invisible | cds_changed | 0.832 | 0.788 | 1.06 | 0.000999 |
| go_invisible | coding_consequence | 0.832 | 0.788 | 1.06 | 0.000999 |
| go_invisible | coding_status_change | 0.347 | 0.411 | 0.84 | 1 |
| go_invisible | first_exon_changed | 0.963 | 0.972 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.928 | 0.945 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.926 | 0.965 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.185 | 0.200 | 0.93 | 1 |
| go_invisible | utr_changed | 0.539 | 0.404 | 1.33 | 0.000999 |
| go_visible | biotype_switch | 0.421 | 0.476 | 0.88 | 1 |
| go_visible | cds_changed | 0.715 | 0.691 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.715 | 0.691 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.290 | 0.340 | 0.85 | 1 |
| go_visible | first_exon_changed | 0.967 | 0.970 | 1.00 | 0.951 |
| go_visible | internal_exon_difference | 0.906 | 0.929 | 0.97 | 1 |
| go_visible | last_exon_changed | 0.908 | 0.949 | 0.96 | 1 |
| go_visible | nmd_switch | 0.151 | 0.160 | 0.94 | 0.997 |
| go_visible | utr_changed | 0.481 | 0.383 | 1.26 | 0.000999 |
