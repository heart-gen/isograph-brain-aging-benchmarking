# Switch coding-consequence enrichment — gtex_hypothalamus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **7863** (GO-invisible 3256, GO-visible 4607)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.445 | 0.511 | 0.87 | 1 |
| all | cds_changed | 0.778 | 0.751 | 1.04 | 0.000999 |
| all | coding_consequence | 0.778 | 0.751 | 1.04 | 0.000999 |
| all | coding_status_change | 0.303 | 0.364 | 0.83 | 1 |
| all | first_exon_changed | 0.953 | 0.960 | 0.99 | 0.999 |
| all | internal_exon_difference | 0.899 | 0.924 | 0.97 | 1 |
| all | last_exon_changed | 0.906 | 0.949 | 0.95 | 1 |
| all | nmd_switch | 0.172 | 0.185 | 0.93 | 0.999 |
| all | utr_changed | 0.538 | 0.423 | 1.27 | 0.000999 |
| go_invisible | biotype_switch | 0.462 | 0.560 | 0.83 | 1 |
| go_invisible | cds_changed | 0.822 | 0.807 | 1.02 | 0.00799 |
| go_invisible | coding_consequence | 0.822 | 0.807 | 1.02 | 0.00799 |
| go_invisible | coding_status_change | 0.302 | 0.390 | 0.77 | 1 |
| go_invisible | first_exon_changed | 0.953 | 0.960 | 0.99 | 0.985 |
| go_invisible | internal_exon_difference | 0.905 | 0.937 | 0.97 | 1 |
| go_invisible | last_exon_changed | 0.885 | 0.947 | 0.94 | 1 |
| go_invisible | nmd_switch | 0.191 | 0.222 | 0.86 | 1 |
| go_invisible | utr_changed | 0.593 | 0.451 | 1.32 | 0.000999 |
| go_visible | biotype_switch | 0.433 | 0.476 | 0.91 | 1 |
| go_visible | cds_changed | 0.747 | 0.712 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.747 | 0.712 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.303 | 0.346 | 0.88 | 1 |
| go_visible | first_exon_changed | 0.954 | 0.959 | 0.99 | 0.991 |
| go_visible | internal_exon_difference | 0.895 | 0.915 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.920 | 0.950 | 0.97 | 1 |
| go_visible | nmd_switch | 0.159 | 0.158 | 1.00 | 0.46 |
| go_visible | utr_changed | 0.499 | 0.403 | 1.24 | 0.000999 |
