# Switch coding-consequence enrichment — gtex_hippocampus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **14894** (GO-invisible 577, GO-visible 14317)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.503 | 0.553 | 0.91 | 1 |
| all | cds_changed | 0.808 | 0.781 | 1.04 | 0.000999 |
| all | coding_consequence | 0.808 | 0.781 | 1.04 | 0.000999 |
| all | coding_status_change | 0.367 | 0.407 | 0.90 | 1 |
| all | first_exon_changed | 0.966 | 0.971 | 0.99 | 1 |
| all | internal_exon_difference | 0.915 | 0.931 | 0.98 | 1 |
| all | last_exon_changed | 0.921 | 0.959 | 0.96 | 1 |
| all | nmd_switch | 0.161 | 0.174 | 0.93 | 1 |
| all | utr_changed | 0.495 | 0.406 | 1.22 | 0.000999 |
| go_invisible | biotype_switch | 0.404 | 0.470 | 0.86 | 1 |
| go_invisible | cds_changed | 0.695 | 0.691 | 1.01 | 0.391 |
| go_invisible | coding_consequence | 0.695 | 0.691 | 1.01 | 0.391 |
| go_invisible | coding_status_change | 0.269 | 0.340 | 0.79 | 0.999 |
| go_invisible | first_exon_changed | 0.977 | 0.978 | 1.00 | 0.613 |
| go_invisible | internal_exon_difference | 0.882 | 0.924 | 0.95 | 1 |
| go_invisible | last_exon_changed | 0.958 | 0.982 | 0.98 | 1 |
| go_invisible | nmd_switch | 0.166 | 0.180 | 0.92 | 0.847 |
| go_invisible | utr_changed | 0.497 | 0.379 | 1.31 | 0.000999 |
| go_visible | biotype_switch | 0.507 | 0.556 | 0.91 | 1 |
| go_visible | cds_changed | 0.813 | 0.785 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.813 | 0.785 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.371 | 0.410 | 0.90 | 1 |
| go_visible | first_exon_changed | 0.966 | 0.971 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.916 | 0.931 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.920 | 0.959 | 0.96 | 1 |
| go_visible | nmd_switch | 0.161 | 0.174 | 0.93 | 1 |
| go_visible | utr_changed | 0.495 | 0.407 | 1.22 | 0.000999 |
