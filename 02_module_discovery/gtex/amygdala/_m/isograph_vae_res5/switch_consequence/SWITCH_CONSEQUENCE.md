# Switch coding-consequence enrichment — gtex_amygdala

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **37560** (GO-invisible 8851, GO-visible 28709)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.514 | 0.562 | 0.91 | 1 |
| all | cds_changed | 0.803 | 0.768 | 1.05 | 0.000999 |
| all | coding_consequence | 0.803 | 0.768 | 1.05 | 0.000999 |
| all | coding_status_change | 0.372 | 0.412 | 0.90 | 1 |
| all | first_exon_changed | 0.971 | 0.974 | 1.00 | 1 |
| all | internal_exon_difference | 0.920 | 0.931 | 0.99 | 1 |
| all | last_exon_changed | 0.919 | 0.960 | 0.96 | 1 |
| all | nmd_switch | 0.171 | 0.178 | 0.96 | 1 |
| all | utr_changed | 0.479 | 0.386 | 1.24 | 0.000999 |
| go_invisible | biotype_switch | 0.506 | 0.532 | 0.95 | 1 |
| go_invisible | cds_changed | 0.717 | 0.689 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.717 | 0.689 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.374 | 0.390 | 0.96 | 0.999 |
| go_invisible | first_exon_changed | 0.972 | 0.975 | 1.00 | 0.972 |
| go_invisible | internal_exon_difference | 0.906 | 0.917 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.918 | 0.957 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.162 | 0.170 | 0.95 | 0.997 |
| go_invisible | utr_changed | 0.385 | 0.325 | 1.19 | 0.000999 |
| go_visible | biotype_switch | 0.516 | 0.572 | 0.90 | 1 |
| go_visible | cds_changed | 0.829 | 0.793 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.829 | 0.793 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.371 | 0.418 | 0.89 | 1 |
| go_visible | first_exon_changed | 0.970 | 0.974 | 1.00 | 1 |
| go_visible | internal_exon_difference | 0.925 | 0.935 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.920 | 0.961 | 0.96 | 1 |
| go_visible | nmd_switch | 0.174 | 0.180 | 0.97 | 0.997 |
| go_visible | utr_changed | 0.508 | 0.405 | 1.26 | 0.000999 |
