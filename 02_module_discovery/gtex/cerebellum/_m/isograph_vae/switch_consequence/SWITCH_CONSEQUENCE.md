# Switch coding-consequence enrichment — gtex_cerebellum

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **29561** (GO-invisible 8864, GO-visible 20697)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.526 | 0.566 | 0.93 | 1 |
| all | cds_changed | 0.792 | 0.759 | 1.04 | 0.000999 |
| all | coding_consequence | 0.792 | 0.759 | 1.04 | 0.000999 |
| all | coding_status_change | 0.392 | 0.415 | 0.94 | 1 |
| all | first_exon_changed | 0.970 | 0.976 | 0.99 | 1 |
| all | internal_exon_difference | 0.922 | 0.932 | 0.99 | 1 |
| all | last_exon_changed | 0.923 | 0.961 | 0.96 | 1 |
| all | nmd_switch | 0.170 | 0.182 | 0.93 | 1 |
| all | utr_changed | 0.437 | 0.365 | 1.20 | 0.000999 |
| go_invisible | biotype_switch | 0.523 | 0.550 | 0.95 | 1 |
| go_invisible | cds_changed | 0.744 | 0.728 | 1.02 | 0.000999 |
| go_invisible | coding_consequence | 0.744 | 0.728 | 1.02 | 0.000999 |
| go_invisible | coding_status_change | 0.390 | 0.400 | 0.97 | 0.986 |
| go_invisible | first_exon_changed | 0.969 | 0.975 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.916 | 0.926 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.928 | 0.958 | 0.97 | 1 |
| go_invisible | nmd_switch | 0.173 | 0.188 | 0.92 | 1 |
| go_invisible | utr_changed | 0.393 | 0.349 | 1.13 | 0.000999 |
| go_visible | biotype_switch | 0.528 | 0.573 | 0.92 | 1 |
| go_visible | cds_changed | 0.812 | 0.772 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.812 | 0.772 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.393 | 0.421 | 0.93 | 1 |
| go_visible | first_exon_changed | 0.970 | 0.976 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.925 | 0.935 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.920 | 0.962 | 0.96 | 1 |
| go_visible | nmd_switch | 0.168 | 0.180 | 0.94 | 1 |
| go_visible | utr_changed | 0.455 | 0.372 | 1.22 | 0.000999 |
