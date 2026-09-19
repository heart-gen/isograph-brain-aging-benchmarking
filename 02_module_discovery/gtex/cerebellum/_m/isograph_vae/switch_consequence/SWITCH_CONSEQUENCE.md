# Switch coding-consequence enrichment — gtex_cerebellum

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **40031** (GO-invisible 14096, GO-visible 25935)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.523 | 0.562 | 0.93 | 1 |
| all | cds_changed | 0.785 | 0.756 | 1.04 | 0.000999 |
| all | coding_consequence | 0.785 | 0.756 | 1.04 | 0.000999 |
| all | coding_status_change | 0.389 | 0.411 | 0.95 | 1 |
| all | first_exon_changed | 0.971 | 0.976 | 0.99 | 1 |
| all | internal_exon_difference | 0.916 | 0.928 | 0.99 | 1 |
| all | last_exon_changed | 0.925 | 0.961 | 0.96 | 1 |
| all | nmd_switch | 0.168 | 0.181 | 0.93 | 1 |
| all | utr_changed | 0.436 | 0.369 | 1.18 | 0.000999 |
| go_invisible | biotype_switch | 0.527 | 0.556 | 0.95 | 1 |
| go_invisible | cds_changed | 0.746 | 0.732 | 1.02 | 0.000999 |
| go_invisible | coding_consequence | 0.746 | 0.732 | 1.02 | 0.000999 |
| go_invisible | coding_status_change | 0.397 | 0.405 | 0.98 | 0.982 |
| go_invisible | first_exon_changed | 0.972 | 0.976 | 1.00 | 0.997 |
| go_invisible | internal_exon_difference | 0.911 | 0.923 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.925 | 0.957 | 0.97 | 1 |
| go_invisible | nmd_switch | 0.166 | 0.184 | 0.91 | 1 |
| go_invisible | utr_changed | 0.391 | 0.352 | 1.11 | 0.000999 |
| go_visible | biotype_switch | 0.521 | 0.565 | 0.92 | 1 |
| go_visible | cds_changed | 0.806 | 0.769 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.806 | 0.769 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.385 | 0.414 | 0.93 | 1 |
| go_visible | first_exon_changed | 0.970 | 0.976 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.919 | 0.930 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.925 | 0.963 | 0.96 | 1 |
| go_visible | nmd_switch | 0.168 | 0.179 | 0.94 | 1 |
| go_visible | utr_changed | 0.460 | 0.378 | 1.21 | 0.000999 |
