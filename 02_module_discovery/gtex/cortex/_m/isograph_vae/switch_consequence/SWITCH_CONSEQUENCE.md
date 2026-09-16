# Switch coding-consequence enrichment — gtex_cortex

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **34061** (GO-invisible 7330, GO-visible 26731)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.498 | 0.555 | 0.90 | 1 |
| all | cds_changed | 0.801 | 0.768 | 1.04 | 0.000999 |
| all | coding_consequence | 0.801 | 0.768 | 1.04 | 0.000999 |
| all | coding_status_change | 0.359 | 0.405 | 0.89 | 1 |
| all | first_exon_changed | 0.968 | 0.974 | 0.99 | 1 |
| all | internal_exon_difference | 0.923 | 0.935 | 0.99 | 1 |
| all | last_exon_changed | 0.916 | 0.958 | 0.96 | 1 |
| all | nmd_switch | 0.165 | 0.177 | 0.93 | 1 |
| all | utr_changed | 0.489 | 0.391 | 1.25 | 0.000999 |
| go_invisible | biotype_switch | 0.457 | 0.528 | 0.87 | 1 |
| go_invisible | cds_changed | 0.773 | 0.749 | 1.03 | 0.000999 |
| go_invisible | coding_consequence | 0.773 | 0.749 | 1.03 | 0.000999 |
| go_invisible | coding_status_change | 0.319 | 0.377 | 0.85 | 1 |
| go_invisible | first_exon_changed | 0.959 | 0.971 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.916 | 0.935 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.884 | 0.947 | 0.93 | 1 |
| go_invisible | nmd_switch | 0.164 | 0.184 | 0.89 | 1 |
| go_invisible | utr_changed | 0.514 | 0.410 | 1.25 | 0.000999 |
| go_visible | biotype_switch | 0.510 | 0.562 | 0.91 | 1 |
| go_visible | cds_changed | 0.808 | 0.773 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.808 | 0.773 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.371 | 0.413 | 0.90 | 1 |
| go_visible | first_exon_changed | 0.971 | 0.974 | 1.00 | 1 |
| go_visible | internal_exon_difference | 0.924 | 0.935 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.925 | 0.962 | 0.96 | 1 |
| go_visible | nmd_switch | 0.165 | 0.175 | 0.94 | 1 |
| go_visible | utr_changed | 0.482 | 0.386 | 1.25 | 0.000999 |
