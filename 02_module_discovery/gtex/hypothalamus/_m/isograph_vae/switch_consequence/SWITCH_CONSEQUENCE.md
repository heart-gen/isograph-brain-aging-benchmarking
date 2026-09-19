# Switch coding-consequence enrichment — gtex_hypothalamus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **16830** (GO-invisible 15793, GO-visible 1037)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.457 | 0.527 | 0.87 | 1 |
| all | cds_changed | 0.776 | 0.749 | 1.04 | 0.000999 |
| all | coding_consequence | 0.776 | 0.749 | 1.04 | 0.000999 |
| all | coding_status_change | 0.308 | 0.374 | 0.83 | 1 |
| all | first_exon_changed | 0.964 | 0.968 | 1.00 | 1 |
| all | internal_exon_difference | 0.909 | 0.930 | 0.98 | 1 |
| all | last_exon_changed | 0.889 | 0.943 | 0.94 | 1 |
| all | nmd_switch | 0.177 | 0.189 | 0.93 | 1 |
| all | utr_changed | 0.531 | 0.411 | 1.29 | 0.000999 |
| go_invisible | biotype_switch | 0.454 | 0.527 | 0.86 | 1 |
| go_invisible | cds_changed | 0.779 | 0.752 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.779 | 0.752 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.306 | 0.373 | 0.82 | 1 |
| go_invisible | first_exon_changed | 0.963 | 0.968 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.910 | 0.932 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.887 | 0.942 | 0.94 | 1 |
| go_invisible | nmd_switch | 0.179 | 0.192 | 0.93 | 1 |
| go_invisible | utr_changed | 0.538 | 0.414 | 1.30 | 0.000999 |
| go_visible | biotype_switch | 0.503 | 0.527 | 0.95 | 0.978 |
| go_visible | cds_changed | 0.728 | 0.714 | 1.02 | 0.0809 |
| go_visible | coding_consequence | 0.728 | 0.714 | 1.02 | 0.0809 |
| go_visible | coding_status_change | 0.345 | 0.381 | 0.91 | 0.993 |
| go_visible | first_exon_changed | 0.968 | 0.971 | 1.00 | 0.757 |
| go_visible | internal_exon_difference | 0.898 | 0.903 | 0.99 | 0.831 |
| go_visible | last_exon_changed | 0.916 | 0.947 | 0.97 | 1 |
| go_visible | nmd_switch | 0.151 | 0.153 | 0.99 | 0.594 |
| go_visible | utr_changed | 0.427 | 0.357 | 1.20 | 0.000999 |
