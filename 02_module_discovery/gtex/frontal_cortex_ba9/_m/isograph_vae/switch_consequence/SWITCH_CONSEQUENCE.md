# Switch coding-consequence enrichment — gtex_frontal_cortex_ba9

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **80065** (GO-invisible 36518, GO-visible 43547)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.472 | 0.526 | 0.90 | 1 |
| all | cds_changed | 0.766 | 0.736 | 1.04 | 0.000999 |
| all | coding_consequence | 0.766 | 0.736 | 1.04 | 0.000999 |
| all | coding_status_change | 0.335 | 0.384 | 0.87 | 1 |
| all | first_exon_changed | 0.966 | 0.971 | 0.99 | 1 |
| all | internal_exon_difference | 0.910 | 0.925 | 0.98 | 1 |
| all | last_exon_changed | 0.918 | 0.955 | 0.96 | 1 |
| all | nmd_switch | 0.161 | 0.169 | 0.95 | 1 |
| all | utr_changed | 0.483 | 0.382 | 1.26 | 0.000999 |
| go_invisible | biotype_switch | 0.453 | 0.505 | 0.90 | 1 |
| go_invisible | cds_changed | 0.723 | 0.691 | 1.05 | 0.000999 |
| go_invisible | coding_consequence | 0.723 | 0.691 | 1.05 | 0.000999 |
| go_invisible | coding_status_change | 0.315 | 0.362 | 0.87 | 1 |
| go_invisible | first_exon_changed | 0.964 | 0.971 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.907 | 0.926 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.914 | 0.955 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.164 | 0.173 | 0.95 | 1 |
| go_invisible | utr_changed | 0.457 | 0.354 | 1.29 | 0.000999 |
| go_visible | biotype_switch | 0.488 | 0.544 | 0.90 | 1 |
| go_visible | cds_changed | 0.802 | 0.773 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.802 | 0.773 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.352 | 0.402 | 0.87 | 1 |
| go_visible | first_exon_changed | 0.967 | 0.971 | 1.00 | 1 |
| go_visible | internal_exon_difference | 0.911 | 0.924 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.921 | 0.956 | 0.96 | 1 |
| go_visible | nmd_switch | 0.159 | 0.166 | 0.96 | 1 |
| go_visible | utr_changed | 0.505 | 0.405 | 1.25 | 0.000999 |
