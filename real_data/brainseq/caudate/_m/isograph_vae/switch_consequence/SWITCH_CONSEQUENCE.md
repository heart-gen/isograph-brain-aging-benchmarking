# Switch coding-consequence enrichment — brainseq_caudate

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **12703** (GO-invisible 12703, GO-visible 0)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.475 | 0.529 | 0.90 | 1 |
| all | cds_changed | 0.756 | 0.719 | 1.05 | 0.000999 |
| all | coding_consequence | 0.756 | 0.719 | 1.05 | 0.000999 |
| all | coding_status_change | 0.347 | 0.390 | 0.89 | 1 |
| all | first_exon_changed | 0.968 | 0.974 | 0.99 | 1 |
| all | internal_exon_difference | 0.912 | 0.924 | 0.99 | 1 |
| all | last_exon_changed | 0.948 | 0.962 | 0.98 | 1 |
| all | nmd_switch | 0.148 | 0.160 | 0.92 | 1 |
| all | utr_changed | 0.452 | 0.356 | 1.27 | 0.000999 |
| go_invisible | biotype_switch | 0.475 | 0.529 | 0.90 | 1 |
| go_invisible | cds_changed | 0.756 | 0.719 | 1.05 | 0.000999 |
| go_invisible | coding_consequence | 0.756 | 0.719 | 1.05 | 0.000999 |
| go_invisible | coding_status_change | 0.347 | 0.390 | 0.89 | 1 |
| go_invisible | first_exon_changed | 0.968 | 0.974 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.912 | 0.924 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.948 | 0.963 | 0.98 | 1 |
| go_invisible | nmd_switch | 0.148 | 0.160 | 0.92 | 1 |
| go_invisible | utr_changed | 0.452 | 0.356 | 1.27 | 0.000999 |
