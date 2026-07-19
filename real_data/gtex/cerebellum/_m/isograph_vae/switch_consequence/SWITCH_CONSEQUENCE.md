# Switch coding-consequence enrichment — gtex_cerebellum

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **11829** (GO-invisible 7621, GO-visible 4208)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.475 | 0.512 | 0.93 | 1 |
| all | cds_changed | 0.742 | 0.714 | 1.04 | 0.000999 |
| all | coding_consequence | 0.742 | 0.714 | 1.04 | 0.000999 |
| all | coding_status_change | 0.346 | 0.375 | 0.92 | 1 |
| all | first_exon_changed | 0.965 | 0.967 | 1.00 | 0.917 |
| all | internal_exon_difference | 0.924 | 0.944 | 0.98 | 1 |
| all | last_exon_changed | 0.918 | 0.957 | 0.96 | 1 |
| all | nmd_switch | 0.160 | 0.172 | 0.93 | 1 |
| all | utr_changed | 0.442 | 0.362 | 1.22 | 0.000999 |
| go_invisible | biotype_switch | 0.502 | 0.539 | 0.93 | 1 |
| go_invisible | cds_changed | 0.752 | 0.725 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.752 | 0.725 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.365 | 0.394 | 0.93 | 1 |
| go_invisible | first_exon_changed | 0.967 | 0.971 | 1.00 | 0.979 |
| go_invisible | internal_exon_difference | 0.934 | 0.953 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.914 | 0.957 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.175 | 0.191 | 0.92 | 1 |
| go_invisible | utr_changed | 0.430 | 0.350 | 1.23 | 0.000999 |
| go_visible | biotype_switch | 0.427 | 0.464 | 0.92 | 1 |
| go_visible | cds_changed | 0.725 | 0.696 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.725 | 0.696 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.311 | 0.341 | 0.91 | 1 |
| go_visible | first_exon_changed | 0.962 | 0.961 | 1.00 | 0.382 |
| go_visible | internal_exon_difference | 0.907 | 0.930 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.927 | 0.956 | 0.97 | 1 |
| go_visible | nmd_switch | 0.131 | 0.137 | 0.95 | 0.915 |
| go_visible | utr_changed | 0.464 | 0.385 | 1.21 | 0.000999 |
