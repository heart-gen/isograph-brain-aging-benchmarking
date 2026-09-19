# Switch coding-consequence enrichment — brainseq_caudate_sczd

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **22483** (GO-invisible 18263, GO-visible 4220)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.498 | 0.530 | 0.94 | 1 |
| all | cds_changed | 0.715 | 0.687 | 1.04 | 0.000999 |
| all | coding_consequence | 0.715 | 0.687 | 1.04 | 0.000999 |
| all | coding_status_change | 0.368 | 0.384 | 0.96 | 1 |
| all | first_exon_changed | 0.969 | 0.974 | 0.99 | 1 |
| all | internal_exon_difference | 0.917 | 0.928 | 0.99 | 1 |
| all | last_exon_changed | 0.941 | 0.958 | 0.98 | 1 |
| all | nmd_switch | 0.157 | 0.176 | 0.90 | 1 |
| all | utr_changed | 0.390 | 0.329 | 1.19 | 0.000999 |
| go_invisible | biotype_switch | 0.494 | 0.527 | 0.94 | 1 |
| go_invisible | cds_changed | 0.713 | 0.682 | 1.05 | 0.000999 |
| go_invisible | coding_consequence | 0.713 | 0.682 | 1.05 | 0.000999 |
| go_invisible | coding_status_change | 0.360 | 0.380 | 0.95 | 1 |
| go_invisible | first_exon_changed | 0.969 | 0.974 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.921 | 0.931 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.943 | 0.961 | 0.98 | 1 |
| go_invisible | nmd_switch | 0.162 | 0.177 | 0.91 | 1 |
| go_invisible | utr_changed | 0.394 | 0.326 | 1.21 | 0.000999 |
| go_visible | biotype_switch | 0.517 | 0.542 | 0.95 | 1 |
| go_visible | cds_changed | 0.722 | 0.711 | 1.02 | 0.024 |
| go_visible | coding_consequence | 0.722 | 0.711 | 1.02 | 0.024 |
| go_visible | coding_status_change | 0.402 | 0.405 | 0.99 | 0.665 |
| go_visible | first_exon_changed | 0.967 | 0.970 | 1.00 | 0.821 |
| go_visible | internal_exon_difference | 0.902 | 0.918 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.935 | 0.949 | 0.99 | 1 |
| go_visible | nmd_switch | 0.137 | 0.169 | 0.81 | 1 |
| go_visible | utr_changed | 0.374 | 0.344 | 1.09 | 0.000999 |
