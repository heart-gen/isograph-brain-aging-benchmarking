# Switch coding-consequence enrichment — brainseq_caudate_sczd

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **5440** (GO-invisible 5089, GO-visible 351)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.477 | 0.529 | 0.90 | 1 |
| all | cds_changed | 0.762 | 0.737 | 1.03 | 0.000999 |
| all | coding_consequence | 0.762 | 0.737 | 1.03 | 0.000999 |
| all | coding_status_change | 0.351 | 0.396 | 0.89 | 1 |
| all | first_exon_changed | 0.970 | 0.976 | 0.99 | 0.997 |
| all | internal_exon_difference | 0.926 | 0.931 | 0.99 | 0.956 |
| all | last_exon_changed | 0.929 | 0.954 | 0.97 | 1 |
| all | nmd_switch | 0.153 | 0.160 | 0.96 | 0.937 |
| all | utr_changed | 0.462 | 0.371 | 1.24 | 0.000999 |
| go_invisible | biotype_switch | 0.476 | 0.527 | 0.90 | 1 |
| go_invisible | cds_changed | 0.762 | 0.734 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.762 | 0.734 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.348 | 0.394 | 0.88 | 1 |
| go_invisible | first_exon_changed | 0.971 | 0.977 | 0.99 | 0.998 |
| go_invisible | internal_exon_difference | 0.927 | 0.932 | 1.00 | 0.929 |
| go_invisible | last_exon_changed | 0.934 | 0.957 | 0.98 | 1 |
| go_invisible | nmd_switch | 0.158 | 0.158 | 1.00 | 0.513 |
| go_invisible | utr_changed | 0.464 | 0.368 | 1.26 | 0.000999 |
| go_visible | biotype_switch | 0.499 | 0.563 | 0.88 | 0.998 |
| go_visible | cds_changed | 0.752 | 0.778 | 0.97 | 0.898 |
| go_visible | coding_consequence | 0.752 | 0.778 | 0.97 | 0.898 |
| go_visible | coding_status_change | 0.402 | 0.410 | 0.98 | 0.657 |
| go_visible | first_exon_changed | 0.954 | 0.965 | 0.99 | 0.897 |
| go_visible | internal_exon_difference | 0.906 | 0.920 | 0.98 | 0.899 |
| go_visible | last_exon_changed | 0.863 | 0.909 | 0.95 | 1 |
| go_visible | nmd_switch | 0.088 | 0.188 | 0.47 | 1 |
| go_visible | utr_changed | 0.427 | 0.425 | 1.01 | 0.485 |
