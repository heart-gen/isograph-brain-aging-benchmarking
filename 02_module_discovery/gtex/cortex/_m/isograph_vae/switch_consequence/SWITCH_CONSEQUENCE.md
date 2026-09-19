# Switch coding-consequence enrichment — gtex_cortex

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **36748** (GO-invisible 11951, GO-visible 24797)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.483 | 0.545 | 0.89 | 1 |
| all | cds_changed | 0.793 | 0.761 | 1.04 | 0.000999 |
| all | coding_consequence | 0.793 | 0.761 | 1.04 | 0.000999 |
| all | coding_status_change | 0.346 | 0.396 | 0.87 | 1 |
| all | first_exon_changed | 0.965 | 0.972 | 0.99 | 1 |
| all | internal_exon_difference | 0.920 | 0.934 | 0.98 | 1 |
| all | last_exon_changed | 0.907 | 0.954 | 0.95 | 1 |
| all | nmd_switch | 0.166 | 0.177 | 0.94 | 1 |
| all | utr_changed | 0.499 | 0.394 | 1.27 | 0.000999 |
| go_invisible | biotype_switch | 0.458 | 0.529 | 0.87 | 1 |
| go_invisible | cds_changed | 0.771 | 0.745 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.771 | 0.745 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.328 | 0.382 | 0.86 | 1 |
| go_invisible | first_exon_changed | 0.965 | 0.973 | 0.99 | 1 |
| go_invisible | internal_exon_difference | 0.914 | 0.935 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.883 | 0.944 | 0.94 | 1 |
| go_invisible | nmd_switch | 0.164 | 0.184 | 0.89 | 1 |
| go_invisible | utr_changed | 0.502 | 0.395 | 1.27 | 0.000999 |
| go_visible | biotype_switch | 0.494 | 0.552 | 0.89 | 1 |
| go_visible | cds_changed | 0.803 | 0.768 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.803 | 0.768 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.354 | 0.403 | 0.88 | 1 |
| go_visible | first_exon_changed | 0.965 | 0.972 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.922 | 0.934 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.918 | 0.959 | 0.96 | 1 |
| go_visible | nmd_switch | 0.167 | 0.174 | 0.96 | 1 |
| go_visible | utr_changed | 0.497 | 0.394 | 1.26 | 0.000999 |
