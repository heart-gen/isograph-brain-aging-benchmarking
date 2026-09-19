# Switch coding-consequence enrichment — brainseq_caudate

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **9303** (GO-invisible 0, GO-visible 9303)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.484 | 0.535 | 0.91 | 1 |
| all | cds_changed | 0.738 | 0.714 | 1.03 | 0.000999 |
| all | coding_consequence | 0.738 | 0.714 | 1.03 | 0.000999 |
| all | coding_status_change | 0.355 | 0.389 | 0.91 | 1 |
| all | first_exon_changed | 0.968 | 0.975 | 0.99 | 1 |
| all | internal_exon_difference | 0.907 | 0.921 | 0.98 | 1 |
| all | last_exon_changed | 0.938 | 0.961 | 0.98 | 1 |
| all | nmd_switch | 0.154 | 0.171 | 0.90 | 1 |
| all | utr_changed | 0.437 | 0.349 | 1.25 | 0.000999 |
| go_visible | biotype_switch | 0.484 | 0.534 | 0.91 | 1 |
| go_visible | cds_changed | 0.738 | 0.713 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.738 | 0.713 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.355 | 0.389 | 0.91 | 1 |
| go_visible | first_exon_changed | 0.968 | 0.975 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.907 | 0.921 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.938 | 0.961 | 0.98 | 1 |
| go_visible | nmd_switch | 0.154 | 0.171 | 0.90 | 1 |
| go_visible | utr_changed | 0.437 | 0.349 | 1.25 | 0.000999 |
