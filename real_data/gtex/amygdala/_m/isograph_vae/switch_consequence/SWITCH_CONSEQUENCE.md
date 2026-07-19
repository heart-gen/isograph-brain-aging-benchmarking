# Switch coding-consequence enrichment — gtex_amygdala

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **18395** (GO-invisible 7508, GO-visible 10887)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.419 | 0.497 | 0.84 | 1 |
| all | cds_changed | 0.756 | 0.714 | 1.06 | 0.000999 |
| all | coding_consequence | 0.756 | 0.714 | 1.06 | 0.000999 |
| all | coding_status_change | 0.294 | 0.359 | 0.82 | 1 |
| all | first_exon_changed | 0.969 | 0.971 | 1.00 | 0.978 |
| all | internal_exon_difference | 0.922 | 0.938 | 0.98 | 1 |
| all | last_exon_changed | 0.911 | 0.960 | 0.95 | 1 |
| all | nmd_switch | 0.143 | 0.155 | 0.92 | 1 |
| all | utr_changed | 0.504 | 0.378 | 1.33 | 0.000999 |
| go_invisible | biotype_switch | 0.445 | 0.534 | 0.83 | 1 |
| go_invisible | cds_changed | 0.795 | 0.746 | 1.07 | 0.000999 |
| go_invisible | coding_consequence | 0.795 | 0.746 | 1.07 | 0.000999 |
| go_invisible | coding_status_change | 0.307 | 0.383 | 0.80 | 1 |
| go_invisible | first_exon_changed | 0.968 | 0.973 | 1.00 | 0.995 |
| go_invisible | internal_exon_difference | 0.929 | 0.943 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.917 | 0.962 | 0.95 | 1 |
| go_invisible | nmd_switch | 0.158 | 0.169 | 0.94 | 0.999 |
| go_invisible | utr_changed | 0.531 | 0.385 | 1.38 | 0.000999 |
| go_visible | biotype_switch | 0.402 | 0.472 | 0.85 | 1 |
| go_visible | cds_changed | 0.729 | 0.692 | 1.05 | 0.000999 |
| go_visible | coding_consequence | 0.729 | 0.692 | 1.05 | 0.000999 |
| go_visible | coding_status_change | 0.285 | 0.342 | 0.83 | 1 |
| go_visible | first_exon_changed | 0.969 | 0.970 | 1.00 | 0.738 |
| go_visible | internal_exon_difference | 0.917 | 0.935 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.906 | 0.958 | 0.95 | 1 |
| go_visible | nmd_switch | 0.133 | 0.147 | 0.91 | 1 |
| go_visible | utr_changed | 0.484 | 0.372 | 1.30 | 0.000999 |
