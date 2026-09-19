# Switch coding-consequence enrichment — gtex_anterior_cingulate_cortex_ba24

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **35958** (GO-invisible 5980, GO-visible 29978)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.500 | 0.555 | 0.90 | 1 |
| all | cds_changed | 0.795 | 0.763 | 1.04 | 0.000999 |
| all | coding_consequence | 0.795 | 0.763 | 1.04 | 0.000999 |
| all | coding_status_change | 0.360 | 0.404 | 0.89 | 1 |
| all | first_exon_changed | 0.967 | 0.971 | 1.00 | 1 |
| all | internal_exon_difference | 0.914 | 0.929 | 0.98 | 1 |
| all | last_exon_changed | 0.924 | 0.961 | 0.96 | 1 |
| all | nmd_switch | 0.167 | 0.180 | 0.93 | 1 |
| all | utr_changed | 0.484 | 0.388 | 1.25 | 0.000999 |
| go_invisible | biotype_switch | 0.470 | 0.505 | 0.93 | 1 |
| go_invisible | cds_changed | 0.689 | 0.661 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.689 | 0.661 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.346 | 0.368 | 0.94 | 1 |
| go_invisible | first_exon_changed | 0.970 | 0.973 | 1.00 | 0.935 |
| go_invisible | internal_exon_difference | 0.901 | 0.918 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.924 | 0.952 | 0.97 | 1 |
| go_invisible | nmd_switch | 0.145 | 0.164 | 0.89 | 1 |
| go_invisible | utr_changed | 0.380 | 0.316 | 1.20 | 0.000999 |
| go_visible | biotype_switch | 0.506 | 0.565 | 0.90 | 1 |
| go_visible | cds_changed | 0.816 | 0.784 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.816 | 0.784 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.362 | 0.411 | 0.88 | 1 |
| go_visible | first_exon_changed | 0.966 | 0.971 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.916 | 0.931 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.924 | 0.963 | 0.96 | 1 |
| go_visible | nmd_switch | 0.171 | 0.184 | 0.93 | 1 |
| go_visible | utr_changed | 0.505 | 0.403 | 1.25 | 0.000999 |
