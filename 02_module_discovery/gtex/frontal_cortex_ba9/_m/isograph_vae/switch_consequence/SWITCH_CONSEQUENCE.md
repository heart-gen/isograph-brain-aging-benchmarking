# Switch coding-consequence enrichment — gtex_frontal_cortex_ba9

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **61434** (GO-invisible 9290, GO-visible 52144)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.479 | 0.536 | 0.89 | 1 |
| all | cds_changed | 0.781 | 0.748 | 1.04 | 0.000999 |
| all | coding_consequence | 0.781 | 0.748 | 1.04 | 0.000999 |
| all | coding_status_change | 0.345 | 0.392 | 0.88 | 1 |
| all | first_exon_changed | 0.967 | 0.972 | 0.99 | 1 |
| all | internal_exon_difference | 0.913 | 0.927 | 0.98 | 1 |
| all | last_exon_changed | 0.917 | 0.957 | 0.96 | 1 |
| all | nmd_switch | 0.160 | 0.170 | 0.94 | 1 |
| all | utr_changed | 0.486 | 0.385 | 1.26 | 0.000999 |
| go_invisible | biotype_switch | 0.409 | 0.477 | 0.86 | 1 |
| go_invisible | cds_changed | 0.722 | 0.692 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.722 | 0.692 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.280 | 0.344 | 0.81 | 1 |
| go_invisible | first_exon_changed | 0.969 | 0.971 | 1.00 | 0.93 |
| go_invisible | internal_exon_difference | 0.905 | 0.926 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.887 | 0.940 | 0.94 | 1 |
| go_invisible | nmd_switch | 0.155 | 0.165 | 0.94 | 0.999 |
| go_invisible | utr_changed | 0.500 | 0.384 | 1.30 | 0.000999 |
| go_visible | biotype_switch | 0.491 | 0.546 | 0.90 | 1 |
| go_visible | cds_changed | 0.792 | 0.758 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.792 | 0.758 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.356 | 0.401 | 0.89 | 1 |
| go_visible | first_exon_changed | 0.967 | 0.972 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.915 | 0.928 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.922 | 0.960 | 0.96 | 1 |
| go_visible | nmd_switch | 0.161 | 0.171 | 0.94 | 1 |
| go_visible | utr_changed | 0.483 | 0.385 | 1.26 | 0.000999 |
