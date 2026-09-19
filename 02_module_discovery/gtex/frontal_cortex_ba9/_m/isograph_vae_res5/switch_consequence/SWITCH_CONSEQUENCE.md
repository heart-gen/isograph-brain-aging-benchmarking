# Switch coding-consequence enrichment — gtex_frontal_cortex_ba9

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **56355** (GO-invisible 25925, GO-visible 30430)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.493 | 0.549 | 0.90 | 1 |
| all | cds_changed | 0.792 | 0.757 | 1.05 | 0.000999 |
| all | coding_consequence | 0.792 | 0.757 | 1.05 | 0.000999 |
| all | coding_status_change | 0.356 | 0.401 | 0.89 | 1 |
| all | first_exon_changed | 0.968 | 0.973 | 1.00 | 1 |
| all | internal_exon_difference | 0.918 | 0.930 | 0.99 | 1 |
| all | last_exon_changed | 0.923 | 0.960 | 0.96 | 1 |
| all | nmd_switch | 0.163 | 0.174 | 0.94 | 1 |
| all | utr_changed | 0.482 | 0.383 | 1.26 | 0.000999 |
| go_invisible | biotype_switch | 0.492 | 0.545 | 0.90 | 1 |
| go_invisible | cds_changed | 0.786 | 0.745 | 1.05 | 0.000999 |
| go_invisible | coding_consequence | 0.786 | 0.745 | 1.05 | 0.000999 |
| go_invisible | coding_status_change | 0.356 | 0.398 | 0.90 | 1 |
| go_invisible | first_exon_changed | 0.971 | 0.973 | 1.00 | 0.993 |
| go_invisible | internal_exon_difference | 0.921 | 0.930 | 0.99 | 1 |
| go_invisible | last_exon_changed | 0.920 | 0.959 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.162 | 0.174 | 0.93 | 1 |
| go_invisible | utr_changed | 0.473 | 0.373 | 1.27 | 0.000999 |
| go_visible | biotype_switch | 0.494 | 0.551 | 0.90 | 1 |
| go_visible | cds_changed | 0.798 | 0.766 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.798 | 0.766 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.356 | 0.403 | 0.88 | 1 |
| go_visible | first_exon_changed | 0.966 | 0.973 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.916 | 0.930 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.926 | 0.962 | 0.96 | 1 |
| go_visible | nmd_switch | 0.163 | 0.173 | 0.94 | 1 |
| go_visible | utr_changed | 0.490 | 0.391 | 1.25 | 0.000999 |
