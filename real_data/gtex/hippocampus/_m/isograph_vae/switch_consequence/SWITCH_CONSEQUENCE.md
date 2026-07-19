# Switch coding-consequence enrichment — gtex_hippocampus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **35463** (GO-invisible 13521, GO-visible 21942)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.420 | 0.479 | 0.88 | 1 |
| all | cds_changed | 0.707 | 0.681 | 1.04 | 0.000999 |
| all | coding_consequence | 0.707 | 0.681 | 1.04 | 0.000999 |
| all | coding_status_change | 0.302 | 0.350 | 0.86 | 1 |
| all | first_exon_changed | 0.968 | 0.971 | 1.00 | 1 |
| all | internal_exon_difference | 0.903 | 0.924 | 0.98 | 1 |
| all | last_exon_changed | 0.913 | 0.954 | 0.96 | 1 |
| all | nmd_switch | 0.138 | 0.151 | 0.91 | 1 |
| all | utr_changed | 0.452 | 0.356 | 1.27 | 0.000999 |
| go_invisible | biotype_switch | 0.448 | 0.514 | 0.87 | 1 |
| go_invisible | cds_changed | 0.736 | 0.711 | 1.04 | 0.000999 |
| go_invisible | coding_consequence | 0.736 | 0.711 | 1.04 | 0.000999 |
| go_invisible | coding_status_change | 0.310 | 0.368 | 0.84 | 1 |
| go_invisible | first_exon_changed | 0.969 | 0.973 | 1.00 | 1 |
| go_invisible | internal_exon_difference | 0.898 | 0.922 | 0.97 | 1 |
| go_invisible | last_exon_changed | 0.911 | 0.953 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.159 | 0.173 | 0.92 | 1 |
| go_invisible | utr_changed | 0.477 | 0.371 | 1.29 | 0.000999 |
| go_visible | biotype_switch | 0.403 | 0.458 | 0.88 | 1 |
| go_visible | cds_changed | 0.690 | 0.662 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.690 | 0.662 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.296 | 0.339 | 0.87 | 1 |
| go_visible | first_exon_changed | 0.967 | 0.970 | 1.00 | 1 |
| go_visible | internal_exon_difference | 0.906 | 0.925 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.914 | 0.955 | 0.96 | 1 |
| go_visible | nmd_switch | 0.125 | 0.137 | 0.91 | 1 |
| go_visible | utr_changed | 0.436 | 0.346 | 1.26 | 0.000999 |
