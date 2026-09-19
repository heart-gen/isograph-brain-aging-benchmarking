# Switch coding-consequence enrichment — gtex_amygdala

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **811** (GO-invisible 0, GO-visible 811)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.475 | 0.517 | 0.92 | 0.998 |
| all | cds_changed | 0.776 | 0.759 | 1.02 | 0.105 |
| all | coding_consequence | 0.776 | 0.759 | 1.02 | 0.105 |
| all | coding_status_change | 0.340 | 0.379 | 0.90 | 0.994 |
| all | first_exon_changed | 0.957 | 0.960 | 1.00 | 0.73 |
| all | internal_exon_difference | 0.930 | 0.934 | 1.00 | 0.732 |
| all | last_exon_changed | 0.871 | 0.934 | 0.93 | 1 |
| all | nmd_switch | 0.164 | 0.165 | 0.99 | 0.553 |
| all | utr_changed | 0.522 | 0.441 | 1.18 | 0.000999 |
| go_visible | biotype_switch | 0.475 | 0.516 | 0.92 | 0.999 |
| go_visible | cds_changed | 0.776 | 0.759 | 1.02 | 0.0979 |
| go_visible | coding_consequence | 0.776 | 0.759 | 1.02 | 0.0979 |
| go_visible | coding_status_change | 0.340 | 0.379 | 0.90 | 0.996 |
| go_visible | first_exon_changed | 0.957 | 0.960 | 1.00 | 0.726 |
| go_visible | internal_exon_difference | 0.930 | 0.934 | 1.00 | 0.734 |
| go_visible | last_exon_changed | 0.871 | 0.934 | 0.93 | 1 |
| go_visible | nmd_switch | 0.164 | 0.164 | 1.00 | 0.553 |
| go_visible | utr_changed | 0.522 | 0.441 | 1.18 | 0.000999 |
