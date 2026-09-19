# Switch coding-consequence enrichment — gtex_anterior_cingulate_cortex_ba24

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **34132** (GO-invisible 0, GO-visible 34132)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.490 | 0.548 | 0.89 | 1 |
| all | cds_changed | 0.798 | 0.767 | 1.04 | 0.000999 |
| all | coding_consequence | 0.798 | 0.767 | 1.04 | 0.000999 |
| all | coding_status_change | 0.350 | 0.399 | 0.88 | 1 |
| all | first_exon_changed | 0.966 | 0.971 | 0.99 | 1 |
| all | internal_exon_difference | 0.913 | 0.929 | 0.98 | 1 |
| all | last_exon_changed | 0.920 | 0.960 | 0.96 | 1 |
| all | nmd_switch | 0.165 | 0.177 | 0.93 | 1 |
| all | utr_changed | 0.498 | 0.398 | 1.25 | 0.000999 |
| go_visible | biotype_switch | 0.490 | 0.548 | 0.89 | 1 |
| go_visible | cds_changed | 0.798 | 0.767 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.798 | 0.767 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.350 | 0.399 | 0.88 | 1 |
| go_visible | first_exon_changed | 0.966 | 0.971 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.913 | 0.929 | 0.98 | 1 |
| go_visible | last_exon_changed | 0.920 | 0.960 | 0.96 | 1 |
| go_visible | nmd_switch | 0.165 | 0.177 | 0.93 | 1 |
| go_visible | utr_changed | 0.498 | 0.397 | 1.25 | 0.000999 |
