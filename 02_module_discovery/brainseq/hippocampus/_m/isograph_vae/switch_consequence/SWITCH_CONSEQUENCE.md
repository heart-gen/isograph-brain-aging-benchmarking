# Switch coding-consequence enrichment — brainseq_hippocampus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **995** (GO-invisible 995, GO-visible 0)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.432 | 0.562 | 0.77 | 1 |
| all | cds_changed | 0.863 | 0.829 | 1.04 | 0.000999 |
| all | coding_consequence | 0.863 | 0.829 | 1.04 | 0.000999 |
| all | coding_status_change | 0.286 | 0.389 | 0.74 | 1 |
| all | first_exon_changed | 0.970 | 0.977 | 0.99 | 0.951 |
| all | internal_exon_difference | 0.898 | 0.933 | 0.96 | 1 |
| all | last_exon_changed | 0.914 | 0.955 | 0.96 | 1 |
| all | nmd_switch | 0.185 | 0.207 | 0.89 | 0.975 |
| all | utr_changed | 0.652 | 0.474 | 1.38 | 0.000999 |
| go_invisible | biotype_switch | 0.432 | 0.563 | 0.77 | 1 |
| go_invisible | cds_changed | 0.863 | 0.829 | 1.04 | 0.003 |
| go_invisible | coding_consequence | 0.863 | 0.829 | 1.04 | 0.003 |
| go_invisible | coding_status_change | 0.286 | 0.390 | 0.73 | 1 |
| go_invisible | first_exon_changed | 0.970 | 0.977 | 0.99 | 0.949 |
| go_invisible | internal_exon_difference | 0.898 | 0.934 | 0.96 | 1 |
| go_invisible | last_exon_changed | 0.914 | 0.955 | 0.96 | 1 |
| go_invisible | nmd_switch | 0.185 | 0.208 | 0.89 | 0.974 |
| go_invisible | utr_changed | 0.652 | 0.473 | 1.38 | 0.000999 |
