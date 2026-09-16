# Switch coding-consequence enrichment — gtex_hypothalamus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **5184** (GO-invisible 5184, GO-visible 0)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.444 | 0.545 | 0.81 | 1 |
| all | cds_changed | 0.807 | 0.783 | 1.03 | 0.000999 |
| all | coding_consequence | 0.807 | 0.783 | 1.03 | 0.000999 |
| all | coding_status_change | 0.293 | 0.386 | 0.76 | 1 |
| all | first_exon_changed | 0.957 | 0.964 | 0.99 | 0.999 |
| all | internal_exon_difference | 0.915 | 0.937 | 0.98 | 1 |
| all | last_exon_changed | 0.868 | 0.931 | 0.93 | 1 |
| all | nmd_switch | 0.186 | 0.206 | 0.91 | 1 |
| all | utr_changed | 0.591 | 0.437 | 1.35 | 0.000999 |
| go_invisible | biotype_switch | 0.444 | 0.545 | 0.81 | 1 |
| go_invisible | cds_changed | 0.807 | 0.782 | 1.03 | 0.000999 |
| go_invisible | coding_consequence | 0.807 | 0.782 | 1.03 | 0.000999 |
| go_invisible | coding_status_change | 0.293 | 0.386 | 0.76 | 1 |
| go_invisible | first_exon_changed | 0.957 | 0.964 | 0.99 | 0.996 |
| go_invisible | internal_exon_difference | 0.915 | 0.937 | 0.98 | 1 |
| go_invisible | last_exon_changed | 0.868 | 0.931 | 0.93 | 1 |
| go_invisible | nmd_switch | 0.186 | 0.206 | 0.91 | 1 |
| go_invisible | utr_changed | 0.591 | 0.436 | 1.35 | 0.000999 |
