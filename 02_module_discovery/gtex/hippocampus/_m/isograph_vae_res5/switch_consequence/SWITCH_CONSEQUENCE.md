# Switch coding-consequence enrichment — gtex_hippocampus

Per consequence class: fraction of IsoGraph switch pairs with the change (`obs_rate`) vs a within-gene permutation null of random transcript pairs (`null_mean`); `enrichment` = obs/null, `p_emp` = one-sided permutation p. `coding_consequence` = CDS change or coding-status change. Stratified by whether the switch gene's phenotype-significant module is GO-invisible.

- switch pairs scored: **18503** (GO-invisible 552, GO-visible 17951)

| stratum | consequence | obs rate | null mean | enrichment | p |
|---------|-------------|----------|-----------|------------|---|
| all | biotype_switch | 0.502 | 0.550 | 0.91 | 1 |
| all | cds_changed | 0.788 | 0.760 | 1.04 | 0.000999 |
| all | coding_consequence | 0.788 | 0.760 | 1.04 | 0.000999 |
| all | coding_status_change | 0.366 | 0.404 | 0.91 | 1 |
| all | first_exon_changed | 0.966 | 0.972 | 0.99 | 1 |
| all | internal_exon_difference | 0.914 | 0.927 | 0.99 | 1 |
| all | last_exon_changed | 0.916 | 0.955 | 0.96 | 1 |
| all | nmd_switch | 0.158 | 0.172 | 0.92 | 1 |
| all | utr_changed | 0.469 | 0.385 | 1.22 | 0.000999 |
| go_invisible | biotype_switch | 0.469 | 0.531 | 0.88 | 0.999 |
| go_invisible | cds_changed | 0.801 | 0.793 | 1.01 | 0.322 |
| go_invisible | coding_consequence | 0.801 | 0.793 | 1.01 | 0.322 |
| go_invisible | coding_status_change | 0.362 | 0.400 | 0.91 | 0.974 |
| go_invisible | first_exon_changed | 0.953 | 0.966 | 0.99 | 0.966 |
| go_invisible | internal_exon_difference | 0.918 | 0.926 | 0.99 | 0.817 |
| go_invisible | last_exon_changed | 0.861 | 0.916 | 0.94 | 1 |
| go_invisible | nmd_switch | 0.125 | 0.144 | 0.87 | 0.934 |
| go_invisible | utr_changed | 0.489 | 0.423 | 1.16 | 0.000999 |
| go_visible | biotype_switch | 0.503 | 0.550 | 0.91 | 1 |
| go_visible | cds_changed | 0.787 | 0.759 | 1.04 | 0.000999 |
| go_visible | coding_consequence | 0.787 | 0.759 | 1.04 | 0.000999 |
| go_visible | coding_status_change | 0.367 | 0.404 | 0.91 | 1 |
| go_visible | first_exon_changed | 0.966 | 0.972 | 0.99 | 1 |
| go_visible | internal_exon_difference | 0.914 | 0.927 | 0.99 | 1 |
| go_visible | last_exon_changed | 0.918 | 0.956 | 0.96 | 1 |
| go_visible | nmd_switch | 0.159 | 0.173 | 0.92 | 1 |
| go_visible | utr_changed | 0.468 | 0.384 | 1.22 | 0.000999 |
