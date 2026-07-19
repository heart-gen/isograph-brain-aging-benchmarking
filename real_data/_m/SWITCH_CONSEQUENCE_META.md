# Switch coding-consequence enrichment — cross-region rollup

Per consequence class and GO-stratum, aggregated over the switch-layer regions: how many regions show enrichment > 1 at empirical p < 0.05 (vs the within-gene random-pair null), how many show depletion, the median enrichment, and a Fisher-combined p. `coding_consequence` = CDS change or coding-status change.

| stratum | consequence | regions | enriched (p<.05) | depleted | median enrich | median obs rate | Fisher p |
|---------|-------------|---------|------------------|----------|---------------|-----------------|----------|
| all | biotype_switch | 10 | 0 | 0 | 0.87 | 0.449 | 1.00e+00 |
| all | cds_changed | 10 | 10 | 0 | 1.04 | 0.756 | 2.05e-19 |
| all | coding_consequence | 10 | 10 | 0 | 1.04 | 0.756 | 2.05e-19 |
| all | coding_status_change | 10 | 0 | 0 | 0.86 | 0.312 | 1.00e+00 |
| all | first_exon_changed | 10 | 0 | 0 | 0.99 | 0.965 | 1.00e+00 |
| all | internal_exon_difference | 10 | 0 | 0 | 0.98 | 0.912 | 1.00e+00 |
| all | last_exon_changed | 10 | 0 | 0 | 0.96 | 0.918 | 1.00e+00 |
| all | nmd_switch | 10 | 0 | 0 | 0.92 | 0.151 | 1.00e+00 |
| all | utr_changed | 10 | 10 | 0 | 1.27 | 0.471 | 1.12e-19 |
| go_invisible | biotype_switch | 10 | 0 | 0 | 0.87 | 0.453 | 1.00e+00 |
| go_invisible | cds_changed | 10 | 10 | 0 | 1.04 | 0.764 | 1.73e-17 |
| go_invisible | coding_consequence | 10 | 10 | 0 | 1.04 | 0.764 | 1.73e-17 |
| go_invisible | coding_status_change | 10 | 0 | 0 | 0.84 | 0.315 | 1.00e+00 |
| go_invisible | first_exon_changed | 10 | 0 | 0 | 0.99 | 0.964 | 1.00e+00 |
| go_invisible | internal_exon_difference | 10 | 0 | 0 | 0.98 | 0.910 | 1.00e+00 |
| go_invisible | last_exon_changed | 10 | 0 | 0 | 0.96 | 0.915 | 1.00e+00 |
| go_invisible | nmd_switch | 10 | 0 | 0 | 0.92 | 0.161 | 1.00e+00 |
| go_invisible | utr_changed | 10 | 10 | 0 | 1.29 | 0.491 | 1.12e-19 |
| go_visible | biotype_switch | 9 | 0 | 0 | 0.89 | 0.421 | 1.00e+00 |
| go_visible | cds_changed | 9 | 9 | 0 | 1.04 | 0.725 | 7.79e-17 |
| go_visible | coding_consequence | 9 | 9 | 0 | 1.04 | 0.725 | 7.79e-17 |
| go_visible | coding_status_change | 9 | 0 | 0 | 0.87 | 0.296 | 1.00e+00 |
| go_visible | first_exon_changed | 9 | 0 | 0 | 1.00 | 0.967 | 1.00e+00 |
| go_visible | internal_exon_difference | 9 | 0 | 0 | 0.98 | 0.907 | 1.00e+00 |
| go_visible | last_exon_changed | 9 | 0 | 0 | 0.96 | 0.914 | 1.00e+00 |
| go_visible | nmd_switch | 9 | 0 | 0 | 0.95 | 0.137 | 1.00e+00 |
| go_visible | utr_changed | 9 | 9 | 0 | 1.26 | 0.481 | 6.29e-18 |
