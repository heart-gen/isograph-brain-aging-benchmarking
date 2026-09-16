# Switch coding-consequence enrichment — cross-region rollup

Per consequence class and GO-stratum, pooled separately for the aging and disease analyses. The primary estimate is the **random-effects pooled enrichment ratio** (DerSimonian-Laird) over the per-analysis log enrichments, with a 95% CI and I² for between-analysis heterogeneity. The Fisher-combined p is secondary and should not be read as an effect size: it grows more significant with every added analysis regardless of magnitude. `coding_consequence` = CDS change or coding-status change.

| class | stratum | consequence | k | pooled ratio (RE) | 95% CI | p (RE) | I² | median obs rate | Fisher p |
|-------|---------|-------------|---|-------------------|--------|--------|----|-----------------|----------|
| aging | all | biotype_switch | 9 | 0.900 | 0.89–0.91 | 2.81e-45 | 0.80 | 0.498 | 1.00e+00 |
| aging | all | cds_changed | 9 | 1.044 | 1.04–1.05 | 2.44e-163 | 0.00 | 0.792 | 1.15e-17 |
| aging | all | coding_consequence | 9 | 1.044 | 1.04–1.05 | 2.44e-163 | 0.00 | 0.792 | 1.15e-17 |
| aging | all | coding_status_change | 9 | 0.890 | 0.87–0.91 | 2.07e-24 | 0.84 | 0.359 | 1.00e+00 |
| aging | all | first_exon_changed | 9 | 0.995 | 0.99–1.00 | 1.55e-23 | 0.00 | 0.968 | 1.00e+00 |
| aging | all | internal_exon_difference | 9 | 0.987 | 0.99–0.99 | 1.63e-54 | 0.11 | 0.918 | 1.00e+00 |
| aging | all | last_exon_changed | 9 | 0.959 | 0.95–0.97 | 2.26e-32 | 0.92 | 0.919 | 1.00e+00 |
| aging | all | nmd_switch | 9 | 0.935 | 0.92–0.95 | 2.24e-20 | 0.00 | 0.165 | 1.00e+00 |
| aging | all | utr_changed | 9 | 1.248 | 1.23–1.27 | 5.50e-161 | 0.70 | 0.479 | 6.29e-18 |
| aging | go_invisible | biotype_switch | 9 | 0.899 | 0.87–0.93 | 9.92e-12 | 0.85 | 0.470 | 1.00e+00 |
| aging | go_invisible | cds_changed | 9 | 1.040 | 1.03–1.05 | 9.08e-18 | 0.50 | 0.771 | 2.26e-14 |
| aging | go_invisible | coding_consequence | 9 | 1.040 | 1.03–1.05 | 9.08e-18 | 0.50 | 0.771 | 2.26e-14 |
| aging | go_invisible | coding_status_change | 9 | 0.894 | 0.86–0.93 | 6.13e-07 | 0.87 | 0.347 | 1.00e+00 |
| aging | go_invisible | first_exon_changed | 9 | 0.995 | 0.99–1.00 | 2.11e-05 | 0.37 | 0.968 | 1.00e+00 |
| aging | go_invisible | internal_exon_difference | 9 | 0.986 | 0.98–0.99 | 1.48e-14 | 0.29 | 0.915 | 1.00e+00 |
| aging | go_invisible | last_exon_changed | 9 | 0.957 | 0.95–0.97 | 2.95e-12 | 0.92 | 0.918 | 1.00e+00 |
| aging | go_invisible | nmd_switch | 9 | 0.920 | 0.90–0.94 | 1.49e-10 | 0.00 | 0.162 | 1.00e+00 |
| aging | go_invisible | utr_changed | 9 | 1.233 | 1.19–1.28 | 3.05e-29 | 0.80 | 0.473 | 6.29e-18 |
| aging | go_visible | biotype_switch | 7 | 0.905 | 0.90–0.91 | 2.85e-106 | 0.31 | 0.506 | 1.00e+00 |
| aging | go_visible | cds_changed | 7 | 1.044 | 1.04–1.05 | 1.55e-113 | 0.00 | 0.808 | 3.12e-13 |
| aging | go_visible | coding_consequence | 7 | 1.044 | 1.04–1.05 | 1.55e-113 | 0.00 | 0.808 | 3.12e-13 |
| aging | go_visible | coding_status_change | 7 | 0.897 | 0.88–0.91 | 1.92e-39 | 0.60 | 0.367 | 1.00e+00 |
| aging | go_visible | first_exon_changed | 7 | 0.995 | 0.99–1.00 | 1.28e-16 | 0.00 | 0.966 | 1.00e+00 |
| aging | go_visible | internal_exon_difference | 7 | 0.987 | 0.99–0.99 | 5.69e-38 | 0.10 | 0.924 | 1.00e+00 |
| aging | go_visible | last_exon_changed | 7 | 0.959 | 0.96–0.96 | 5.88e-254 | 0.12 | 0.920 | 1.00e+00 |
| aging | go_visible | nmd_switch | 7 | 0.942 | 0.93–0.96 | 1.19e-11 | 0.00 | 0.165 | 1.00e+00 |
| aging | go_visible | utr_changed | 7 | 1.246 | 1.23–1.26 | 0.00e+00 | 0.00 | 0.482 | 2.01e-14 |
| disease | all | biotype_switch | 1 | 0.902 | 0.87–0.94 | 6.89e-08 | 0.00 | 0.477 | 1.00e+00 |
| disease | all | cds_changed | 1 | 1.034 | 1.01–1.06 | 2.31e-03 | 0.00 | 0.762 | 9.99e-04 |
| disease | all | coding_consequence | 1 | 1.034 | 1.01–1.06 | 2.31e-03 | 0.00 | 0.762 | 9.99e-04 |
| disease | all | coding_status_change | 1 | 0.888 | 0.84–0.94 | 9.50e-06 | 0.00 | 0.351 | 1.00e+00 |
| disease | all | first_exon_changed | 1 | 0.994 | 0.99–1.00 | 3.72e-02 | 0.00 | 0.970 | 9.97e-01 |
| disease | all | internal_exon_difference | 1 | 0.995 | 0.99–1.00 | 2.20e-01 | 0.00 | 0.926 | 9.56e-01 |
| disease | all | last_exon_changed | 1 | 0.974 | 0.96–0.98 | 1.91e-06 | 0.00 | 0.929 | 1.00e+00 |
| disease | all | nmd_switch | 1 | 0.959 | 0.87–1.06 | 3.97e-01 | 0.00 | 0.153 | 9.37e-01 |
| disease | all | utr_changed | 1 | 1.244 | 1.18–1.31 | 7.26e-18 | 0.00 | 0.462 | 9.99e-04 |
| disease | go_invisible | biotype_switch | 1 | 0.903 | 0.87–0.94 | 1.95e-07 | 0.00 | 0.476 | 1.00e+00 |
| disease | go_invisible | cds_changed | 1 | 1.039 | 1.02–1.06 | 1.21e-03 | 0.00 | 0.762 | 9.99e-04 |
| disease | go_invisible | coding_consequence | 1 | 1.039 | 1.02–1.06 | 1.21e-03 | 0.00 | 0.762 | 9.99e-04 |
| disease | go_invisible | coding_status_change | 1 | 0.881 | 0.83–0.93 | 6.82e-06 | 0.00 | 0.348 | 1.00e+00 |
| disease | go_invisible | first_exon_changed | 1 | 0.994 | 0.99–1.00 | 4.55e-02 | 0.00 | 0.971 | 9.98e-01 |
| disease | go_invisible | internal_exon_difference | 1 | 0.995 | 0.99–1.00 | 3.03e-01 | 0.00 | 0.927 | 9.29e-01 |
| disease | go_invisible | last_exon_changed | 1 | 0.976 | 0.96–0.99 | 1.04e-05 | 0.00 | 0.934 | 1.00e+00 |
| disease | go_invisible | nmd_switch | 1 | 1.000 | 0.91–1.10 | 9.92e-01 | 0.00 | 0.158 | 5.13e-01 |
| disease | go_invisible | utr_changed | 1 | 1.263 | 1.20–1.33 | 1.27e-17 | 0.00 | 0.464 | 9.99e-04 |
| disease | go_visible | biotype_switch | 1 | 0.885 | 0.76–1.02 | 1.02e-01 | 0.00 | 0.499 | 9.98e-01 |
| disease | go_visible | cds_changed | 1 | 0.967 | 0.88–1.06 | 4.62e-01 | 0.00 | 0.752 | 8.98e-01 |
| disease | go_visible | coding_consequence | 1 | 0.967 | 0.88–1.06 | 4.62e-01 | 0.00 | 0.752 | 8.98e-01 |
| disease | go_visible | coding_status_change | 1 | 0.979 | 0.82–1.17 | 8.14e-01 | 0.00 | 0.402 | 6.57e-01 |
| disease | go_visible | first_exon_changed | 1 | 0.989 | 0.96–1.01 | 4.06e-01 | 0.00 | 0.954 | 8.97e-01 |
| disease | go_visible | internal_exon_difference | 1 | 0.985 | 0.95–1.02 | 3.82e-01 | 0.00 | 0.906 | 8.99e-01 |
| disease | go_visible | last_exon_changed | 1 | 0.950 | 0.89–1.01 | 1.15e-01 | 0.00 | 0.863 | 1.00e+00 |
| disease | go_visible | nmd_switch | 1 | 0.469 | 0.24–0.93 | 3.13e-02 | 0.00 | 0.088 | 1.00e+00 |
| disease | go_visible | utr_changed | 1 | 1.007 | 0.83–1.22 | 9.46e-01 | 0.00 | 0.427 | 4.85e-01 |
