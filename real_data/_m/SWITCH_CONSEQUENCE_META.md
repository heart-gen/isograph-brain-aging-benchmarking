# Switch coding-consequence enrichment — cross-region rollup

Per consequence class and GO-stratum, pooled separately for the aging and disease analyses. The primary estimate is the **random-effects pooled enrichment ratio** (DerSimonian-Laird) over the per-analysis log enrichments, with a 95% CI and I² for between-analysis heterogeneity. The Fisher-combined p is secondary and should not be read as an effect size: it grows more significant with every added analysis regardless of magnitude. `coding_consequence` = CDS change or coding-status change.

| class | stratum | consequence | k | pooled ratio (RE) | 95% CI | p (RE) | I² | median obs rate | Fisher p |
|-------|---------|-------------|---|-------------------|--------|--------|----|-----------------|----------|
| aging | all | biotype_switch | 9 | 0.881 | 0.87–0.90 | 7.21e-50 | 0.80 | 0.453 | 1.00e+00 |
| aging | all | cds_changed | 9 | 1.045 | 1.04–1.05 | 1.11e-81 | 0.36 | 0.756 | 1.15e-17 |
| aging | all | coding_consequence | 9 | 1.045 | 1.04–1.05 | 1.11e-81 | 0.36 | 0.756 | 1.15e-17 |
| aging | all | coding_status_change | 9 | 0.861 | 0.84–0.88 | 1.96e-45 | 0.75 | 0.322 | 1.00e+00 |
| aging | all | first_exon_changed | 9 | 0.995 | 0.99–1.00 | 7.20e-17 | 0.12 | 0.966 | 1.00e+00 |
| aging | all | internal_exon_difference | 9 | 0.981 | 0.98–0.98 | 1.57e-59 | 0.31 | 0.913 | 1.00e+00 |
| aging | all | last_exon_changed | 9 | 0.960 | 0.95–0.97 | 2.16e-29 | 0.92 | 0.918 | 1.00e+00 |
| aging | all | nmd_switch | 9 | 0.939 | 0.92–0.95 | 4.34e-16 | 0.00 | 0.153 | 1.00e+00 |
| aging | all | utr_changed | 9 | 1.280 | 1.26–1.30 | 4.43e-224 | 0.62 | 0.483 | 6.29e-18 |
| aging | go_invisible | biotype_switch | 9 | 0.874 | 0.85–0.89 | 9.75e-32 | 0.78 | 0.453 | 1.00e+00 |
| aging | go_invisible | cds_changed | 9 | 1.047 | 1.04–1.05 | 5.53e-40 | 0.45 | 0.771 | 6.38e-16 |
| aging | go_invisible | coding_consequence | 9 | 1.047 | 1.04–1.05 | 5.53e-40 | 0.45 | 0.771 | 6.38e-16 |
| aging | go_invisible | coding_status_change | 9 | 0.850 | 0.82–0.88 | 2.22e-25 | 0.77 | 0.314 | 1.00e+00 |
| aging | go_invisible | first_exon_changed | 9 | 0.993 | 0.99–0.99 | 4.09e-18 | 0.00 | 0.964 | 1.00e+00 |
| aging | go_invisible | internal_exon_difference | 9 | 0.980 | 0.98–0.98 | 6.00e-28 | 0.48 | 0.912 | 1.00e+00 |
| aging | go_invisible | last_exon_changed | 9 | 0.958 | 0.95–0.97 | 2.88e-19 | 0.91 | 0.914 | 1.00e+00 |
| aging | go_invisible | nmd_switch | 9 | 0.931 | 0.91–0.95 | 1.45e-11 | 0.00 | 0.163 | 1.00e+00 |
| aging | go_invisible | utr_changed | 9 | 1.303 | 1.28–1.33 | 3.59e-144 | 0.59 | 0.504 | 6.29e-18 |
| aging | go_visible | biotype_switch | 8 | 0.888 | 0.88–0.90 | 9.31e-63 | 0.40 | 0.424 | 1.00e+00 |
| aging | go_visible | cds_changed | 8 | 1.041 | 1.04–1.05 | 5.88e-70 | 0.00 | 0.727 | 5.46e-15 |
| aging | go_visible | coding_consequence | 8 | 1.041 | 1.04–1.05 | 5.88e-70 | 0.00 | 0.727 | 5.46e-15 |
| aging | go_visible | coding_status_change | 8 | 0.871 | 0.86–0.88 | 1.50e-103 | 0.00 | 0.300 | 1.00e+00 |
| aging | go_visible | first_exon_changed | 8 | 0.997 | 1.00–1.00 | 2.63e-06 | 0.00 | 0.967 | 1.00e+00 |
| aging | go_visible | internal_exon_difference | 8 | 0.983 | 0.98–0.99 | 1.08e-21 | 0.40 | 0.909 | 1.00e+00 |
| aging | go_visible | last_exon_changed | 8 | 0.960 | 0.95–0.96 | 2.10e-56 | 0.63 | 0.916 | 1.00e+00 |
| aging | go_visible | nmd_switch | 8 | 0.949 | 0.93–0.97 | 2.46e-06 | 0.00 | 0.141 | 1.00e+00 |
| aging | go_visible | utr_changed | 8 | 1.255 | 1.24–1.27 | 0.00e+00 | 0.00 | 0.483 | 3.54e-16 |
| disease | all | biotype_switch | 1 | 0.898 | 0.85–0.95 | 2.23e-04 | 0.00 | 0.412 | 1.00e+00 |
| disease | all | cds_changed | 1 | 1.043 | 1.01–1.08 | 1.33e-02 | 0.00 | 0.669 | 9.99e-04 |
| disease | all | coding_consequence | 1 | 1.043 | 1.01–1.08 | 1.33e-02 | 0.00 | 0.669 | 9.99e-04 |
| disease | all | coding_status_change | 1 | 0.892 | 0.83–0.96 | 2.15e-03 | 0.00 | 0.300 | 1.00e+00 |
| disease | all | first_exon_changed | 1 | 0.988 | 0.98–1.00 | 2.78e-02 | 0.00 | 0.956 | 1.00e+00 |
| disease | all | internal_exon_difference | 1 | 0.973 | 0.96–0.99 | 1.14e-03 | 0.00 | 0.894 | 1.00e+00 |
| disease | all | last_exon_changed | 1 | 0.975 | 0.96–0.99 | 3.25e-03 | 0.00 | 0.928 | 1.00e+00 |
| disease | all | nmd_switch | 1 | 0.907 | 0.78–1.05 | 2.02e-01 | 0.00 | 0.131 | 9.95e-01 |
| disease | all | utr_changed | 1 | 1.264 | 1.18–1.36 | 1.00e-10 | 0.00 | 0.430 | 9.99e-04 |
| disease | go_invisible | biotype_switch | 1 | 0.902 | 0.84–0.97 | 4.65e-03 | 0.00 | 0.434 | 1.00e+00 |
| disease | go_invisible | cds_changed | 1 | 1.033 | 0.99–1.08 | 1.54e-01 | 0.00 | 0.670 | 9.99e-04 |
| disease | go_invisible | coding_consequence | 1 | 1.033 | 0.99–1.08 | 1.54e-01 | 0.00 | 0.670 | 9.99e-04 |
| disease | go_invisible | coding_status_change | 1 | 0.900 | 0.82–0.99 | 3.21e-02 | 0.00 | 0.318 | 1.00e+00 |
| disease | go_invisible | first_exon_changed | 1 | 0.990 | 0.98–1.00 | 1.32e-01 | 0.00 | 0.958 | 9.94e-01 |
| disease | go_invisible | internal_exon_difference | 1 | 0.967 | 0.95–0.99 | 3.49e-03 | 0.00 | 0.904 | 1.00e+00 |
| disease | go_invisible | last_exon_changed | 1 | 0.975 | 0.95–1.00 | 2.11e-02 | 0.00 | 0.938 | 1.00e+00 |
| disease | go_invisible | nmd_switch | 1 | 0.885 | 0.73–1.07 | 2.08e-01 | 0.00 | 0.136 | 9.97e-01 |
| disease | go_invisible | utr_changed | 1 | 1.280 | 1.16–1.41 | 1.04e-06 | 0.00 | 0.415 | 9.99e-04 |
| disease | go_visible | biotype_switch | 1 | 0.892 | 0.81–0.99 | 2.50e-02 | 0.00 | 0.377 | 1.00e+00 |
| disease | go_visible | cds_changed | 1 | 1.059 | 1.01–1.11 | 1.44e-02 | 0.00 | 0.667 | 9.99e-04 |
| disease | go_visible | coding_consequence | 1 | 1.059 | 1.01–1.11 | 1.44e-02 | 0.00 | 0.667 | 9.99e-04 |
| disease | go_visible | coding_status_change | 1 | 0.876 | 0.78–0.99 | 2.72e-02 | 0.00 | 0.271 | 1.00e+00 |
| disease | go_visible | first_exon_changed | 1 | 0.985 | 0.97–1.00 | 9.63e-02 | 0.00 | 0.952 | 1.00e+00 |
| disease | go_visible | internal_exon_difference | 1 | 0.982 | 0.96–1.01 | 1.60e-01 | 0.00 | 0.879 | 9.94e-01 |
| disease | go_visible | last_exon_changed | 1 | 0.976 | 0.95–1.00 | 9.43e-02 | 0.00 | 0.911 | 1.00e+00 |
| disease | go_visible | nmd_switch | 1 | 0.948 | 0.75–1.20 | 6.62e-01 | 0.00 | 0.122 | 7.94e-01 |
| disease | go_visible | utr_changed | 1 | 1.243 | 1.12–1.38 | 6.65e-05 | 0.00 | 0.454 | 9.99e-04 |
