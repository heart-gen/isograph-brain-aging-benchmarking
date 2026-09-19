# Switch coding-consequence enrichment — cross-region rollup

Per consequence class and GO-stratum, pooled separately for the aging and disease analyses. The primary estimate is the **random-effects pooled enrichment ratio** (DerSimonian-Laird) over the per-analysis log enrichments, with a 95% CI and I² for between-analysis heterogeneity. The Fisher-combined p is secondary and should not be read as an effect size: it grows more significant with every added analysis regardless of magnitude. `coding_consequence` = CDS change or coding-status change.

| class | stratum | consequence | k | pooled ratio (RE) | 95% CI | p (RE) | I² | median obs rate | Fisher p |
|-------|---------|-------------|---|-------------------|--------|--------|----|-----------------|----------|
| aging | all | biotype_switch | 9 | 0.896 | 0.88–0.91 | 2.79e-34 | 0.83 | 0.483 | 1.00e+00 |
| aging | all | cds_changed | 9 | 1.040 | 1.04–1.04 | 1.45e-126 | 0.00 | 0.785 | 3.58e-16 |
| aging | all | coding_consequence | 9 | 1.040 | 1.04–1.04 | 1.45e-126 | 0.00 | 0.785 | 3.58e-16 |
| aging | all | coding_status_change | 9 | 0.884 | 0.86–0.91 | 4.93e-15 | 0.90 | 0.346 | 1.00e+00 |
| aging | all | first_exon_changed | 9 | 0.994 | 0.99–1.00 | 2.40e-25 | 0.00 | 0.966 | 1.00e+00 |
| aging | all | internal_exon_difference | 9 | 0.984 | 0.98–0.99 | 4.78e-52 | 0.27 | 0.913 | 1.00e+00 |
| aging | all | last_exon_changed | 9 | 0.958 | 0.95–0.96 | 2.79e-47 | 0.86 | 0.917 | 1.00e+00 |
| aging | all | nmd_switch | 9 | 0.933 | 0.92–0.95 | 5.25e-20 | 0.00 | 0.165 | 1.00e+00 |
| aging | all | utr_changed | 9 | 1.247 | 1.22–1.28 | 5.05e-81 | 0.83 | 0.498 | 6.29e-18 |
| aging | go_invisible | biotype_switch | 6 | 0.869 | 0.83–0.91 | 3.08e-08 | 0.91 | 0.443 | 1.00e+00 |
| aging | go_invisible | cds_changed | 6 | 1.033 | 1.03–1.04 | 2.30e-17 | 0.13 | 0.759 | 5.44e-10 |
| aging | go_invisible | coding_consequence | 6 | 1.033 | 1.03–1.04 | 2.30e-17 | 0.13 | 0.759 | 5.44e-10 |
| aging | go_invisible | coding_status_change | 6 | 0.842 | 0.77–0.92 | 1.55e-04 | 0.94 | 0.296 | 1.00e+00 |
| aging | go_invisible | first_exon_changed | 6 | 0.995 | 0.99–1.00 | 1.83e-05 | 0.00 | 0.969 | 1.00e+00 |
| aging | go_invisible | internal_exon_difference | 6 | 0.979 | 0.97–0.98 | 4.72e-16 | 0.42 | 0.908 | 1.00e+00 |
| aging | go_invisible | last_exon_changed | 6 | 0.952 | 0.94–0.97 | 4.21e-11 | 0.91 | 0.900 | 1.00e+00 |
| aging | go_invisible | nmd_switch | 6 | 0.915 | 0.89–0.94 | 1.36e-08 | 0.00 | 0.166 | 1.00e+00 |
| aging | go_invisible | utr_changed | 6 | 1.266 | 1.19–1.35 | 4.05e-13 | 0.91 | 0.501 | 1.15e-12 |
| aging | go_visible | biotype_switch | 8 | 0.905 | 0.90–0.91 | 4.09e-83 | 0.38 | 0.493 | 1.00e+00 |
| aging | go_visible | cds_changed | 8 | 1.043 | 1.04–1.05 | 1.71e-108 | 0.00 | 0.795 | 8.35e-13 |
| aging | go_visible | coding_consequence | 8 | 1.043 | 1.04–1.05 | 1.71e-108 | 0.00 | 0.795 | 8.35e-13 |
| aging | go_visible | coding_status_change | 8 | 0.898 | 0.88–0.91 | 8.69e-33 | 0.60 | 0.355 | 1.00e+00 |
| aging | go_visible | first_exon_changed | 8 | 0.994 | 0.99–1.00 | 1.71e-21 | 0.00 | 0.966 | 1.00e+00 |
| aging | go_visible | internal_exon_difference | 8 | 0.986 | 0.98–0.99 | 8.61e-53 | 0.00 | 0.915 | 1.00e+00 |
| aging | go_visible | last_exon_changed | 8 | 0.962 | 0.96–0.97 | 1.53e-71 | 0.67 | 0.920 | 1.00e+00 |
| aging | go_visible | nmd_switch | 8 | 0.940 | 0.92–0.96 | 1.23e-12 | 0.00 | 0.162 | 1.00e+00 |
| aging | go_visible | utr_changed | 8 | 1.242 | 1.23–1.26 | 1.26e-247 | 0.31 | 0.489 | 3.54e-16 |
| disease | all | biotype_switch | 1 | 0.941 | 0.92–0.96 | 7.04e-12 | 0.00 | 0.498 | 1.00e+00 |
| disease | all | cds_changed | 1 | 1.040 | 1.03–1.05 | 1.21e-12 | 0.00 | 0.715 | 9.99e-04 |
| disease | all | coding_consequence | 1 | 1.040 | 1.03–1.05 | 1.21e-12 | 0.00 | 0.715 | 9.99e-04 |
| disease | all | coding_status_change | 1 | 0.958 | 0.94–0.98 | 2.50e-04 | 0.00 | 0.368 | 1.00e+00 |
| disease | all | first_exon_changed | 1 | 0.995 | 0.99–1.00 | 9.93e-04 | 0.00 | 0.969 | 1.00e+00 |
| disease | all | internal_exon_difference | 1 | 0.988 | 0.98–0.99 | 3.67e-06 | 0.00 | 0.917 | 1.00e+00 |
| disease | all | last_exon_changed | 1 | 0.982 | 0.98–0.99 | 9.29e-15 | 0.00 | 0.941 | 1.00e+00 |
| disease | all | nmd_switch | 1 | 0.895 | 0.85–0.94 | 7.82e-06 | 0.00 | 0.157 | 1.00e+00 |
| disease | all | utr_changed | 1 | 1.186 | 1.15–1.22 | 2.63e-31 | 0.00 | 0.390 | 9.99e-04 |
| disease | go_invisible | biotype_switch | 1 | 0.938 | 0.92–0.96 | 4.08e-11 | 0.00 | 0.494 | 1.00e+00 |
| disease | go_invisible | cds_changed | 1 | 1.045 | 1.03–1.06 | 5.98e-13 | 0.00 | 0.713 | 9.99e-04 |
| disease | go_invisible | coding_consequence | 1 | 1.045 | 1.03–1.06 | 5.98e-13 | 0.00 | 0.713 | 9.99e-04 |
| disease | go_invisible | coding_status_change | 1 | 0.949 | 0.92–0.97 | 8.85e-05 | 0.00 | 0.360 | 1.00e+00 |
| disease | go_invisible | first_exon_changed | 1 | 0.994 | 0.99–1.00 | 1.10e-03 | 0.00 | 0.969 | 1.00e+00 |
| disease | go_invisible | internal_exon_difference | 1 | 0.989 | 0.98–0.99 | 1.95e-04 | 0.00 | 0.921 | 1.00e+00 |
| disease | go_invisible | last_exon_changed | 1 | 0.981 | 0.98–0.99 | 1.86e-12 | 0.00 | 0.943 | 1.00e+00 |
| disease | go_invisible | nmd_switch | 1 | 0.914 | 0.87–0.96 | 9.95e-04 | 0.00 | 0.162 | 1.00e+00 |
| disease | go_invisible | utr_changed | 1 | 1.210 | 1.17–1.25 | 2.77e-32 | 0.00 | 0.394 | 9.99e-04 |
| disease | go_visible | biotype_switch | 1 | 0.953 | 0.92–0.99 | 1.15e-02 | 0.00 | 0.517 | 1.00e+00 |
| disease | go_visible | cds_changed | 1 | 1.016 | 0.99–1.04 | 1.78e-01 | 0.00 | 0.722 | 2.40e-02 |
| disease | go_visible | coding_consequence | 1 | 1.016 | 0.99–1.04 | 1.78e-01 | 0.00 | 0.722 | 2.40e-02 |
| disease | go_visible | coding_status_change | 1 | 0.993 | 0.95–1.04 | 7.81e-01 | 0.00 | 0.402 | 6.65e-01 |
| disease | go_visible | first_exon_changed | 1 | 0.998 | 0.99–1.01 | 5.57e-01 | 0.00 | 0.967 | 8.21e-01 |
| disease | go_visible | internal_exon_difference | 1 | 0.983 | 0.97–0.99 | 4.91e-03 | 0.00 | 0.902 | 1.00e+00 |
| disease | go_visible | last_exon_changed | 1 | 0.985 | 0.98–1.00 | 4.52e-03 | 0.00 | 0.935 | 1.00e+00 |
| disease | go_visible | nmd_switch | 1 | 0.812 | 0.72–0.91 | 5.37e-04 | 0.00 | 0.137 | 1.00e+00 |
| disease | go_visible | utr_changed | 1 | 1.086 | 1.02–1.16 | 1.00e-02 | 0.00 | 0.374 | 9.99e-04 |
