# Is estimated cell-type composition itself age-coupled?

Spearman correlation of each MuSiC cell-type proportion with donor age, per analysis, read from the committed fractions and the bundle sample tables. Nothing is re-deconvolved and no model is refit.

29 of 94 analysis x cell-type tests reach FDR < 0.05 (Benjamini-Hochberg within analysis) across 12 deconvolved analyses.

## Strongest age-coupled cell type vs switch-unique persistence

| label | class | top_cell_type | top_rho_age | median_abs_rho_age | n_celltypes_age_fdr05 | n_celltypes_tested | comp_unique_base | n_overlap | persistence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SCZD (caudate) | Disease (SCZD) | Immune | 0.367 | 0.0931 | 3 | 8 | 60 | 5 | 0.0833 |
| aging DLPFC | Cortical | Micro | -0.341 | 0.192 | 4 | 8 | 22 | 8 | 0.364 |
| aging caudate | Limbic / striatal | Immune | 0.303 | 0.0969 | 1 | 8 | 45 | 10 | 0.222 |
| GTEx frontal_cortex_ba9 | Cortical | Inhib | -0.299 | 0.17 | 5 | 8 | 797 | 1 | 0.00125 |
| GTEx anterior_cingulate_cortex_ba24 | Cortical | Micro | 0.296 | 0.163 | 3 | 6 | 928 | 293 | 0.316 |
| GTEx cortex | Cortical | Astro | 0.26 | 0.13 | 4 | 8 | 1169 | 0 | 0 |
| GTEx hippocampus | Limbic / striatal | Micro | 0.254 | 0.0965 | 3 | 8 | 44 | 21 | 0.477 |
| GTEx amygdala | Limbic / striatal | Micro | 0.253 | 0.118 | 2 | 8 | 33 | 4 | 0.121 |
| GTEx putamen_basal_ganglia | Limbic / striatal | Immune | 0.195 | 0.0792 | 1 | 8 | 2 | 2 | 1 |
| GTEx caudate_basal_ganglia | Limbic / striatal | Immune | 0.184 | 0.105 | 2 | 8 | 87 | 19 | 0.218 |
| GTEx nucleus_accumbens_basal_ganglia | Limbic / striatal | Inhib | -0.183 | 0.0743 | 1 | 8 | 2 | 1 | 0.5 |
| aging hippocampus | Limbic / striatal | Micro | -0.128 | 0.0716 | 0 | 8 | 1 | 1 | 1 |

## Interpretation

- `persistence` is `n_overlap / comp_unique_base`: the share of a region's unadjusted switch-unique genes that are still switch-unique after adjustment. It is a real fraction; the ratio of the two counts is not, because adjustment adds genes as well as dropping them.
- **A strong `top_rho_age` marks a region where adjustment and the age term compete for the same variance.** In those regions a collapse in switch-unique genes is equally consistent with composition confounding the signal and with over-adjustment removing real signal; this table bounds the risk rather than resolving it.
- **Age-coupling alone does not predict the collapse.** Cortical regions do carry the most age-coupled reference panels, but ACC BA24 and BrainSEQ DLPFC are cortical and retain a third of their switch-unique genes, while GTEx BA9 and cortex retain almost none. Composition-age coupling is a necessary part of the over-adjustment story, not a sufficient one; the reference match matters too (GTEx cortex/BA9 are deconvolved against a DLPFC snRNA panel).
- Correlation is not causation in either direction. Read this panel next to the marker-depletion cut in `COMPOSITION_ADJUSTMENT_SUMMARY.md`, which does not depend on the covariate model at all.