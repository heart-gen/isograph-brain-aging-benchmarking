# Cell-type composition adjustment of the DTU-without-DGE layer

Re-running the de-confounded gene-level switch-vs-abundance test with MuSiC cell-type fractions added as inference covariates (reference cell type dropped for the simplex). `comp_unique` = genes whose isoform composition is phenotype-associated beyond their own abundance (IsoGraph's unique signal).

## BrainSEQ (disease + primary aging)

| region | n_tested | comp_unique_base | comp_unique_adj | retained_frac | both_base | both_adj | n_cell_types | samples_covered | marker_depleted | marker_enriched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SCZD (caudate) | 17193 | 34 | 2 | 0.06 | 16 | 1 | 8 | 390 | 0 | 1 |
| aging caudate | 11648 | 43 | 17 | 0.4 | 34 | 22 | 8 | 238 | 0 | 0 |
| aging hippocampus | 11114 | 0 | 0 |  | 0 | 0 | 8 | 238 | 0 | 0 |
| aging DLPFC | 11554 | 8 | 15 | 1.88 | 4 | 18 | 8 | 222 | 0 | 1 |

## GTEx (aging replication arm)

In-repo MuSiC re-run; 6/8 deconvolved regions retain a composition-robust switch signal and 2/8 collapse to zero — a **region-dependent, partial** replication: limbic/striatal aging DTU reproduces as composition-robust, while the two cortical regions collapse (composition-entangled and/or over-adjusted). Full breakdown + caveats in `real_data/gtex/_m/composition/GTEX_COMPOSITION_SUMMARY.md`.

| region | n_tested | comp_unique_base | comp_unique_adj | retained_frac | both_base | both_adj | n_cell_types | samples_covered | marker_depleted | marker_enriched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GTEx amygdala | 16351 | 7 | 6 | 0.86 | 2 | 0 | 8 | 181 | 0 | 4 |
| GTEx anterior_cingulate_cortex_ba24 | 16250 | 545 | 88 | 0.16 | 107 | 2 | 6 | 233 | 0 | 6 |
| GTEx frontal_cortex_ba9 | 16362 | 531 | 0 | 0.0 | 175 | 0 | 8 | 269 | 0 | 6 |
| GTEx cortex | 16480 | 438 | 0 | 0.0 | 72 | 0 | 8 | 270 | 0 | 1 |
| GTEx hippocampus | 16500 | 61 | 21 | 0.34 | 14 | 3 | 8 | 255 | 0 | 3 |
| GTEx caudate_basal_ganglia | 16543 | 32 | 10 | 0.31 | 7 | 1 | 8 | 300 | 0 | 4 |
| GTEx putamen_basal_ganglia | 16158 | 4 | 3 | 0.75 | 1 | 1 | 8 | 254 | 0 | 5 |
| GTEx nucleus_accumbens_basal_ganglia | 16529 | 1 | 4 | 4.0 | 0 | 0 | 8 | 285 | 0 | 2 |

## Interpretation

- The composition-adjusted `comp_unique_adj` is the honest, composition-robust count that should anchor the DTU-without-DGE claim; report it alongside the unadjusted number rather than in place of it.
- A large base→adj drop (e.g. SCZD) means much of that switch signal co-varies with cell-type proportion. **Confounder vs mediator matters:** if disease/age *causes* the composition shift that drives the switch, covariate adjustment *removes real signal* (over-adjustment); if composition varies for technical/sampling reasons, adjustment is the correct control. State this both-ways and lean on the layers that do not depend on it (module-level genetic anchoring, GO-invisible gate).
- The aging switch layer survives composition adjustment in BrainSEQ (caudate, DLPFC) and **partially replicates in GTEx** (limbic/striatal robust; the two cortical regions collapse — composition-entangled and/or over-adjusted), whereas the SCZD disease signal is largely composition-confounded. Lead the DTU claim with the composition-robust aging layer, and disclose the GTEx cortical collapse as the boundary of that robustness rather than burying it.
- The marker cut (`marker_depleted`/`marker_enriched`) reports whether module member genes are over/under-represented for cell-type marker genes; near-null argues the modules are not simply bags of cell-type markers.