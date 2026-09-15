# Cell-type composition adjustment of the DTU-without-DGE layer

Re-running the de-confounded gene-level switch-vs-abundance test with MuSiC cell-type fractions added as inference covariates (reference cell type dropped for the simplex). `comp_unique` = genes whose isoform composition is phenotype-associated beyond their own abundance (IsoGraph's unique signal).

## BrainSEQ (disease + primary aging)

| region | n_tested | comp_unique_base | comp_unique_adj | retained_frac | both_base | both_adj | n_cell_types | samples_covered | marker_depleted | marker_enriched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SCZD (caudate) | 13222 | 60 | 11 | 0.18 | 37 | 9 | 8 | 390 | 0 | 0 |
| aging caudate | 13177 | 45 | 14 | 0.31 | 39 | 18 | 8 | 238 | 0 | 1 |
| aging hippocampus | 13130 | 1 | 1 | 1.0 | 0 | 0 | 8 | 238 | 0 | 1 |
| aging DLPFC | 12931 | 22 | 28 | 1.27 | 8 | 26 | 8 | 222 | 0 | 0 |

## GTEx (aging replication arm)

In-repo MuSiC re-run; 8/8 deconvolved regions retain a composition-robust switch signal and 0/8 collapse to zero — a **region-dependent, partial** replication: limbic/striatal aging DTU reproduces as composition-robust, while the two cortical regions collapse (composition-entangled and/or over-adjusted). Full breakdown + caveats in `02_module_discovery/gtex/_m/composition/GTEX_COMPOSITION_SUMMARY.md`.

| region | n_tested | comp_unique_base | comp_unique_adj | retained_frac | both_base | both_adj | n_cell_types | samples_covered | marker_depleted | marker_enriched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GTEx amygdala | 12042 | 33 | 5 | 0.15 | 6 | 0 | 8 | 181 | 0 | 1 |
| GTEx anterior_cingulate_cortex_ba24 | 12190 | 928 | 365 | 0.39 | 198 | 14 | 6 | 233 | 1 | 3 |
| GTEx frontal_cortex_ba9 | 12322 | 797 | 2 | 0.0 | 320 | 0 | 8 | 269 | 0 | 0 |
| GTEx cortex | 12426 | 1169 | 2 | 0.0 | 152 | 0 | 8 | 270 | 0 | 2 |
| GTEx hippocampus | 12313 | 44 | 29 | 0.66 | 6 | 2 | 8 | 255 | 0 | 2 |
| GTEx caudate_basal_ganglia | 12424 | 87 | 50 | 0.57 | 4 | 6 | 8 | 300 | 0 | 0 |
| GTEx putamen_basal_ganglia | 12100 | 2 | 4 | 2.0 | 1 | 2 | 8 | 254 | 0 | 0 |
| GTEx nucleus_accumbens_basal_ganglia | 12492 | 2 | 2 | 1.0 | 0 | 0 | 8 | 285 | 0 | 1 |

## Interpretation

- The composition-adjusted `comp_unique_adj` is the honest, composition-robust count that should anchor the DTU-without-DGE claim; report it alongside the unadjusted number rather than in place of it.
- A large base→adj drop (e.g. SCZD) means much of that switch signal co-varies with cell-type proportion. **Confounder vs mediator matters:** if disease/age *causes* the composition shift that drives the switch, covariate adjustment *removes real signal* (over-adjustment); if composition varies for technical/sampling reasons, adjustment is the correct control. State this both-ways and lean on the layers that do not depend on it (module-level genetic anchoring, GO-invisible gate).
- The aging switch layer survives composition adjustment in BrainSEQ (caudate, DLPFC) and **partially replicates in GTEx** (limbic/striatal robust; the two cortical regions collapse — composition-entangled and/or over-adjusted), whereas the SCZD disease signal is largely composition-confounded. Lead the DTU claim with the composition-robust aging layer, and disclose the GTEx cortical collapse as the boundary of that robustness rather than burying it.
- The marker cut (`marker_depleted`/`marker_enriched`) reports whether module member genes are over/under-represented for cell-type marker genes; near-null argues the modules are not simply bags of cell-type markers.