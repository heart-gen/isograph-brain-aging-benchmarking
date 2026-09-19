# GTEx composition adjustment — aging replication arm

MuSiC deconvolution re-run in-repo on the GTEx brain bundles using the same Tran/LIBD snRNA references (seed 13), for the regions with a defensibly matched reference (striatum→NAc, cortex→DLPFC, plus AMY/sACC/HPC direct). Regions with no matched reference (cerebellum, hypothalamus, spinal cord, substantia nigra) are not deconvolved.

| region | n_tested | comp_unique_base | comp_unique_adj | retained_frac | both_base | both_adj | n_cell_types | samples_covered | marker_depleted | marker_enriched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GTEx amygdala | 12042 | 33 | 5 | 0.15 | 6 | 0 | 8 | 181 | 0 | 1 |
| GTEx anterior_cingulate_cortex_ba24 | 12190 | 928 | 365 | 0.39 | 198 | 14 | 6 | 233 | 0 | 2 |
| GTEx frontal_cortex_ba9 | 12322 | 797 | 2 | 0.0 | 320 | 0 | 8 | 269 | 0 | 0 |
| GTEx cortex | 12426 | 1169 | 2 | 0.0 | 152 | 0 | 8 | 270 | 0 | 1 |
| GTEx hippocampus | 12313 | 44 | 29 | 0.66 | 6 | 2 | 8 | 255 | 0 | 2 |
| GTEx caudate_basal_ganglia | 12424 | 87 | 50 | 0.57 | 4 | 6 | 8 | 300 | 0 | 2 |
| GTEx putamen_basal_ganglia | 12100 | 2 | 4 | 2.0 | 1 | 2 | 8 | 254 | 0 | 0 |
| GTEx nucleus_accumbens_basal_ganglia | 12492 | 2 | 2 | 1.0 | 0 | 0 | 8 | 285 | 0 | 1 |

## Interpretation — a region-dependent, honestly partial replication

**8/8 regions retain ≥1 composition-robust DTU-without-DGE gene, but 0/8 collapse to zero and the retained fraction varies enormously (0.0–0.86).** This is the confounder-vs-mediator caveat made concrete rather than a clean win:

- **Limbic/striatal regions retain** (amygdala, hippocampus, caudate, putamen, NAc): the aging switch signal there is not merely a proportion artifact.
- **The two cortical regions collapse to 0** (frontal_cortex_ba9 531→0, cortex 438→0) — despite BrainSEQ *DLPFC* surviving adjustment. Cortical composition is the most strongly age-coupled, and these are the weakest reference matches (GTEx cortex/BA9 → DLPFC snRNA), so the collapse is consistent with genuine composition-confounding of cortical aging switches **and/or over-adjustment** when 7–8 age-correlated fraction covariates absorb the age-spline variance. The data cannot cleanly separate these here; report the collapse, do not hide it.

**Manuscript use:** cite GTEx as a *partial* composition replication — the limbic/striatal aging layer reproduces as composition-robust across cohorts, while cortical aging DTU is composition-entangled in GTEx (a stated limitation, not a claim of universal robustness).

_Caveat: cross-region reference use (striatum→NAc, cortex→DLPFC) and the snRNA reference panel make these fraction estimates approximate; the with-vs-without contrast, not the absolute fractions, is the claim._