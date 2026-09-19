# Table 3 -- Splicing-led colocalized genes (per-gene resolution)

**Table 3. Genes where a disease-GWAS-colocalizing sQTL resolves onto an IsoGraph switch pair (splicing-led genes; eCAVIAR CLPP layer).** This table is the eCAVIAR colocalization layer, retained as orthogonal sensitivity evidence; locus nominations rest on the signal-level hierarchy (coloc.susie > coloc.abf), reported separately, whose gene list overlaps this one only in part. The 12 genes at which a brain sQTL colocalizing with a GWAS credible set (eCAVIAR CLPP) maps onto an IsoGraph switch pair -- the DTU-without-DGE class; all 12 lie in GO-invisible modules. Max CLPP confidence tiers: * suggestive (>0.01, the coloc inclusion bar), ** moderate (>0.05), *** high-confidence (>0.10). **Colocalization posteriors are individually modest** (1 high-confidence, 3 moderate; CTSH is the sole *** case at 0.39): the defensible claim is the *set-level* coherence -- splicing-led equivalent to GO-invisible, cross-disease concordance (SNCA in LBD+PD) -- not any single locus. The set-level statistical support is Table 2 (contrast), not per-gene CLPP; the S-LDSC splicing arm does **not** support it -- the sQTL coefficient clears nominal p<0.05 in 1 of 6 trait-contexts (AD, p=0.0355) with no multiple-testing correction, while eQTL clears 4 of 6. Risk allele aligned to the GWAS trait; LOEUF is gnomAD constraint; Status = established disease isoform biology vs novel candidate. Verbatim from 05_genetic_anchoring/_m/deep_dive/ (deep_dive_panel.parquet, deep_dive_literature.parquet); see Tables S8-S12 for per-event/RBP/clinical/literature layers.

| Gene | Ensembl | Trait(s) | Concordant traits | Lead variant | Risk allele | Top tissue | Max CLPP | n resolved events | GO-invisible | LOEUF | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CTSH | ENSG00000103811 | AD | AD | rs12148472 | T | Hippocampus | 0.386 *** | 1 | yes | 1.178 | known isoform biology |
| SNCA | ENSG00000145335 | LBD,PD | LBD,PD | rs7680557 | A | Cortex | 0.038 * | 2 | yes | 0.397 | known isoform biology |
| DLG1 | ENSG00000075711 | SCZ | SCZ | rs9843908 | T | Cerebellar Hemisphere | 0.026 * | 2 | yes | 0.441 | known isoform biology |
| ARVCF | ENSG00000099889 | SCZ | SCZ | rs4819527 | C | Cerebellar Hemisphere | 0.010 * | 1 | yes | 0.999 | known isoform biology |
| PGS1 | ENSG00000087157 | ALS | ALS | rs11656568 | A | Cerebellar Hemisphere | 0.093 ** | 1 | yes | 0.998 | novel candidate |
| PPP6R2 | ENSG00000100239 | ALS,SCZ | ALS,SCZ | rs76300267 | G | Frontal Cortex BA9 | 0.065 ** | 3 | yes | 0.766 | novel candidate |
| GGNBP2 | ENSG00000278311 | ALS,SCZ | ALS | rs9903355 | C | Spinal cord cervical c-1 | 0.058 ** | 1 | yes | 0.201 | novel candidate |
| CDIP1 | ENSG00000089486 | SCZ | SCZ | rs4786493 | C | Caudate basal ganglia | 0.028 * | 2 | yes | 1.132 | novel candidate |
| TBC1D15 | ENSG00000121749 | PD | PD | rs61754230 | T | Frontal Cortex BA9 | 0.017 * | 1 | yes | 0.502 | novel candidate |
| PRRC2B | ENSG00000288701 | SCZ | SCZ | rs113866326 | C | Anterior cingulate cortex BA24 | 0.017 * | 1 | yes | 0.338 | novel candidate |
| RTEL1 | ENSG00000258366 | AD,SCZ | SCZ | rs61753459 | C | Cerebellar Hemisphere | 0.016 * | 1 | yes | 0.638 | novel candidate |
| TPCN1 | ENSG00000186815 | AD | AD | rs3815889 | G | Cerebellar Hemisphere | 0.011 * | 1 | yes | 0.558 | novel candidate |
