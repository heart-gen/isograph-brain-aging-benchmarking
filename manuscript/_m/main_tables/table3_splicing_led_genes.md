# Table 3 -- Splicing-led colocalized genes (per-gene resolution)

**Table 3. Genes where a disease-GWAS-colocalizing sQTL resolves onto an IsoGraph switch pair (splicing-led genes; eCAVIAR CLPP layer).** This table is the eCAVIAR colocalization layer, retained as orthogonal sensitivity evidence; locus nominations rest on the signal-level hierarchy (coloc.susie > coloc.abf), reported separately, whose gene list overlaps this one only in part. The 12 genes at which a brain sQTL colocalizing with a GWAS credible set (eCAVIAR CLPP) maps onto an IsoGraph switch pair -- the DTU-without-DGE class; all 12 lie in GO-invisible modules. Max CLPP confidence tiers: * suggestive (>0.01, the coloc inclusion bar), ** moderate (>0.05), *** high-confidence (>0.10). **Colocalization posteriors are individually modest** (2 high-confidence, 4 moderate; CTSH is the sole *** case at 0.39): the defensible claim is the *set-level* coherence -- splicing-led equivalent to GO-invisible, cross-disease concordance (SNCA in LBD+PD) -- not any single locus. The set-level statistical support is Table 2 (contrast), not per-gene CLPP; the S-LDSC splicing arm does **not** support it -- the sQTL coefficient clears nominal p<0.05 in 1 of 6 trait-contexts (AD, p=0.0355) with no multiple-testing correction, while eQTL clears 4 of 6. Risk allele aligned to the GWAS trait; LOEUF is gnomAD constraint; Status = established disease isoform biology vs novel candidate. Verbatim from 08_integration/_m/deep_dive/ (deep_dive_panel.parquet, deep_dive_literature.parquet); see Tables S8-S12 for per-event/RBP/clinical/literature layers.

| Gene | Ensembl | Trait(s) | Concordant traits | Lead variant | Risk allele | Top tissue | Max CLPP | n resolved events | GO-invisible | LOEUF | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLG1 | ENSG00000075711 | SCZ | SCZ | rs9843908 | T | Cerebellar Hemisphere | 0.026 * | 2 | no | 0.441 | known isoform biology |
| TMEM175 | ENSG00000127419 | AD,ALS,PD | AD,ALS,PD | rs873786 | C | Cerebellar Hemisphere | 0.483 *** | 6 | yes | 1.295 | nan |
| TPP1 | ENSG00000166340 | ALS | ALS | rs2072651 | G | Cerebellum | 0.480 *** | 1 | yes | 0.746 | nan |
| CRELD2 | ENSG00000184164 | SCZ | SCZ | rs7284417 | G | Cerebellum | 0.062 ** | 5 | yes | 1.09 | nan |
| PPP6R2 | ENSG00000100239 | SCZ | SCZ | rs76300267 | G | Frontal Cortex BA9 | 0.058 ** | 1 | no | 0.766 | novel candidate |
| PRDM2 | ENSG00000116731 | ALS | ALS | rs2744682 | A | Cortex | 0.056 ** | 1 | no | 0.185 | nan |
| CTC1 | ENSG00000178971 | ALS | ALS | rs3027235 | G | Cerebellum | 0.053 ** | 2 | no | 0.776 | nan |
| DOC2A | ENSG00000149927 | AD,SCZ | SCZ | rs11150580 | C | Cerebellar Hemisphere | 0.048 * | 1 | yes | 0.732 | nan |
| SNAP91 | ENSG00000065609 | SCZ | SCZ | rs217342 | T | Cerebellum | 0.047 * | 1 | yes | 0.362 | nan |
| SPG7 | ENSG00000197912 | SCZ | SCZ | rs12919215 | G | Cerebellum | 0.037 * | 10 | yes | 1.433 | nan |
| INO80E | ENSG00000169592 | SCZ | SCZ | rs4788204 | G | Cerebellar Hemisphere | 0.034 * | 3 | yes | 1.282 | nan |
| IFNAR2 | ENSG00000159110 | AD | AD | rs9975538 | C | Cerebellum | 0.029 * | 2 | no | 0.761 | nan |
| DNAJA3 | ENSG00000103423 | SCZ | SCZ | rs6500605 | G | Cerebellar Hemisphere | 0.025 * | 2 | yes | 0.733 | nan |
| THAP3 | ENSG00000041988 | SCZ | SCZ | rs3789572 | T | Cerebellar Hemisphere | 0.023 * | 1 | no | 1.037 | nan |
| PCGF3 | ENSG00000185619 | AD,PD | PD | rs13116048 | A | Cerebellar Hemisphere | 0.021 * | 2 | no | 0.445 | nan |
| PTPRN | ENSG00000054356 | ALS | ALS | rs6436132 | G | Cerebellum | 0.021 * | 1 | yes | 0.735 | nan |
| GSTO2 | ENSG00000065621 | SCZ | SCZ | rs10883990 | A | Frontal Cortex BA9 | 0.021 * | 4 | yes | 1.238 | nan |
| RPAIN | ENSG00000129197 | SCZ | SCZ | rs2189336 | G | Cortex | 0.020 * | 2 | no | 1.32 | nan |
| TBC1D15 | ENSG00000121749 | PD | PD | rs61754230 | T | Frontal Cortex BA9 | 0.017 * | 1 | no | 0.502 | novel candidate |
| NT5C2 | ENSG00000076685 | SCZ | SCZ | rs12412038 | G | Cerebellum | 0.017 * | 1 | no | 1.128 | nan |
| NADSYN1 | ENSG00000172890 | SCZ | SCZ | rs12803256 | G | Cerebellar Hemisphere | 0.016 * | 8 | yes | 0.929 | nan |
| B3GAT1 | ENSG00000109956 | SCZ | SCZ | rs1440480 | A | Caudate basal ganglia | 0.015 * | 2 | yes | 0.749 | nan |
| PIGQ | ENSG00000007541 | ALS | ALS | rs2891650 | G | Putamen basal ganglia | 0.015 * | 3 | yes | 1.275 | nan |
| FLCN | ENSG00000154803 | AD,SCZ | AD,SCZ | rs1708618 | C | Cerebellar Hemisphere | 0.014 * | 6 | yes | 0.487 | nan |
| VAMP2 | ENSG00000220205 | ALS | ALS | rs8066511 | C | Cortex | 0.013 * | 1 | no | 0.176 | nan |
| RBFA | ENSG00000101546 | SCZ | SCZ | rs8098075 | C | Cerebellum | 0.013 * | 1 | yes | 0.977 | nan |
| BAIAP3 | ENSG00000007516 | ALS | ALS | rs185060782 | G | Hypothalamus | 0.013 * | 1 | no | 1.383 | nan |
| CTSB | ENSG00000164733 | PD | PD | rs1293295 | C | Amygdala | 0.011 * | 1 | no | 1.705 | nan |
| GPR135 | ENSG00000181619 | SCZ | SCZ | rs77789667 | T | Cerebellum | 0.010 * | 1 | yes | 1.944 | nan |
| TARBP1 | ENSG00000059588 | SCZ | SCZ | rs1039996 | T | Cerebellum | 0.010 * | 3 | yes | 0.95 | nan |
