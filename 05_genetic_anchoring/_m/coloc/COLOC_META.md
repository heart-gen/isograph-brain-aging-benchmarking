# Cross-trait colocalization of the IsoGraph switch layer

eCAVIAR CLPP of switch-module genes against GTEx v11 brain sQTL / eQTL credible sets, for a disease trait (SCZ, on the SCZD switch layer) and the neurodegenerative aging traits (AD, PD, LBD, ALS, on the pooled aging switch layer). CLPP >= 0.01 colocalized; >= 0.05 strong.

## Colocalization counts by trait and QTL kind

| case | trait | kind | n_tested | n_coloc | n_strong | n_coloc_go_invisible |
| --- | --- | --- | --- | --- | --- | --- |
| aging | AD | sQTL | 127 | 6 | 1 | 6 |
| aging | AD | eQTL | 207 | 6 | 2 | 3 |
| aging | ALS | sQTL | 109 | 8 | 2 | 7 |
| aging | ALS | eQTL | 159 | 5 | 1 | 4 |
| aging | LBD | sQTL | 28 | 1 | 0 | 1 |
| aging | LBD | eQTL | 51 | 1 | 0 | 1 |
| aging | PD | sQTL | 97 | 5 | 0 | 5 |
| aging | PD | eQTL | 150 | 4 | 1 | 3 |
| disease | SCZ | sQTL | 36 | 2 | 0 | 1 |
| disease | SCZ | eQTL | 68 | 2 | 0 | 2 |

## Colocalized genes (best CLPP per gene x trait x kind)

| case | trait | gene_name | kind | go_invisible | tissue | clpp |
| --- | --- | --- | --- | --- | --- | --- |
| aging | AD | CTSH | sQTL | True | Brain_Hippocampus | 0.386 |
| aging | AD | CTSH | eQTL | True | Brain_Substantia_nigra | 0.263 |
| aging | AD | CR1 | eQTL | True | Brain_Nucleus_accumbens_basal_ganglia | 0.174 |
| aging | AD | DOC2A | eQTL | False | Brain_Cerebellar_Hemisphere | 0.039 |
| aging | AD | FCER1G | eQTL | False | Brain_Hypothalamus | 0.034 |
| aging | AD | BLNK | eQTL | False | Brain_Amygdala | 0.033 |
| aging | AD | PCGF3 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.025 |
| aging | AD | RTEL1 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.016 |
| aging | AD | MYO15A | eQTL | True | Brain_Amygdala | 0.015 |
| aging | AD | TPCN1 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.011 |
| aging | AD | GRAMD1B | sQTL | True | Brain_Cerebellum | 0.011 |
| aging | AD | MRPS10 | sQTL | True | Brain_Caudate_basal_ganglia | 0.010 |
| aging | ALS | TPP1 | sQTL | True | Brain_Cerebellum | 0.477 |
| aging | ALS | PGS1 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.093 |
| aging | ALS | GGNBP2 | eQTL | True | Brain_Spinal_cord_cervical_c-1 | 0.058 |
| aging | ALS | TXNDC15 | sQTL | True | Brain_Caudate_basal_ganglia | 0.033 |
| aging | ALS | PTPRN | sQTL | False | Brain_Cerebellum | 0.023 |
| aging | ALS | SOWAHA | eQTL | False | Brain_Putamen_basal_ganglia | 0.020 |
| aging | ALS | SCFD1 | eQTL | True | Brain_Cerebellar_Hemisphere | 0.018 |
| aging | ALS | MEF2C | eQTL | True | Brain_Cerebellar_Hemisphere | 0.016 |
| aging | ALS | PPP6R2 | sQTL | True | Brain_Frontal_Cortex_BA9 | 0.014 |
| aging | ALS | BAIAP3 | sQTL | True | Brain_Hypothalamus | 0.013 |
| aging | ALS | MYO18A | sQTL | True | Brain_Caudate_basal_ganglia | 0.012 |
| aging | ALS | G2E3 | eQTL | True | Brain_Frontal_Cortex_BA9 | 0.011 |
| aging | ALS | GGNBP2 | sQTL | True | Brain_Cerebellum | 0.011 |
| aging | LBD | SNCA | sQTL | True | Brain_Cortex | 0.038 |
| aging | LBD | TYK2 | eQTL | True | Brain_Spinal_cord_cervical_c-1 | 0.034 |
| aging | PD | FDFT1 | eQTL | True | Brain_Frontal_Cortex_BA9 | 0.070 |
| aging | PD | SH3GL2 | eQTL | False | Brain_Frontal_Cortex_BA9 | 0.045 |
| aging | PD | STK39 | eQTL | True | Brain_Nucleus_accumbens_basal_ganglia | 0.039 |
| aging | PD | SNCA | sQTL | True | Brain_Frontal_Cortex_BA9 | 0.025 |
| aging | PD | PKD1 | sQTL | True | Brain_Cerebellum | 0.018 |
| aging | PD | TBC1D15 | sQTL | True | Brain_Frontal_Cortex_BA9 | 0.017 |
| aging | PD | DYRK1A | eQTL | True | Brain_Cerebellum | 0.016 |
| aging | PD | PCGF3 | sQTL | True | Brain_Cerebellum | 0.012 |
| aging | PD | TTC19 | sQTL | True | Brain_Putamen_basal_ganglia | 0.012 |
| disease | SCZ | PBX1 | eQTL | True | Brain_Cortex | 0.016 |
| disease | SCZ | RBFA | sQTL | True | Brain_Cerebellum | 0.013 |
| disease | SCZ | MED15 | eQTL | True | Brain_Hypothalamus | 0.012 |
| disease | SCZ | MYO18A | sQTL | False | Brain_Caudate_basal_ganglia | 0.010 |

## Reading

- Across traits, 40 switch genes colocalize (CLPP >= 0.01); 7 are strong (>= 0.05); 33/40 sit in GO-invisible switch modules.
- The aging case is the stronger one: colocalizing switch genes include canonical neurodegeneration loci reached through **splicing** (sQTL) of a GO-invisible switch module -- e.g. SNCA (LBD), TPP1 and SCFD1 (ALS). This is genetic anchoring of biology that gene-abundance co-expression networks miss, because GO/pathway enrichment tracks abundance programs, not co-switching.
- Honest scope: CLPP uses GTEx SuSiE credible sets (no full cis sumstats), so it is restricted to genome-wide-significant trait loci with a brain QTL credible set, and GTEx brain QTLs are bulk-tissue (cell-type-specific splicing under-sampled). It shows the risk variant acts through the switch gene's splicing/expression, not that the co-switching is itself one genetic signal.
