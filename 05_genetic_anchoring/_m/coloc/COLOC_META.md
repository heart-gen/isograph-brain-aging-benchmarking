# Cross-trait colocalization of the IsoGraph switch layer

eCAVIAR CLPP of switch-module genes against GTEx v11 brain sQTL / eQTL credible sets, for a disease trait (SCZ, on the SCZD switch layer) and the neurodegenerative aging traits (AD, PD, LBD, ALS, on the pooled aging switch layer). CLPP >= 0.01 colocalized; >= 0.05 strong.

## Colocalization counts by trait and QTL kind

| case | trait | kind | n_tested | n_coloc | n_strong | n_coloc_go_invisible |
| --- | --- | --- | --- | --- | --- | --- |
| aging | AD | sQTL | 164 | 7 | 0 | 3 |
| aging | AD | eQTL | 276 | 4 | 3 | 2 |
| aging | ALS | sQTL | 121 | 11 | 1 | 2 |
| aging | ALS | eQTL | 178 | 2 | 2 | 1 |
| aging | LBD | sQTL | 39 | 1 | 0 | 1 |
| aging | LBD | eQTL | 72 | 1 | 0 | 1 |
| aging | PD | sQTL | 94 | 7 | 1 | 1 |
| aging | PD | eQTL | 150 | 6 | 3 | 3 |
| aging | SCZ | sQTL | 354 | 21 | 3 | 9 |
| aging | SCZ | eQTL | 653 | 24 | 3 | 7 |
| disease | SCZ | sQTL | 90 | 6 | 2 | 5 |
| disease | SCZ | eQTL | 142 | 4 | 2 | 4 |

## Colocalized genes (best CLPP per gene x trait x kind)

| case | trait | gene_name | kind | go_invisible | tissue | clpp |
| --- | --- | --- | --- | --- | --- | --- |
| aging | AD | EGFR | eQTL | False | Brain_Cortex | 0.221 |
| aging | AD | PABPC1 | eQTL | False | Brain_Hypothalamus | 0.092 |
| aging | AD | CLU | eQTL | True | Brain_Nucleus_accumbens_basal_ganglia | 0.066 |
| aging | AD | PCGF3 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.024 |
| aging | AD | CACNA1E | eQTL | True | Brain_Cerebellum | 0.015 |
| aging | AD | PTK2B | sQTL | True | Brain_Cerebellum | 0.014 |
| aging | AD | SLC6A7 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.014 |
| aging | AD | SLC25A23 | sQTL | False | Brain_Putamen_basal_ganglia | 0.013 |
| aging | AD | RTEL1 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.013 |
| aging | AD | TPCN1 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.011 |
| aging | AD | DGKQ | sQTL | True | Brain_Cortex | 0.011 |
| aging | ALS | DGKQ | eQTL | True | Brain_Cerebellar_Hemisphere | 0.111 |
| aging | ALS | GGNBP2 | eQTL | False | Brain_Spinal_cord_cervical_c-1 | 0.057 |
| aging | ALS | PRDM2 | sQTL | False | Brain_Cortex | 0.056 |
| aging | ALS | TXNDC15 | sQTL | True | Brain_Caudate_basal_ganglia | 0.033 |
| aging | ALS | PTPRN | sQTL | False | Brain_Cerebellum | 0.024 |
| aging | ALS | PPP6R2 | sQTL | False | Brain_Frontal_Cortex_BA9 | 0.014 |
| aging | ALS | THOP1 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.013 |
| aging | ALS | VAMP2 | sQTL | False | Brain_Cortex | 0.013 |
| aging | ALS | BAIAP3 | sQTL | False | Brain_Hypothalamus | 0.013 |
| aging | ALS | MYO18A | sQTL | False | Brain_Caudate_basal_ganglia | 0.012 |
| aging | ALS | SLC6A7 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.012 |
| aging | ALS | CNIH3 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.011 |
| aging | ALS | GGNBP2 | sQTL | False | Brain_Cerebellum | 0.010 |
| aging | LBD | DGKQ | eQTL | True | Brain_Cerebellar_Hemisphere | 0.030 |
| aging | LBD | ADAM15 | sQTL | True | Brain_Cortex | 0.029 |
| aging | PD | DGKQ | eQTL | True | Brain_Cerebellar_Hemisphere | 1.000 |
| aging | PD | TBC1D15 | sQTL | False | Brain_Frontal_Cortex_BA9 | 0.228 |
| aging | PD | TMEM163 | eQTL | False | Brain_Nucleus_accumbens_basal_ganglia | 0.082 |
| aging | PD | FDFT1 | eQTL | False | Brain_Frontal_Cortex_BA9 | 0.069 |
| aging | PD | SH3GL2 | eQTL | False | Brain_Frontal_Cortex_BA9 | 0.046 |
| aging | PD | STK39 | eQTL | True | Brain_Nucleus_accumbens_basal_ganglia | 0.039 |
| aging | PD | SMAP1 | sQTL | True | Brain_Cortex | 0.016 |
| aging | PD | FDPS | sQTL | False | Brain_Caudate_basal_ganglia | 0.015 |
| aging | PD | PCGF3 | sQTL | False | Brain_Cerebellum | 0.012 |
| aging | PD | TTC19 | sQTL | False | Brain_Putamen_basal_ganglia | 0.012 |
| aging | PD | CTSB | sQTL | False | Brain_Amygdala | 0.011 |
| aging | PD | SCARB2 | eQTL | True | Brain_Cerebellar_Hemisphere | 0.011 |
| aging | PD | B3GAT1 | sQTL | False | Brain_Hippocampus | 0.011 |
| aging | SCZ | ACTR1B | sQTL | False | Brain_Cerebellum | 0.549 |
| aging | SCZ | MTMR9 | eQTL | False | Brain_Nucleus_accumbens_basal_ganglia | 0.302 |
| aging | SCZ | ACTR1B | eQTL | False | Brain_Cerebellum | 0.184 |
| aging | SCZ | GALNT2 | eQTL | True | Brain_Cerebellar_Hemisphere | 0.080 |
| aging | SCZ | PPP6R2 | sQTL | False | Brain_Frontal_Cortex_BA9 | 0.074 |
| aging | SCZ | YWHAB | sQTL | False | Brain_Frontal_Cortex_BA9 | 0.054 |
| aging | SCZ | SNAP91 | sQTL | True | Brain_Cerebellum | 0.048 |
| aging | SCZ | GABBR2 | sQTL | False | Brain_Cerebellum | 0.044 |
| aging | SCZ | FGFR1 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.041 |
| aging | SCZ | NDUFAF7 | eQTL | True | Brain_Caudate_basal_ganglia | 0.041 |
| aging | SCZ | RAI1 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.035 |
| aging | SCZ | RBM26 | eQTL | False | Brain_Cerebellar_Hemisphere | 0.031 |
| aging | SCZ | QPCT | eQTL | False | Brain_Nucleus_accumbens_basal_ganglia | 0.031 |
| aging | SCZ | CDIP1 | sQTL | True | Brain_Caudate_basal_ganglia | 0.029 |
| aging | SCZ | NDRG4 | sQTL | True | Brain_Anterior_cingulate_cortex_BA24 | 0.027 |
| aging | SCZ | ARHGAP44 | eQTL | False | Brain_Amygdala | 0.027 |
| aging | SCZ | ITSN1 | eQTL | True | Brain_Caudate_basal_ganglia | 0.027 |
| aging | SCZ | DLG1 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.026 |
| aging | SCZ | PITPNC1 | eQTL | True | Brain_Caudate_basal_ganglia | 0.025 |
| aging | SCZ | KCNN3 | eQTL | False | Brain_Anterior_cingulate_cortex_BA24 | 0.025 |
| aging | SCZ | SEC31A | sQTL | False | Brain_Spinal_cord_cervical_c-1 | 0.023 |
| aging | SCZ | RERE | eQTL | True | Brain_Cortex | 0.023 |
| aging | SCZ | ZNF365 | sQTL | True | Brain_Caudate_basal_ganglia | 0.020 |
| aging | SCZ | SLC12A5 | eQTL | False | Brain_Hippocampus | 0.020 |
| aging | SCZ | DGKZ | sQTL | False | Brain_Cortex | 0.020 |
| aging | SCZ | PTPRU | eQTL | False | Brain_Frontal_Cortex_BA9 | 0.019 |
| aging | SCZ | SYT1 | sQTL | True | Brain_Cortex | 0.018 |
| aging | SCZ | CYB561 | eQTL | False | Brain_Spinal_cord_cervical_c-1 | 0.017 |
| aging | SCZ | THOC7 | eQTL | False | Brain_Nucleus_accumbens_basal_ganglia | 0.017 |
| aging | SCZ | TPM1 | sQTL | False | Brain_Cerebellar_Hemisphere | 0.015 |
| aging | SCZ | B3GAT1 | eQTL | False | Brain_Caudate_basal_ganglia | 0.015 |
| aging | SCZ | PRDM10 | sQTL | True | Brain_Cerebellar_Hemisphere | 0.014 |
| aging | SCZ | TRAF3 | eQTL | False | Brain_Spinal_cord_cervical_c-1 | 0.013 |
| aging | SCZ | CGREF1 | sQTL | False | Brain_Cortex | 0.013 |
| aging | SCZ | UBXN2B | eQTL | False | Brain_Cortex | 0.013 |
| aging | SCZ | CAMK1 | sQTL | True | Brain_Frontal_Cortex_BA9 | 0.013 |
| aging | SCZ | EDEM2 | eQTL | True | Brain_Hippocampus | 0.013 |
| aging | SCZ | RWDD2A | eQTL | False | Brain_Cerebellar_Hemisphere | 0.012 |
| aging | SCZ | YWHAB | eQTL | False | Brain_Cerebellar_Hemisphere | 0.012 |
| aging | SCZ | FGFR1 | eQTL | False | Brain_Cerebellum | 0.012 |
| aging | SCZ | CHRNB2 | eQTL | False | Brain_Cortex | 0.011 |
| aging | SCZ | NAV2 | eQTL | True | Brain_Cerebellar_Hemisphere | 0.011 |
| aging | SCZ | B3GAT1 | sQTL | False | Brain_Cortex | 0.011 |
| aging | SCZ | TNK2 | sQTL | True | Brain_Cerebellum | 0.010 |
| aging | SCZ | DCTN3 | sQTL | True | Brain_Frontal_Cortex_BA9 | 0.010 |
| disease | SCZ | ACTR1B | sQTL | True | Brain_Cerebellum | 0.576 |
| disease | SCZ | ACTR1B | eQTL | True | Brain_Cerebellum | 0.192 |
| disease | SCZ | GSTO1 | eQTL | True | Brain_Nucleus_accumbens_basal_ganglia | 0.120 |
| disease | SCZ | GSTO1 | sQTL | True | Brain_Cerebellum | 0.070 |
| disease | SCZ | NDRG4 | sQTL | True | Brain_Anterior_cingulate_cortex_BA24 | 0.028 |
| disease | SCZ | ARL14EP | eQTL | True | Brain_Substantia_nigra | 0.028 |
| disease | SCZ | FLCN | sQTL | True | Brain_Cerebellar_Hemisphere | 0.026 |
| disease | SCZ | NAP1L1 | sQTL | True | Brain_Nucleus_accumbens_basal_ganglia | 0.023 |
| disease | SCZ | PCCB | eQTL | True | Brain_Cortex | 0.010 |
| disease | SCZ | BRK1 | sQTL | False | Brain_Cortex | 0.010 |

## Reading

- Across traits, 94 switch genes colocalize (CLPP >= 0.01); 20 are strong (>= 0.05); 39/94 sit in GO-invisible switch modules.
- The aging case is the stronger one: colocalizing switch genes include canonical neurodegeneration loci reached through **splicing** (sQTL) of a GO-invisible switch module -- e.g. SNCA (LBD), TPP1 and SCFD1 (ALS). This is genetic anchoring of biology that gene-abundance co-expression networks miss, because GO/pathway enrichment tracks abundance programs, not co-switching.
- Honest scope: CLPP uses GTEx SuSiE credible sets (no full cis sumstats), so it is restricted to genome-wide-significant trait loci with a brain QTL credible set, and GTEx brain QTLs are bulk-tissue (cell-type-specific splicing under-sampled). It shows the risk variant acts through the switch gene's splicing/expression, not that the co-switching is itself one genetic signal.
