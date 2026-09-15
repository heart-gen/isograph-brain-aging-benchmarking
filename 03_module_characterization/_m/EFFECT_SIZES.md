# Conditional effect sizes for the unique-switch-information test

The de-confounded gene-level test asks whether the switch channel carries phenotype signal the abundance channel does not, and vice versa. Reporting only the FDR-significant count makes the answer a function of sample size. These are the effect-size distributions behind those counts: `partial_r2_switch_given_abund` is the proportion of residual variance the switch block explains after the covariates and abundance are already in the model (`partial_r2_abund_given_switch` is the mirror image).

Read the FDR-significant categories against `neither`, which is the same analysis's own noise floor. `composition_unique` = switch-significant only; `abundance_unique` = abundance-significant only; `both` = both.

## `partial_r2_abund_given_switch`

| analysis | region | adj | n | median (all) | q90 (all) | median: composition_unique | median: neither | ratio vs neither |
|---|---|---|---|---|---|---|---|---|
| brainseq-aging | caudate | no | 13,177 | 0.0166 | 0.0591 | 0.0121 | 0.0113 | 1.1x |
| brainseq-aging | dlpfc | no | 12,931 | 0.0144 | 0.0651 | 0.0119 | 0.0104 | 1.1x |
| brainseq-aging | hippocampus | no | 13,130 | 0.0133 | 0.0400 | 0.0023 | 0.0127 | 0.2x |
| brainseq-sczd | caudate | no | 13,222 | 0.0043 | 0.0260 | 0.0030 | 0.0025 | 1.2x |
| gtex-aging | amygdala | no | 12,042 | 0.0178 | 0.0545 | 0.0097 | 0.0147 | 0.7x |
| gtex-aging | anterior_cingulate_cortex_ba24 | no | 12,190 | 0.0166 | 0.0573 | 0.0094 | 0.0104 | 0.9x |
| gtex-aging | caudate_basal_ganglia | no | 12,424 | 0.0067 | 0.0241 | 0.0083 | 0.0063 | 1.3x |
| gtex-aging | cerebellar_hemisphere | no | 12,416 | 0.0074 | 0.0245 | 0.0158 | 0.0071 | 2.2x |
| gtex-aging | cerebellum | no | 12,480 | 0.0135 | 0.0444 | 0.0095 | 0.0094 | 1.0x |
| gtex-aging | cortex | no | 12,426 | 0.0100 | 0.0398 | 0.0072 | 0.0079 | 0.9x |
| gtex-aging | frontal_cortex_ba9 | no | 12,322 | 0.0169 | 0.0460 | 0.0111 | 0.0106 | 1.1x |
| gtex-aging | hippocampus | no | 12,313 | 0.0115 | 0.0422 | 0.0067 | 0.0090 | 0.7x |
| gtex-aging | hypothalamus | no | 12,648 | 0.0067 | 0.0267 | 0.0093 | 0.0063 | 1.5x |
| gtex-aging | nucleus_accumbens_basal_ganglia | no | 12,492 | 0.0063 | 0.0227 | 0.0165 | 0.0062 | 2.7x |
| gtex-aging | putamen_basal_ganglia | no | 12,100 | 0.0130 | 0.0366 | 0.0172 | 0.0116 | 1.5x |
| gtex-aging | spinal_cord_cervical_c_1 | no | 12,305 | 0.0081 | 0.0276 | 0.0022 | 0.0081 | 0.3x |
| gtex-aging | substantia_nigra | no | 12,265 | 0.0116 | 0.0352 | 0.0082 | 0.0116 | 0.7x |

## `partial_r2_switch_given_abund`

| analysis | region | adj | n | median (all) | q90 (all) | median: composition_unique | median: neither | ratio vs neither |
|---|---|---|---|---|---|---|---|---|
| brainseq-aging | caudate | no | 13,177 | 0.0080 | 0.0277 | 0.0815 | 0.0079 | 10.3x |
| brainseq-aging | dlpfc | no | 12,931 | 0.0077 | 0.0267 | 0.1071 | 0.0074 | 14.4x |
| brainseq-aging | hippocampus | no | 13,130 | 0.0079 | 0.0267 | 0.1209 | 0.0079 | 15.3x |
| brainseq-sczd | caudate | no | 13,222 | 0.0019 | 0.0117 | 0.0377 | 0.0018 | 20.4x |
| gtex-aging | amygdala | no | 12,042 | 0.0105 | 0.0349 | 0.0990 | 0.0107 | 9.3x |
| gtex-aging | anterior_cingulate_cortex_ba24 | no | 12,190 | 0.0110 | 0.0393 | 0.0547 | 0.0101 | 5.4x |
| gtex-aging | caudate_basal_ganglia | no | 12,424 | 0.0060 | 0.0213 | 0.0563 | 0.0060 | 9.4x |
| gtex-aging | cerebellar_hemisphere | no | 12,416 | 0.0061 | 0.0200 | 0.0782 | 0.0061 | 12.9x |
| gtex-aging | cerebellum | no | 12,480 | 0.0071 | 0.0238 | 0.0610 | 0.0072 | 8.5x |
| gtex-aging | cortex | no | 12,426 | 0.0089 | 0.0357 | 0.0480 | 0.0077 | 6.3x |
| gtex-aging | frontal_cortex_ba9 | no | 12,322 | 0.0098 | 0.0338 | 0.0461 | 0.0091 | 5.1x |
| gtex-aging | hippocampus | no | 12,313 | 0.0075 | 0.0247 | 0.0723 | 0.0076 | 9.5x |
| gtex-aging | hypothalamus | no | 12,648 | 0.0084 | 0.0291 | 0.0536 | 0.0080 | 6.7x |
| gtex-aging | nucleus_accumbens_basal_ganglia | no | 12,492 | 0.0054 | 0.0182 | 0.0832 | 0.0054 | 15.5x |
| gtex-aging | putamen_basal_ganglia | no | 12,100 | 0.0063 | 0.0210 | 0.1116 | 0.0061 | 18.1x |
| gtex-aging | spinal_cord_cervical_c_1 | no | 12,305 | 0.0083 | 0.0264 | 0.1198 | 0.0083 | 14.4x |
| gtex-aging | substantia_nigra | no | 12,265 | 0.0108 | 0.0341 | 0.1088 | 0.0108 | 10.1x |

