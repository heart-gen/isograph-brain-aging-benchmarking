# Conditional effect sizes for the unique-switch-information test

The de-confounded gene-level test asks whether the switch channel carries phenotype signal the abundance channel does not, and vice versa. Reporting only the FDR-significant count makes the answer a function of sample size. These are the effect-size distributions behind those counts: `partial_r2_switch_given_abund` is the proportion of residual variance the switch block explains after the covariates and abundance are already in the model (`partial_r2_abund_given_switch` is the mirror image).

Read the FDR-significant categories against `neither`, which is the same analysis's own noise floor. `composition_unique` = switch-significant only; `abundance_unique` = abundance-significant only; `both` = both.

## `partial_r2_abund_given_switch`

| analysis | region | adj | n | median (all) | q90 (all) | median: composition_unique | median: neither | ratio vs neither |
|---|---|---|---|---|---|---|---|---|
| brainseq-aging | caudate | no | 11,648 | 0.0158 | 0.0577 | 0.0206 | 0.0111 | 1.8x |
| brainseq-aging | dlpfc | no | 11,554 | 0.0140 | 0.0658 | 0.0210 | 0.0102 | 2.1x |
| brainseq-aging | hippocampus | no | 11,114 | 0.0128 | 0.0389 | n/a | 0.0122 | n/a |
| brainseq-sczd | caudate | no | 17,193 | 0.0045 | 0.0267 | 0.0034 | 0.0025 | 1.4x |
| gtex-aging | amygdala | no | 16,351 | 0.0181 | 0.0543 | 0.0117 | 0.0146 | 0.8x |
| gtex-aging | anterior_cingulate_cortex_ba24 | no | 16,250 | 0.0163 | 0.0562 | 0.0094 | 0.0102 | 0.9x |
| gtex-aging | caudate_basal_ganglia | no | 16,543 | 0.0066 | 0.0233 | 0.0091 | 0.0064 | 1.4x |
| gtex-aging | cerebellar_hemisphere | no | 16,140 | 0.0073 | 0.0247 | 0.0225 | 0.0071 | 3.2x |
| gtex-aging | cerebellum | no | 16,382 | 0.0137 | 0.0448 | 0.0011 | 0.0097 | 0.1x |
| gtex-aging | cortex | no | 16,480 | 0.0104 | 0.0409 | 0.0072 | 0.0078 | 0.9x |
| gtex-aging | frontal_cortex_ba9 | no | 16,362 | 0.0167 | 0.0464 | 0.0097 | 0.0104 | 0.9x |
| gtex-aging | hippocampus | no | 16,500 | 0.0112 | 0.0404 | 0.0088 | 0.0090 | 1.0x |
| gtex-aging | hypothalamus | no | 16,898 | 0.0067 | 0.0257 | 0.0082 | 0.0064 | 1.3x |
| gtex-aging | nucleus_accumbens_basal_ganglia | no | 16,529 | 0.0061 | 0.0219 | 0.0258 | 0.0060 | 4.3x |
| gtex-aging | putamen_basal_ganglia | no | 16,158 | 0.0123 | 0.0346 | 0.0107 | 0.0114 | 0.9x |
| gtex-aging | spinal_cord_cervical_c_1 | no | 16,385 | 0.0081 | 0.0270 | 0.0216 | 0.0081 | 2.7x |
| gtex-aging | substantia_nigra | no | 16,368 | 0.0111 | 0.0342 | 0.0841 | 0.0111 | 7.6x |

## `partial_r2_switch_given_abund`

| analysis | region | adj | n | median (all) | q90 (all) | median: composition_unique | median: neither | ratio vs neither |
|---|---|---|---|---|---|---|---|---|
| brainseq-aging | caudate | no | 11,648 | 0.0087 | 0.0303 | 0.0858 | 0.0086 | 10.0x |
| brainseq-aging | dlpfc | no | 11,554 | 0.0075 | 0.0257 | 0.1102 | 0.0072 | 15.4x |
| brainseq-aging | hippocampus | no | 11,114 | 0.0079 | 0.0263 | n/a | 0.0078 | n/a |
| brainseq-sczd | caudate | no | 17,193 | 0.0017 | 0.0100 | 0.0446 | 0.0016 | 27.7x |
| gtex-aging | amygdala | no | 16,351 | 0.0098 | 0.0317 | 0.1257 | 0.0098 | 12.8x |
| gtex-aging | anterior_cingulate_cortex_ba24 | no | 16,250 | 0.0095 | 0.0330 | 0.0586 | 0.0095 | 6.2x |
| gtex-aging | caudate_basal_ganglia | no | 16,543 | 0.0056 | 0.0188 | 0.0605 | 0.0055 | 11.0x |
| gtex-aging | cerebellar_hemisphere | no | 16,140 | 0.0056 | 0.0183 | 0.0833 | 0.0056 | 14.9x |
| gtex-aging | cerebellum | no | 16,382 | 0.0064 | 0.0212 | 0.0909 | 0.0064 | 14.2x |
| gtex-aging | cortex | no | 16,480 | 0.0072 | 0.0261 | 0.0540 | 0.0071 | 7.6x |
| gtex-aging | frontal_cortex_ba9 | no | 16,362 | 0.0084 | 0.0292 | 0.0494 | 0.0083 | 6.0x |
| gtex-aging | hippocampus | no | 16,500 | 0.0080 | 0.0255 | 0.0693 | 0.0079 | 8.8x |
| gtex-aging | hypothalamus | no | 16,898 | 0.0076 | 0.0251 | 0.0696 | 0.0076 | 9.2x |
| gtex-aging | nucleus_accumbens_basal_ganglia | no | 16,529 | 0.0053 | 0.0177 | 0.0863 | 0.0053 | 16.4x |
| gtex-aging | putamen_basal_ganglia | no | 16,158 | 0.0062 | 0.0199 | 0.0839 | 0.0061 | 13.8x |
| gtex-aging | spinal_cord_cervical_c_1 | no | 16,385 | 0.0082 | 0.0264 | 0.1208 | 0.0082 | 14.6x |
| gtex-aging | substantia_nigra | no | 16,368 | 0.0107 | 0.0329 | 0.1662 | 0.0107 | 15.6x |

