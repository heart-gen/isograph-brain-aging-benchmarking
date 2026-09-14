# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/anterior_cingulate_cortex_ba24** (RSEM, n=233, 42 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 5.32907e-15 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16258 | 0 | 0.9884 | 0.09823 | 0.9952 | 19 | 20 | 19 | 0 |
| pseudocount | pseudocount=0.25 | False | 16258 | 0 | 0.9973 | 0.04927 | 0.997 | 20 | 20 | 20 | 0 |
| pseudocount | pseudocount=0.5 | True | 16258 | 0 | 1 | 0 | 1 | 20 | 20 | 20 | 0 |
| pseudocount | pseudocount=1 | False | 16258 | 0 | 0.9963 | 0.06003 | 0.9764 | 20 | 20 | 19 | 0 |
| pseudocount | pseudocount=2 | False | 16258 | 0 | 0.982 | 0.1244 | 0.9551 | 20 | 20 | 19 | 0 |
| expression | no filter | True | 16258 | 0 | 1 | 0 | 1 | 20 | 20 | 20 | 0 |
| expression | count>5,frac>=0.5 | False | 13077 | 3181 | 0.8367 | 0.3112 | 0.911 | 22 | 20 | 19 | 0 |
| expression | count>10,frac>=0.5 | False | 11801 | 4457 | 0.7425 | 0.3441 | 0.8972 | 23 | 20 | 19 | 1 |
| expression | count>10,frac>=0.7 | False | 10403 | 5855 | 0.4271 | 0.4214 | 0.815 | 18 | 20 | 16 | 1 |
| expression | count>20,frac>=0.7 | False | 8951 | 7307 | 0.3563 | 0.4367 | 0.7813 | 16 | 20 | 14 | 1 |
| expression | count>10,frac>=0.9 | False | 8215 | 8043 | 0.2527 | 0.4789 | 0.6869 | 14 | 20 | 13 | 1 |
| minor_isoform | min_usage=0 | True | 16258 | 0 | 1 | 0 | 1 | 20 | 20 | 20 | 0 |
| minor_isoform | min_usage=0.01 | False | 14878 | 1380 | 0.9746 | 0.2046 | 0.9667 | 20 | 20 | 19 | 0 |
| minor_isoform | min_usage=0.05 | False | 12756 | 3502 | 0.7571 | 0.3527 | 0.8601 | 20 | 20 | 18 | 0 |
| minor_isoform | min_usage=0.1 | False | 10705 | 5553 | 0.4725 | 0.4261 | 0.8044 | 18 | 20 | 15 | 1 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1122 | 0.07553 | 0.05348 | 0.2701 | gtex | anterior_cingulate_cortex_ba24 |
| [3, 5) | 2469 | 0.06606 | 0.0482 | 0.2592 | gtex | anterior_cingulate_cortex_ba24 |
| [5, 10) | 5377 | 0.06486 | 0.04315 | 0.2784 | gtex | anterior_cingulate_cortex_ba24 |
| [10, 20) | 4658 | 0.07316 | 0.05689 | 0.3669 | gtex | anterior_cingulate_cortex_ba24 |
| [20, 10000) | 2632 | 0.08731 | 0.103 | 0.4696 | gtex | anterior_cingulate_cortex_ba24 |
