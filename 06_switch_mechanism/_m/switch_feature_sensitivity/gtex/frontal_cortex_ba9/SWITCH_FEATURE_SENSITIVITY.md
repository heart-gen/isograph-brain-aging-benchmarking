# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/frontal_cortex_ba9** (RSEM, n=269, 48 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 4.42979e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16371 | 0 | 0.9888 | 0.09981 | 0.9987 | 41 | 40 | 40 | 0 |
| pseudocount | pseudocount=0.25 | False | 16371 | 0 | 0.9974 | 0.05094 | 0.9997 | 41 | 40 | 40 | 0 |
| pseudocount | pseudocount=0.5 | True | 16371 | 0 | 1 | 0 | 1 | 40 | 40 | 40 | 0 |
| pseudocount | pseudocount=1 | False | 16371 | 0 | 0.9964 | 0.06316 | 0.9904 | 38 | 40 | 38 | 1 |
| pseudocount | pseudocount=2 | False | 16371 | 0 | 0.9828 | 0.125 | 0.9762 | 34 | 40 | 34 | 1 |
| expression | no filter | True | 16371 | 0 | 1 | 0 | 1 | 40 | 40 | 40 | 0 |
| expression | count>5,frac>=0.5 | False | 13342 | 3029 | 0.8485 | 0.307 | 0.9567 | 30 | 40 | 30 | 2 |
| expression | count>10,frac>=0.5 | False | 12101 | 4270 | 0.7533 | 0.3427 | 0.9067 | 29 | 40 | 29 | 3 |
| expression | count>10,frac>=0.7 | False | 10756 | 5615 | 0.4258 | 0.4264 | 0.8099 | 25 | 40 | 25 | 10 |
| expression | count>20,frac>=0.7 | False | 9352 | 7019 | 0.3549 | 0.4431 | 0.7442 | 22 | 40 | 22 | 11 |
| expression | count>10,frac>=0.9 | False | 8630 | 7741 | 0.2523 | 0.4766 | 0.5833 | 22 | 40 | 22 | 12 |
| minor_isoform | min_usage=0 | True | 16371 | 0 | 1 | 0 | 1 | 40 | 40 | 40 | 0 |
| minor_isoform | min_usage=0.01 | False | 15039 | 1332 | 0.9758 | 0.2153 | 0.9912 | 34 | 40 | 34 | 0 |
| minor_isoform | min_usage=0.05 | False | 12897 | 3474 | 0.7289 | 0.3617 | 0.8834 | 27 | 40 | 27 | 8 |
| minor_isoform | min_usage=0.1 | False | 10847 | 5524 | 0.4647 | 0.425 | 0.7314 | 27 | 40 | 27 | 11 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1146 | 0.07256 | 0.05934 | 0.4363 | gtex | frontal_cortex_ba9 |
| [3, 5) | 2500 | 0.06606 | 0.048 | 0.4188 | gtex | frontal_cortex_ba9 |
| [5, 10) | 5396 | 0.07035 | 0.05541 | 0.4375 | gtex | frontal_cortex_ba9 |
| [10, 20) | 4680 | 0.07989 | 0.09316 | 0.5013 | gtex | frontal_cortex_ba9 |
| [20, 10000) | 2649 | 0.09982 | 0.1495 | 0.5515 | gtex | frontal_cortex_ba9 |
