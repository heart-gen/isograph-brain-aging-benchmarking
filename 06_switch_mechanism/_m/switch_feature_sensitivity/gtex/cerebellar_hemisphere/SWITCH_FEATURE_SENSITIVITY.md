# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/cerebellar_hemisphere** (RSEM, n=277, 30 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 1.17684e-13 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16149 | 0 | 0.9891 | 0.09518 | 0.9945 | 11 | 12 | 11 | 0 |
| pseudocount | pseudocount=0.25 | False | 16149 | 0 | 0.9975 | 0.04737 | 0.9977 | 12 | 12 | 12 | 0 |
| pseudocount | pseudocount=0.5 | True | 16149 | 0 | 1 | 0 | 1 | 12 | 12 | 12 | 0 |
| pseudocount | pseudocount=1 | False | 16149 | 0 | 0.9965 | 0.06062 | 0.9959 | 11 | 12 | 11 | 0 |
| pseudocount | pseudocount=2 | False | 16149 | 0 | 0.9833 | 0.1208 | 0.9818 | 10 | 12 | 10 | 0 |
| expression | no filter | True | 16149 | 0 | 1 | 0 | 1 | 12 | 12 | 12 | 0 |
| expression | count>5,frac>=0.5 | False | 13370 | 2779 | 0.8668 | 0.3 | 0.924 | 9 | 12 | 9 | 0 |
| expression | count>10,frac>=0.5 | False | 12206 | 3943 | 0.7891 | 0.3379 | 0.8876 | 8 | 12 | 8 | 0 |
| expression | count>10,frac>=0.7 | False | 10956 | 5193 | 0.4724 | 0.4194 | 0.7121 | 8 | 12 | 6 | 1 |
| expression | count>20,frac>=0.7 | False | 9687 | 6462 | 0.408 | 0.4437 | 0.6316 | 14 | 12 | 6 | 3 |
| expression | count>10,frac>=0.9 | False | 9005 | 7144 | 0.2864 | 0.4888 | 0.2949 | 19 | 12 | 9 | 5 |
| minor_isoform | min_usage=0 | True | 16149 | 0 | 1 | 0 | 1 | 12 | 12 | 12 | 0 |
| minor_isoform | min_usage=0.01 | False | 14996 | 1153 | 0.9771 | 0.2042 | 0.9703 | 10 | 12 | 10 | 0 |
| minor_isoform | min_usage=0.05 | False | 13034 | 3115 | 0.7015 | 0.3729 | 0.7479 | 8 | 12 | 6 | 2 |
| minor_isoform | min_usage=0.1 | False | 11027 | 5122 | 0.4554 | 0.4334 | 0.5032 | 12 | 12 | 6 | 4 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1188 | 0.06438 | 0.02104 | 0.1978 | gtex | cerebellar_hemisphere |
| [3, 5) | 2480 | 0.06416 | 0.02419 | 0.2206 | gtex | cerebellar_hemisphere |
| [5, 10) | 5320 | 0.06618 | 0.02744 | 0.2797 | gtex | cerebellar_hemisphere |
| [10, 20) | 4610 | 0.07843 | 0.04534 | 0.3885 | gtex | cerebellar_hemisphere |
| [20, 10000) | 2551 | 0.1015 | 0.07409 | 0.4833 | gtex | cerebellar_hemisphere |
