# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/hypothalamus** (RSEM, n=257, 27 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12670 | 0 | 0.9945 | 0.02644 | 0.9971 | 10 | 9 | 9 | 0 |
| pseudocount | pseudocount=0.25 | False | 12670 | 0 | 0.9987 | 0.01334 | 0.9987 | 9 | 9 | 9 | 0 |
| pseudocount | pseudocount=0.5 | True | 12670 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| pseudocount | pseudocount=1 | False | 12670 | 0 | 0.9983 | 0.01665 | 0.9971 | 8 | 9 | 8 | 0 |
| pseudocount | pseudocount=2 | False | 12670 | 0 | 0.9921 | 0.03994 | 0.9882 | 7 | 9 | 7 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12670 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13701 | 0 | 1 | 0.1455 | 0.9772 | 6 | 9 | 6 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10642 | 2028 | 1 | 0.198 | 0.6275 | 0 | 9 | 0 | 2 |
| expression | switching share>=0.1,frac>=0.05 | False | 13281 | 0 | 1 | 0.09211 | 0.9867 | 8 | 9 | 8 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11259 | 1411 | 1 | 0.155 | 0.8936 | 1 | 9 | 1 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10771 | 2989 | 0.93 | 0.2918 | 0.6913 | 1 | 9 | 1 | 1 |
| expression | no filter | False | 16907 | 0 | 0.7928 | 0.3442 | 0.6263 | 0 | 9 | 0 | 2 |
| minor_isoform | min_usage=0 | True | 12670 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| minor_isoform | min_usage=0.01 | False | 12670 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| minor_isoform | min_usage=0.05 | False | 12514 | 156 | 1 | 0.04091 | 0.9894 | 7 | 9 | 7 | 0 |
| minor_isoform | min_usage=0.1 | False | 11222 | 1448 | 1 | 0.1943 | 0.7254 | 0 | 9 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4883 | 0.07852 | 0.06861 | 0.5486 | gtex | hypothalamus |
| [3, 5) | 5612 | 0.07029 | 0.04241 | 0.5203 | gtex | hypothalamus |
| [5, 10) | 2142 | 0.06768 | 0.03175 | 0.4659 | gtex | hypothalamus |
| [10, 20) | 33 | 0.04014 | 0 | 0.1818 | gtex | hypothalamus |
