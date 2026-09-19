# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/hippocampus** (RSEM, n=255, 30 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12334 | 0 | 0.9944 | 0.03073 | 0.9958 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=0.25 | False | 12334 | 0 | 0.9988 | 0.0163 | 0.9989 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=0.5 | True | 12334 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=1 | False | 12334 | 0 | 0.9983 | 0.01849 | 0.9958 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=2 | False | 12334 | 0 | 0.9921 | 0.04362 | 0.9818 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12334 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13322 | 0 | 1 | 0.1425 | 0.9754 | 3 | 4 | 3 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10458 | 1876 | 1 | 0.1935 | 0.8783 | 4 | 4 | 3 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 12955 | 0 | 1 | 0.09186 | 0.9868 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10985 | 1349 | 1 | 0.1571 | 0.9313 | 3 | 4 | 3 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10180 | 3106 | 0.9248 | 0.2911 | 0.7726 | 4 | 4 | 3 | 1 |
| expression | no filter | False | 16508 | 0 | 0.8166 | 0.3468 | 0.6367 | 10 | 4 | 3 | 1 |
| minor_isoform | min_usage=0 | True | 12334 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.01 | False | 12334 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.05 | False | 12181 | 153 | 1 | 0.04195 | 0.9926 | 3 | 4 | 3 | 0 |
| minor_isoform | min_usage=0.1 | False | 10882 | 1452 | 1 | 0.1984 | 0.8618 | 4 | 4 | 3 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4663 | 0.07041 | 0.04697 | 0.5818 | gtex | hippocampus |
| [3, 5) | 5488 | 0.06752 | 0.04027 | 0.5528 | gtex | hippocampus |
| [5, 10) | 2143 | 0.06976 | 0.04106 | 0.4998 | gtex | hippocampus |
| [10, 20) | 40 | 0.07111 | 0.025 | 0.35 | gtex | hippocampus |
