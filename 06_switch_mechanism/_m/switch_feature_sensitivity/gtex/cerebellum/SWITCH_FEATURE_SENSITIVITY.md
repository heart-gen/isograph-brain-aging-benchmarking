# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/cerebellum** (RSEM, n=266, 21 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12501 | 0 | 0.996 | 0.02128 | 0.9982 | 8 | 7 | 7 | 0 |
| pseudocount | pseudocount=0.25 | False | 12501 | 0 | 0.9991 | 0.01048 | 0.9995 | 7 | 7 | 7 | 0 |
| pseudocount | pseudocount=0.5 | True | 12501 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| pseudocount | pseudocount=1 | False | 12501 | 0 | 0.9988 | 0.01464 | 0.9989 | 7 | 7 | 7 | 0 |
| pseudocount | pseudocount=2 | False | 12501 | 0 | 0.9943 | 0.03064 | 0.9947 | 7 | 7 | 7 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12501 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13532 | 0 | 1 | 0.1636 | 0.9579 | 9 | 7 | 6 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10240 | 2261 | 1 | 0.2079 | 0.8167 | 15 | 7 | 7 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 13032 | 0 | 1 | 0.08823 | 0.9936 | 6 | 7 | 6 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11254 | 1247 | 1 | 0.1441 | 0.8764 | 12 | 7 | 7 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 11092 | 2674 | 0.7686 | 0.3237 | 0.7845 | 13 | 7 | 7 | 0 |
| expression | no filter | False | 16391 | 0 | 0.6424 | 0.3642 | 0.7629 | 16 | 7 | 7 | 0 |
| minor_isoform | min_usage=0 | True | 12501 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| minor_isoform | min_usage=0.01 | False | 12501 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| minor_isoform | min_usage=0.05 | False | 12396 | 105 | 1 | 0.03308 | 0.9769 | 8 | 7 | 7 | 0 |
| minor_isoform | min_usage=0.1 | False | 11294 | 1207 | 1 | 0.175 | 0.8504 | 14 | 7 | 7 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4810 | 0.06733 | 0.03035 | 0.4898 | gtex | cerebellum |
| [3, 5) | 5689 | 0.06387 | 0.02848 | 0.4305 | gtex | cerebellum |
| [5, 10) | 1979 | 0.0601 | 0.02072 | 0.3557 | gtex | cerebellum |
| [10, 20) | 23 | 0.05714 | 0 | 0.04348 | gtex | cerebellum |
