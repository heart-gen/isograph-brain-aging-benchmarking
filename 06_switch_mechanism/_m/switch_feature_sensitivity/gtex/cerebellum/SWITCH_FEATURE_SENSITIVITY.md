# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/cerebellum** (RSEM, n=266, 31 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 5.60663e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16391 | 0 | 0.9896 | 0.09225 | 0.9914 | 7 | 7 | 6 | 0 |
| pseudocount | pseudocount=0.25 | False | 16391 | 0 | 0.9976 | 0.04783 | 0.9965 | 7 | 7 | 7 | 0 |
| pseudocount | pseudocount=0.5 | True | 16391 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| pseudocount | pseudocount=1 | False | 16391 | 0 | 0.9966 | 0.06241 | 0.9905 | 9 | 7 | 7 | 0 |
| pseudocount | pseudocount=2 | False | 16391 | 0 | 0.9839 | 0.1271 | 0.9746 | 10 | 7 | 7 | 0 |
| expression | no filter | True | 16391 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| expression | count>5,frac>=0.5 | False | 13499 | 2892 | 0.845 | 0.3085 | 0.8973 | 14 | 7 | 7 | 0 |
| expression | count>10,frac>=0.5 | False | 12246 | 4145 | 0.705 | 0.3466 | 0.8409 | 14 | 7 | 7 | 0 |
| expression | count>10,frac>=0.7 | False | 11092 | 5299 | 0.3335 | 0.4354 | 0.6256 | 21 | 7 | 6 | 1 |
| expression | count>20,frac>=0.7 | False | 9779 | 6612 | 0.2678 | 0.454 | 0.5618 | 23 | 7 | 6 | 1 |
| expression | count>10,frac>=0.9 | False | 9424 | 6967 | 0.1905 | 0.4949 | 0.429 | 25 | 7 | 7 | 1 |
| minor_isoform | min_usage=0 | True | 16391 | 0 | 1 | 0 | 1 | 7 | 7 | 7 | 0 |
| minor_isoform | min_usage=0.01 | False | 15279 | 1112 | 0.9817 | 0.2017 | 0.9625 | 11 | 7 | 7 | 0 |
| minor_isoform | min_usage=0.05 | False | 13324 | 3067 | 0.5607 | 0.3783 | 0.6452 | 18 | 7 | 6 | 1 |
| minor_isoform | min_usage=0.1 | False | 11311 | 5080 | 0.3139 | 0.4459 | 0.4665 | 25 | 7 | 7 | 1 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1223 | 0.05785 | 0.01881 | 0.1447 | gtex | cerebellum |
| [3, 5) | 2532 | 0.05627 | 0.01343 | 0.1047 | gtex | cerebellum |
| [5, 10) | 5396 | 0.05352 | 0.01334 | 0.1205 | gtex | cerebellum |
| [10, 20) | 4658 | 0.05731 | 0.02125 | 0.1818 | gtex | cerebellum |
| [20, 10000) | 2582 | 0.06454 | 0.02866 | 0.2924 | gtex | cerebellum |
