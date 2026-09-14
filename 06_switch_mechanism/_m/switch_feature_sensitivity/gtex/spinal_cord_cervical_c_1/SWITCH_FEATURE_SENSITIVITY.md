# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/spinal_cord_cervical_c_1** (RSEM, n=204, 32 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 1.24345e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16393 | 0 | 0.9879 | 0.09955 | 0.9894 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 16393 | 0 | 0.9972 | 0.05136 | 0.9959 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 16393 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 16393 | 0 | 0.9961 | 0.06436 | 0.9455 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 16393 | 0 | 0.981 | 0.1309 | 0.8465 | 0 | 0 | 0 | 0 |
| expression | no filter | True | 16393 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | count>5,frac>=0.5 | False | 13059 | 3334 | 0.8092 | 0.32 | 0.769 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.5 | False | 11731 | 4662 | 0.7011 | 0.3524 | 0.7742 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.7 | False | 10368 | 6025 | 0.3883 | 0.4391 | 0.6721 | 0 | 0 | 0 | 0 |
| expression | count>20,frac>=0.7 | False | 8911 | 7482 | 0.3051 | 0.4607 | 0.5439 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.9 | False | 8170 | 8223 | 0.2191 | 0.4933 | 0.4682 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 16393 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 15207 | 1186 | 0.9763 | 0.2073 | 0.8952 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 13085 | 3308 | 0.7181 | 0.3635 | 0.7353 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 10994 | 5399 | 0.4456 | 0.4323 | 0.6018 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1127 | 0.05264 | 0.008873 | 0.2325 | gtex | spinal_cord_cervical_c_1 |
| [3, 5) | 2506 | 0.04919 | 0.005188 | 0.2047 | gtex | spinal_cord_cervical_c_1 |
| [5, 10) | 5438 | 0.04846 | 0.004597 | 0.2196 | gtex | spinal_cord_cervical_c_1 |
| [10, 20) | 4708 | 0.0476 | 0.004885 | 0.2918 | gtex | spinal_cord_cervical_c_1 |
| [20, 10000) | 2614 | 0.04773 | 0.002678 | 0.4457 | gtex | spinal_cord_cervical_c_1 |
