# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/frontal_cortex_ba9** (RSEM, n=269, 25 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12345 | 0 | 0.9951 | 0.02916 | 0.9939 | 9 | 10 | 9 | 0 |
| pseudocount | pseudocount=0.25 | False | 12345 | 0 | 0.9989 | 0.01409 | 0.9959 | 9 | 10 | 9 | 0 |
| pseudocount | pseudocount=0.5 | True | 12345 | 0 | 1 | 0 | 1 | 10 | 10 | 10 | 0 |
| pseudocount | pseudocount=1 | False | 12345 | 0 | 0.9985 | 0.01952 | 0.9993 | 10 | 10 | 10 | 0 |
| pseudocount | pseudocount=2 | False | 12345 | 0 | 0.9932 | 0.03985 | 0.9962 | 8 | 10 | 8 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12345 | 0 | 1 | 0 | 1 | 10 | 10 | 10 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13413 | 0 | 1 | 0.1473 | 0.9908 | 10 | 10 | 9 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10264 | 2081 | 1 | 0.1903 | 0.9106 | 3 | 10 | 3 | 1 |
| expression | switching share>=0.1,frac>=0.05 | False | 12932 | 0 | 1 | 0.08692 | 0.9941 | 10 | 10 | 10 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10976 | 1369 | 1 | 0.154 | 0.9764 | 7 | 10 | 7 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10756 | 2780 | 0.8986 | 0.3004 | 0.9113 | 7 | 10 | 7 | 0 |
| expression | no filter | False | 16371 | 0 | 0.8111 | 0.3469 | 0.9442 | 3 | 10 | 3 | 0 |
| minor_isoform | min_usage=0 | True | 12345 | 0 | 1 | 0 | 1 | 10 | 10 | 10 | 0 |
| minor_isoform | min_usage=0.01 | False | 12345 | 0 | 1 | 0 | 1 | 10 | 10 | 10 | 0 |
| minor_isoform | min_usage=0.05 | False | 12180 | 165 | 1 | 0.04007 | 0.9974 | 12 | 10 | 10 | 0 |
| minor_isoform | min_usage=0.1 | False | 10923 | 1422 | 1 | 0.1894 | 0.9147 | 3 | 10 | 3 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4734 | 0.08444 | 0.0978 | 0.4081 | gtex | frontal_cortex_ba9 |
| [3, 5) | 5482 | 0.07976 | 0.08172 | 0.4172 | gtex | frontal_cortex_ba9 |
| [5, 10) | 2095 | 0.07238 | 0.07399 | 0.3661 | gtex | frontal_cortex_ba9 |
| [10, 20) | 34 | 0.06065 | 0.08824 | 0.2059 | gtex | frontal_cortex_ba9 |
