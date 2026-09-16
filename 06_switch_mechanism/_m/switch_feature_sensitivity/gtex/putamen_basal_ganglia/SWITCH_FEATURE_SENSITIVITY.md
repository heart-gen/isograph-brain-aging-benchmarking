# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/putamen_basal_ganglia** (RSEM, n=254, 22 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12121 | 0 | 0.9944 | 0.02871 | 0.993 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 12121 | 0 | 0.9987 | 0.01378 | 0.9946 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 12121 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 12121 | 0 | 0.9983 | 0.01922 | 0.998 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 12121 | 0 | 0.9921 | 0.04348 | 0.9922 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12121 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13140 | 0 | 1 | 0.1469 | 0.9812 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10195 | 1926 | 1 | 0.1972 | 0.8555 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 12728 | 0 | 1 | 0.09389 | 0.9883 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10777 | 1344 | 1 | 0.1569 | 0.8751 | 1 | 0 | 0 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10105 | 2990 | 0.9281 | 0.2909 | 0.801 | 0 | 0 | 0 | 0 |
| expression | no filter | False | 16166 | 0 | 0.8026 | 0.3556 | 0.702 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 12121 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 12121 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 11970 | 151 | 1 | 0.03818 | 0.9759 | 1 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 10665 | 1456 | 1 | 0.1975 | 0.7973 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4719 | 0.06231 | 0.01229 | 0.3132 | gtex | putamen_basal_ganglia |
| [3, 5) | 5322 | 0.05722 | 0.01033 | 0.3797 | gtex | putamen_basal_ganglia |
| [5, 10) | 2050 | 0.05512 | 0.005366 | 0.3132 | gtex | putamen_basal_ganglia |
| [10, 20) | 30 | 0.04791 | 0 | 0.3 | gtex | putamen_basal_ganglia |
