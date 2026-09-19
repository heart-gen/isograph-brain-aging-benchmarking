# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/caudate_basal_ganglia** (RSEM, n=300, 28 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12445 | 0 | 0.9945 | 0.02796 | 0.991 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 12445 | 0 | 0.9988 | 0.01543 | 0.9974 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 12445 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 12445 | 0 | 0.9983 | 0.02129 | 0.9931 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 12445 | 0 | 0.9923 | 0.04275 | 0.9754 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12445 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13501 | 0 | 1 | 0.1463 | 0.9554 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10429 | 2016 | 1 | 0.1968 | 0.825 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 13067 | 0 | 1 | 0.09048 | 0.9923 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11163 | 1282 | 1 | 0.1596 | 0.9133 | 0 | 0 | 0 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10607 | 2928 | 0.8961 | 0.2994 | 0.7326 | 0 | 0 | 0 | 0 |
| expression | no filter | False | 16551 | 0 | 0.7904 | 0.3483 | 0.6154 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 12445 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 12445 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 12298 | 147 | 1 | 0.03993 | 0.983 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 10979 | 1466 | 1 | 0.1941 | 0.8636 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4725 | 0.05178 | 0.006984 | 0.5742 | gtex | caudate_basal_ganglia |
| [3, 5) | 5532 | 0.04843 | 0.005061 | 0.5703 | gtex | caudate_basal_ganglia |
| [5, 10) | 2152 | 0.04544 | 0.003253 | 0.4995 | gtex | caudate_basal_ganglia |
| [10, 20) | 36 | 0.06577 | 0 | 0.4167 | gtex | caudate_basal_ganglia |
