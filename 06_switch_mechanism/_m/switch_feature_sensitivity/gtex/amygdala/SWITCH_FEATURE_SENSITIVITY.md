# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/amygdala** (RSEM, n=181, 33 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12063 | 0 | 0.9949 | 0.02843 | 0.9948 | 6 | 6 | 6 | 0 |
| pseudocount | pseudocount=0.25 | False | 12063 | 0 | 0.9988 | 0.01351 | 0.9991 | 6 | 6 | 6 | 0 |
| pseudocount | pseudocount=0.5 | True | 12063 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| pseudocount | pseudocount=1 | False | 12063 | 0 | 0.9984 | 0.01915 | 0.9922 | 0 | 6 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 12063 | 0 | 0.9926 | 0.04418 | 0.9787 | 0 | 6 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12063 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13053 | 0 | 1 | 0.1407 | 0.9507 | 6 | 6 | 5 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10128 | 1935 | 1 | 0.1924 | 0.8372 | 0 | 6 | 0 | 1 |
| expression | switching share>=0.1,frac>=0.05 | False | 12651 | 0 | 1 | 0.08862 | 0.9835 | 6 | 6 | 5 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10701 | 1362 | 1 | 0.1589 | 0.8997 | 2 | 6 | 2 | 1 |
| expression | legacy count>10,frac>=0.7 | False | 10031 | 3022 | 0.9182 | 0.2878 | 0.6154 | 5 | 6 | 3 | 2 |
| expression | no filter | False | 16359 | 0 | 0.8221 | 0.355 | 0.6238 | 9 | 6 | 4 | 2 |
| minor_isoform | min_usage=0 | True | 12063 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| minor_isoform | min_usage=0.01 | False | 12063 | 0 | 1 | 0 | 1 | 6 | 6 | 6 | 0 |
| minor_isoform | min_usage=0.05 | False | 11898 | 165 | 1 | 0.03967 | 0.9884 | 2 | 6 | 2 | 0 |
| minor_isoform | min_usage=0.1 | False | 10678 | 1385 | 1 | 0.1916 | 0.8324 | 2 | 6 | 2 | 1 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4681 | 0.07875 | 0.05127 | 0.517 | gtex | amygdala |
| [3, 5) | 5327 | 0.07493 | 0.04468 | 0.5234 | gtex | amygdala |
| [5, 10) | 2022 | 0.06807 | 0.0272 | 0.4585 | gtex | amygdala |
| [10, 20) | 33 | 0.06096 | 0 | 0.3636 | gtex | amygdala |
