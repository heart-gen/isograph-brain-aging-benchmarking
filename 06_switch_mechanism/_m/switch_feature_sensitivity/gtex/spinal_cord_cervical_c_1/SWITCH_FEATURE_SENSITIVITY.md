# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/spinal_cord_cervical_c_1** (RSEM, n=204, 27 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12326 | 0 | 0.9948 | 0.02556 | 0.9836 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 12326 | 0 | 0.9988 | 0.01298 | 0.9919 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 12326 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 12326 | 0 | 0.9984 | 0.01801 | 0.9924 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 12326 | 0 | 0.9924 | 0.03992 | 0.9447 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12326 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13317 | 0 | 1 | 0.1518 | 0.8382 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10281 | 2045 | 1 | 0.1967 | 0.6141 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 12897 | 0 | 1 | 0.08851 | 0.9368 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10966 | 1360 | 1 | 0.1566 | 0.8327 | 0 | 0 | 0 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10368 | 3022 | 0.8579 | 0.3097 | 0.5102 | 0 | 0 | 0 | 0 |
| expression | no filter | False | 16393 | 0 | 0.7936 | 0.3552 | 0.5085 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 12326 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 12326 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 12184 | 142 | 1 | 0.03931 | 0.9812 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 10910 | 1416 | 1 | 0.1902 | 0.6657 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4703 | 0.04419 | 0.00404 | 0.2819 | gtex | spinal_cord_cervical_c_1 |
| [3, 5) | 5441 | 0.0456 | 0.003124 | 0.3251 | gtex | spinal_cord_cervical_c_1 |
| [5, 10) | 2144 | 0.04745 | 0.002799 | 0.2738 | gtex | spinal_cord_cervical_c_1 |
| [10, 20) | 38 | 0.06067 | 0 | 0.2632 | gtex | spinal_cord_cervical_c_1 |
