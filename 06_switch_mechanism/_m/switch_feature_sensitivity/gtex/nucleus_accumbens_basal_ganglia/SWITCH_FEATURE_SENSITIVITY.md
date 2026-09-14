# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/nucleus_accumbens_basal_ganglia** (RSEM, n=285, 32 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 1.58276e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16538 | 0 | 0.9883 | 0.09911 | 0.9872 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 16538 | 0 | 0.9973 | 0.05043 | 0.9954 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 16538 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 16538 | 0 | 0.9962 | 0.06149 | 0.921 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 16538 | 0 | 0.9821 | 0.1272 | 0.84 | 0 | 0 | 0 | 0 |
| expression | no filter | True | 16538 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | count>5,frac>=0.5 | False | 13363 | 3175 | 0.8346 | 0.3061 | 0.7167 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.5 | False | 12135 | 4403 | 0.7498 | 0.3389 | 0.6103 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.7 | False | 10717 | 5821 | 0.4714 | 0.4131 | 0.3717 | 0 | 0 | 0 | 0 |
| expression | count>20,frac>=0.7 | False | 9209 | 7329 | 0.3999 | 0.4353 | 0.3357 | 0 | 0 | 0 | 0 |
| expression | count>10,frac>=0.9 | False | 8206 | 8332 | 0.2905 | 0.4831 | 0.07021 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 16538 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 15213 | 1325 | 0.9729 | 0.2082 | 0.9437 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 13086 | 3452 | 0.7624 | 0.3599 | 0.6618 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 11026 | 5512 | 0.5283 | 0.4282 | 0.2517 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1143 | 0.05785 | 0.005249 | 0.2283 | gtex | nucleus_accumbens_basal_ganglia |
| [3, 5) | 2542 | 0.05778 | 0.009835 | 0.2349 | gtex | nucleus_accumbens_basal_ganglia |
| [5, 10) | 5475 | 0.0572 | 0.009315 | 0.2612 | gtex | nucleus_accumbens_basal_ganglia |
| [10, 20) | 4711 | 0.06686 | 0.01528 | 0.3725 | gtex | nucleus_accumbens_basal_ganglia |
| [20, 10000) | 2667 | 0.07703 | 0.018 | 0.5017 | gtex | nucleus_accumbens_basal_ganglia |
