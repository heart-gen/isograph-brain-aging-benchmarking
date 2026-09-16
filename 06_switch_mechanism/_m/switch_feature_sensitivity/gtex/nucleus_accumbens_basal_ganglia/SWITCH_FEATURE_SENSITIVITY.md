# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/nucleus_accumbens_basal_ganglia** (RSEM, n=285, 21 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12513 | 0 | 0.9946 | 0.02821 | 0.9969 | 6 | 5 | 5 | 0 |
| pseudocount | pseudocount=0.25 | False | 12513 | 0 | 0.9988 | 0.01391 | 0.9984 | 6 | 5 | 5 | 0 |
| pseudocount | pseudocount=0.5 | True | 12513 | 0 | 1 | 0 | 1 | 5 | 5 | 5 | 0 |
| pseudocount | pseudocount=1 | False | 12513 | 0 | 0.9984 | 0.01838 | 0.9973 | 0 | 5 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 12513 | 0 | 0.9925 | 0.04228 | 0.9902 | 0 | 5 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12513 | 0 | 1 | 0 | 1 | 5 | 5 | 5 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13534 | 0 | 1 | 0.1439 | 0.9867 | 0 | 5 | 0 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10515 | 1998 | 1 | 0.1916 | 0.9578 | 0 | 5 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 13087 | 0 | 1 | 0.08879 | 0.9925 | 1 | 5 | 1 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11173 | 1340 | 1 | 0.1548 | 0.9617 | 1 | 5 | 1 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10717 | 2853 | 0.9351 | 0.2925 | 0.8917 | 0 | 5 | 0 | 0 |
| expression | no filter | False | 16538 | 0 | 0.8206 | 0.3476 | 0.8517 | 0 | 5 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 12513 | 0 | 1 | 0 | 1 | 5 | 5 | 5 | 0 |
| minor_isoform | min_usage=0.01 | False | 12513 | 0 | 1 | 0 | 1 | 5 | 5 | 5 | 0 |
| minor_isoform | min_usage=0.05 | False | 12363 | 150 | 1 | 0.03947 | 0.9954 | 0 | 5 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 11080 | 1433 | 1 | 0.1913 | 0.8879 | 0 | 5 | 0 | 1 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4832 | 0.07049 | 0.01718 | 0.3148 | gtex | nucleus_accumbens_basal_ganglia |
| [3, 5) | 5482 | 0.06796 | 0.01332 | 0.3457 | gtex | nucleus_accumbens_basal_ganglia |
| [5, 10) | 2163 | 0.059 | 0.009246 | 0.3107 | gtex | nucleus_accumbens_basal_ganglia |
| [10, 20) | 36 | 0.04793 | 0 | 0.1667 | gtex | nucleus_accumbens_basal_ganglia |
