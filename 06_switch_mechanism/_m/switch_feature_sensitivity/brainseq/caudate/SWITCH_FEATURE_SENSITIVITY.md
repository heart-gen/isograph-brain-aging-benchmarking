# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **brainseq/caudate** (Salmon, n=238, 28 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 13177 | 0 | 0.9968 | 0.02292 | 0.9971 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=0.25 | False | 13177 | 0 | 0.9993 | 0.01184 | 0.9977 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=0.5 | True | 13177 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=1 | False | 13177 | 0 | 0.9989 | 0.01427 | 0.998 | 5 | 4 | 4 | 0 |
| pseudocount | pseudocount=2 | False | 13177 | 0 | 0.9947 | 0.03324 | 0.9959 | 6 | 4 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 13177 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 14415 | 0 | 1 | 0.1697 | 0.9777 | 9 | 4 | 4 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10303 | 2874 | 0.9999 | 0.2172 | 0.9435 | 5 | 4 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 13940 | 0 | 1 | 0.1004 | 0.9926 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11497 | 1680 | 1 | 0.1681 | 0.9912 | 5 | 4 | 4 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 11648 | 3011 | 0.7357 | 0.3439 | 0.9622 | 6 | 4 | 4 | 0 |
| expression | no filter | False | 17111 | 0 | 0.7907 | 0.3338 | 0.9404 | 5 | 4 | 3 | 0 |
| minor_isoform | min_usage=0 | True | 13177 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.01 | False | 13177 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.05 | False | 13047 | 130 | 1 | 0.04124 | 0.9996 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.1 | False | 11539 | 1638 | 1 | 0.1974 | 0.972 | 6 | 4 | 4 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4820 | 0.06159 | 0.04668 | 0.3359 | brainseq | caudate |
| [3, 5) | 5864 | 0.06083 | 0.03871 | 0.2729 | brainseq | caudate |
| [5, 10) | 2446 | 0.05682 | 0.01431 | 0.1721 | brainseq | caudate |
| [10, 20) | 47 | 0.0654 | 0 | 0.08511 | brainseq | caudate |
