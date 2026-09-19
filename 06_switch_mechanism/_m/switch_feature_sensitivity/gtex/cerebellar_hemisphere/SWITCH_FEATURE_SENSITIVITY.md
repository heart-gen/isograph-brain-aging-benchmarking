# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/cerebellar_hemisphere** (RSEM, n=277, 24 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12438 | 0 | 0.9953 | 0.02621 | 0.9977 | 8 | 8 | 8 | 0 |
| pseudocount | pseudocount=0.25 | False | 12438 | 0 | 0.999 | 0.01335 | 0.9997 | 8 | 8 | 8 | 0 |
| pseudocount | pseudocount=0.5 | True | 12438 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| pseudocount | pseudocount=1 | False | 12438 | 0 | 0.9986 | 0.01737 | 0.9931 | 8 | 8 | 8 | 0 |
| pseudocount | pseudocount=2 | False | 12438 | 0 | 0.9935 | 0.03642 | 0.9873 | 8 | 8 | 7 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12438 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13404 | 0 | 1 | 0.1533 | 0.9778 | 8 | 8 | 7 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10353 | 2085 | 1 | 0.2021 | 0.8784 | 9 | 8 | 5 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 13002 | 0 | 1 | 0.09334 | 0.9916 | 7 | 8 | 7 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 11114 | 1324 | 1 | 0.1534 | 0.9029 | 7 | 8 | 5 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10956 | 2663 | 0.8591 | 0.3026 | 0.8534 | 10 | 8 | 5 | 1 |
| expression | no filter | False | 16149 | 0 | 0.7694 | 0.3595 | 0.8844 | 12 | 8 | 7 | 0 |
| minor_isoform | min_usage=0 | True | 12438 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| minor_isoform | min_usage=0.01 | False | 12438 | 0 | 1 | 0 | 1 | 8 | 8 | 8 | 0 |
| minor_isoform | min_usage=0.05 | False | 12303 | 135 | 1 | 0.03479 | 0.9941 | 8 | 8 | 8 | 0 |
| minor_isoform | min_usage=0.1 | False | 11093 | 1345 | 1 | 0.1819 | 0.8532 | 7 | 8 | 5 | 1 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4705 | 0.07793 | 0.03401 | 0.6017 | gtex | cerebellar_hemisphere |
| [3, 5) | 5675 | 0.07758 | 0.03824 | 0.5864 | gtex | cerebellar_hemisphere |
| [5, 10) | 2032 | 0.07385 | 0.03642 | 0.5477 | gtex | cerebellar_hemisphere |
| [10, 20) | 26 | 0.05937 | 0.03846 | 0.3077 | gtex | cerebellar_hemisphere |
