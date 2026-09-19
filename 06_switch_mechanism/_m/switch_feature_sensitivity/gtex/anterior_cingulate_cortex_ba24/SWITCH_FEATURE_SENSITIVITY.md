# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/anterior_cingulate_cortex_ba24** (RSEM, n=233, 25 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12212 | 0 | 0.9949 | 0.0271 | 0.9934 | 3 | 4 | 3 | 0 |
| pseudocount | pseudocount=0.25 | False | 12212 | 0 | 0.9988 | 0.01286 | 0.9992 | 3 | 4 | 3 | 0 |
| pseudocount | pseudocount=0.5 | True | 12212 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=1 | False | 12212 | 0 | 0.9984 | 0.01883 | 0.9969 | 4 | 4 | 4 | 0 |
| pseudocount | pseudocount=2 | False | 12212 | 0 | 0.9927 | 0.04209 | 0.9836 | 3 | 4 | 3 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12212 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13234 | 0 | 1 | 0.1442 | 0.9542 | 4 | 4 | 4 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10219 | 1993 | 1 | 0.1887 | 0.8481 | 0 | 4 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 12809 | 0 | 1 | 0.09171 | 0.9833 | 2 | 4 | 2 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10851 | 1361 | 1 | 0.1514 | 0.8834 | 0 | 4 | 0 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10403 | 2913 | 0.9174 | 0.2975 | 0.8 | 0 | 4 | 0 | 1 |
| expression | no filter | False | 16258 | 0 | 0.8251 | 0.3405 | 0.7037 | 3 | 4 | 1 | 1 |
| minor_isoform | min_usage=0 | True | 12212 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.01 | False | 12212 | 0 | 1 | 0 | 1 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.05 | False | 12047 | 165 | 1 | 0.03943 | 0.9959 | 4 | 4 | 4 | 0 |
| minor_isoform | min_usage=0.1 | False | 10791 | 1421 | 1 | 0.1892 | 0.8526 | 0 | 4 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4736 | 0.07337 | 0.06292 | 0.568 | gtex | anterior_cingulate_cortex_ba24 |
| [3, 5) | 5371 | 0.0746 | 0.05828 | 0.529 | gtex | anterior_cingulate_cortex_ba24 |
| [5, 10) | 2074 | 0.06795 | 0.05063 | 0.4769 | gtex | anterior_cingulate_cortex_ba24 |
| [10, 20) | 31 | 0.07761 | 0.129 | 0.5161 | gtex | anterior_cingulate_cortex_ba24 |
