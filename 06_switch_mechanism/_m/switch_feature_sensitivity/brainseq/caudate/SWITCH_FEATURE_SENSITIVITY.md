# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **brainseq/caudate** (Salmon, n=238, 46 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 10.0, "min_fraction": 0.7, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 11648 | 0 | 0.9946 | 0.02576 | 0.9975 | 18 | 20 | 18 | 0 |
| pseudocount | pseudocount=0.25 | False | 11648 | 0 | 0.9988 | 0.01288 | 0.9992 | 20 | 20 | 20 | 0 |
| pseudocount | pseudocount=0.5 | True | 11648 | 0 | 1 | 0 | 1 | 20 | 20 | 20 | 0 |
| pseudocount | pseudocount=1 | False | 11648 | 0 | 0.9983 | 0.01356 | 0.9981 | 20 | 20 | 20 | 0 |
| pseudocount | pseudocount=2 | False | 11648 | 0 | 0.992 | 0.03666 | 0.9924 | 21 | 20 | 20 | 0 |
| expression | count>5,frac>=0.5 | False | 14080 | 0 | 0.8417 | 0.3125 | 0.8888 | 23 | 20 | 18 | 0 |
| expression | count>10,frac>=0.5 | False | 12807 | 0 | 1 | 0.221 | 0.9315 | 23 | 20 | 19 | 0 |
| expression | count>10,frac>=0.7 | True | 11648 | 0 | 1 | 0 | 1 | 20 | 20 | 20 | 0 |
| expression | count>20,frac>=0.7 | False | 10167 | 1481 | 1 | 0.1606 | 0.97 | 21 | 20 | 19 | 0 |
| expression | count>10,frac>=0.9 | False | 9704 | 1944 | 0.9152 | 0.281 | 0.9033 | 21 | 20 | 18 | 0 |
| minor_isoform | min_usage=0 | True | 11648 | 0 | 1 | 0 | 1 | 20 | 20 | 20 | 0 |
| minor_isoform | min_usage=0.01 | False | 11551 | 97 | 1 | 0.02987 | 0.9961 | 20 | 20 | 20 | 0 |
| minor_isoform | min_usage=0.05 | False | 10450 | 1198 | 1 | 0.1725 | 0.9627 | 20 | 20 | 18 | 0 |
| minor_isoform | min_usage=0.1 | False | 8921 | 2727 | 0.9591 | 0.2671 | 0.8955 | 20 | 20 | 15 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 3253 | 0.06566 | 0.05287 | 0.2804 | brainseq | caudate |
| [3, 5) | 4271 | 0.06908 | 0.05011 | 0.3329 | brainseq | caudate |
| [5, 10) | 3339 | 0.07143 | 0.04073 | 0.3687 | brainseq | caudate |
| [10, 20) | 716 | 0.08424 | 0.05168 | 0.4609 | brainseq | caudate |
| [20, 10000) | 69 | 0.1023 | 0.05797 | 0.6232 | brainseq | caudate |
