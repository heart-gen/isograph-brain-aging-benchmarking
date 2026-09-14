# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/substantia_nigra** (RSEM, n=183, 42 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "min_count": 0.0, "min_fraction": 0.0, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 1.5099e-14 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 16377 | 0 | 0.9873 | 0.1005 | 0.9806 | 8 | 9 | 8 | 0 |
| pseudocount | pseudocount=0.25 | False | 16377 | 0 | 0.9971 | 0.05233 | 0.9808 | 8 | 9 | 8 | 0 |
| pseudocount | pseudocount=0.5 | True | 16377 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| pseudocount | pseudocount=1 | False | 16377 | 0 | 0.9958 | 0.06265 | 0.942 | 11 | 9 | 9 | 0 |
| pseudocount | pseudocount=2 | False | 16377 | 0 | 0.98 | 0.1272 | 0.8699 | 9 | 9 | 8 | 0 |
| expression | no filter | True | 16377 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| expression | count>5,frac>=0.5 | False | 12918 | 3459 | 0.8029 | 0.3184 | 0.7883 | 9 | 9 | 7 | 0 |
| expression | count>10,frac>=0.5 | False | 11618 | 4759 | 0.6963 | 0.3574 | 0.708 | 9 | 9 | 7 | 0 |
| expression | count>10,frac>=0.7 | False | 10184 | 6193 | 0.3827 | 0.4343 | 0.5909 | 15 | 9 | 7 | 0 |
| expression | count>20,frac>=0.7 | False | 8757 | 7620 | 0.3114 | 0.456 | 0.4953 | 12 | 9 | 7 | 1 |
| expression | count>10,frac>=0.9 | False | 8027 | 8350 | 0.2251 | 0.4979 | 0.2893 | 6 | 9 | 5 | 1 |
| minor_isoform | min_usage=0 | True | 16377 | 0 | 1 | 0 | 1 | 9 | 9 | 9 | 0 |
| minor_isoform | min_usage=0.01 | False | 15065 | 1312 | 0.9734 | 0.2083 | 0.9274 | 7 | 9 | 7 | 0 |
| minor_isoform | min_usage=0.05 | False | 12954 | 3423 | 0.7179 | 0.3676 | 0.6721 | 11 | 9 | 8 | 0 |
| minor_isoform | min_usage=0.1 | False | 10894 | 5483 | 0.4425 | 0.442 | 0.4874 | 17 | 9 | 8 | 1 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 1120 | 0.07337 | 0.05714 | 0.2348 | gtex | substantia_nigra |
| [3, 5) | 2492 | 0.06781 | 0.04494 | 0.2063 | gtex | substantia_nigra |
| [5, 10) | 5420 | 0.06958 | 0.04723 | 0.2155 | gtex | substantia_nigra |
| [10, 20) | 4720 | 0.07631 | 0.05614 | 0.3044 | gtex | substantia_nigra |
| [20, 10000) | 2625 | 0.08893 | 0.09333 | 0.4465 | gtex | substantia_nigra |
