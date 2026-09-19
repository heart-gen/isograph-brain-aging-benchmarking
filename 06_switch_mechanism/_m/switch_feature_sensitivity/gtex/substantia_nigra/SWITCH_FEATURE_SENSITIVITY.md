# Preprocessing sensitivity of the switch representation

`switch_feature_sensitivity.py` on **gtex/substantia_nigra** (RSEM, n=183, 24 modules). The cross-region summary and the quantification axis are in `../../SWITCH_FEATURE_SENSITIVITY.md`.

## Scope, stated up front

module partition held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. A full refit per setting would additionally let the network change and is not done here.

Published settings for this cohort: `{"pseudocount": 0.5, "transcript_filter": "switching", "min_tx_prop": 0.1, "min_tx_fraction": 0.1, "min_usage": 0.0}`. Gate passed: the baseline rebuild reproduces the published switch channel with max |diff| = 0 (tolerance 1e-10), so every comparison below is against the published quantity and not a lookalike.

## 1-3. Pseudocount, expression filter, minor-isoform threshold

| axis | setting | is_published | n_switch_genes | n_switch_genes_lost | median_abs_feature_r_vs_published | frac_features_sign_flipped | effect_pearson_vs_published | n_fdr_sig | n_fdr_sig_published | n_published_sig_retained | n_sign_flips_among_published_sig |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pseudocount | pseudocount=0.1 | False | 12286 | 0 | 0.9947 | 0.02849 | 0.9945 | 2 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.25 | False | 12286 | 0 | 0.9988 | 0.01392 | 0.998 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=0.5 | True | 12286 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=1 | False | 12286 | 0 | 0.9984 | 0.02076 | 0.9973 | 0 | 0 | 0 | 0 |
| pseudocount | pseudocount=2 | False | 12286 | 0 | 0.9924 | 0.04371 | 0.9876 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.1 | True | 12286 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.05,frac>=0.1 | False | 13286 | 0 | 1 | 0.1463 | 0.9699 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.2,frac>=0.1 | False | 10307 | 1979 | 1 | 0.1974 | 0.8369 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.05 | False | 12879 | 0 | 1 | 0.08668 | 0.9849 | 0 | 0 | 0 | 0 |
| expression | switching share>=0.1,frac>=0.25 | False | 10869 | 1417 | 1 | 0.163 | 0.918 | 0 | 0 | 0 | 0 |
| expression | legacy count>10,frac>=0.7 | False | 10184 | 3099 | 0.9058 | 0.2942 | 0.7057 | 0 | 0 | 0 | 0 |
| expression | no filter | False | 16377 | 0 | 0.7896 | 0.3502 | 0.6375 | 7 | 0 | 0 | 0 |
| minor_isoform | min_usage=0 | True | 12286 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.01 | False | 12286 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.05 | False | 12118 | 168 | 1 | 0.04151 | 0.9946 | 0 | 0 | 0 | 0 |
| minor_isoform | min_usage=0.1 | False | 10816 | 1470 | 1 | 0.2025 | 0.8081 | 0 | 0 | 0 | 0 |

`median_abs_feature_r_vs_published` is the per-gene correlation of the rebuilt switch coordinate with the published one; `n_published_sig_retained` is how many of the published FDR-significant module-age associations survive. A sign flip in a switch coordinate is not itself a problem -- PC1's sign is arbitrary and sign-stabilised -- but a flip among *published-significant modules* would change the direction of a reported effect, so it is counted separately.

## 4. Identifiability (transcript number)

If the switch signal were a quantification artefact it should grow with the number of annotated isoforms, which is what makes a gene hard to quantify. It is reported by stratum rather than adjusted away:

| n_tx_stratum | n_genes | median_abs_age_r | frac_abs_age_r_above_0_2 | module_membership_rate | cohort | region |
|---|---|---|---|---|---|---|
| [2, 3) | 4740 | 0.07221 | 0.03945 | 0.5738 | gtex | substantia_nigra |
| [3, 5) | 5386 | 0.06849 | 0.03231 | 0.4987 | gtex | substantia_nigra |
| [5, 10) | 2128 | 0.06738 | 0.03618 | 0.4272 | gtex | substantia_nigra |
| [10, 20) | 32 | 0.07564 | 0.0625 | 0.25 | gtex | substantia_nigra |
