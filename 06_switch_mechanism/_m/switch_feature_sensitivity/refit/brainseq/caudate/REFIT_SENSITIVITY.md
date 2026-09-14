# Full refit per preprocessing setting

`switch_feature_refit.py` on **brainseq/caudate**: 12 fits of the production IsoGraph configuration, one per distinct setting of the pseudocount, expression-filter and minor-isoform axes, each compared with the committed fit.

## The noise floor

The `published` row refits the published setting itself. Its agreement with the committed partition is the ceiling for every other row; a perturbed setting that matches it has cost nothing beyond refit noise.

Published-setting refit: ARI **0.966**, NMI **0.944**, 19/20 published age-significant modules retained.

## Every setting

| axis | setting | n_modules_refit | n_genes_both | ari | nmi | n_age_sig_published | n_age_sig_refit | n_age_sig_retained | median_best_jaccard_age_sig | selected_alpha_abundance |
|---|---|---|---|---|---|---|---|---|---|---|
| published | published | 45 | 5551 | 0.966 | 0.944 | 20 | 19 | 19 | 0.695 | 0.95 |
| pseudocount | pseudocount=0.1 | 64 | 4583 | 0.686 | 0.677 | 20 | 30 | 18 | 0.155 | 0.95 |
| pseudocount | pseudocount=0.25 | 56 | 4750 | 0.766 | 0.736 | 20 | 26 | 19 | 0.266 | 0.95 |
| pseudocount | pseudocount=1 | 40 | 5066 | 0.714 | 0.72 | 20 | 20 | 19 | 0.262 | 0.95 |
| pseudocount | pseudocount=2 | 47 | 5039 | 0.659 | 0.642 | 20 | 21 | 19 | 0.16 | 0.9 |
| expression | count>5,frac>=0.5 | 26 | 2628 | 0.552 | 0.54 | 20 | 10 | 18 | 0.0543 | 0.95 |
| expression | count>10,frac>=0.5 | 35 | 3435 | 0.594 | 0.592 | 20 | 21 | 19 | 0.117 | 0.95 |
| expression | count>20,frac>=0.7 | 53 | 4480 | 0.599 | 0.6 | 20 | 28 | 19 | 0.164 | 0.95 |
| expression | count>10,frac>=0.9 | 51 | 4492 | 0.44 | 0.481 | 20 | 29 | 20 | 0.138 | 0.95 |
| minor_isoform | min_usage=0.01 | 47 | 4718 | 0.651 | 0.655 | 20 | 23 | 18 | 0.21 | 0.95 |
| minor_isoform | min_usage=0.05 | 39 | 4742 | 0.427 | 0.453 | 20 | 6 | 15 | 0.0422 | 0.95 |
| minor_isoform | min_usage=0.1 | 35 | 4633 | 0.285 | 0.336 | 20 | 7 | 15 | 0.0204 | 0.95 |

`n_age_sig_retained`: published FDR-significant module–age associations whose best-Jaccard refit counterpart is FDR-significant with the same sign. Agreement metrics are over genes assigned in both partitions.
