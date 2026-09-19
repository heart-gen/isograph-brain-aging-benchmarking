# Full refit per preprocessing setting

`switch_feature_refit.py` on **brainseq/caudate**: 18 fits of the production IsoGraph configuration, one per distinct setting of the pseudocount, expression-filter and minor-isoform axes, each compared with the committed fit.

## The noise floor

The `published` row refits the published setting itself. Its agreement with the committed partition is the ceiling for every other row; a perturbed setting that matches it has cost nothing beyond refit noise.

Published-setting refit: ARI **0.988**, NMI **0.966**, 1/1 published age-significant modules retained.

## Every setting

| axis | setting | n_modules_refit | n_genes_both | ari | nmi | n_age_sig_published | n_age_sig_refit | n_age_sig_retained | median_best_jaccard_age_sig | selected_alpha_abundance |
|---|---|---|---|---|---|---|---|---|---|---|
| published | published | 14 | 6847 | 0.988 | 0.966 | 1 | 1 | 1 | 1 | 0.95 |
| pseudocount | pseudocount=0.1 | 17 | 5844 | 0.792 | 0.714 | 1 | 1 | 1 | 0.761 | 0.95 |
| pseudocount | pseudocount=0.25 | 14 | 6282 | 0.821 | 0.746 | 1 | 1 | 1 | 0.839 | 0.95 |
| pseudocount | pseudocount=1 | 14 | 6585 | 0.821 | 0.745 | 1 | 1 | 1 | 0.836 | 0.95 |
| pseudocount | pseudocount=2 | 17 | 6566 | 0.712 | 0.673 | 1 | 2 | 1 | 0.787 | 0.95 |
| expression | count>5,frac>=0.5 | 26 | 2628 | 0.552 | 0.54 | 20 | 10 | 18 | 0.0543 | 0.95 |
| expression | switching share>=0.05,frac>=0.1 | 12 | 4859 | 0.528 | 0.448 | 1 | 1 | 1 | 0.527 | 0.95 |
| expression | switching share>=0.2,frac>=0.1 | 14 | 5624 | 0.444 | 0.395 | 1 | 5 | 1 | 0.469 | 0.95 |
| expression | count>10,frac>=0.5 | 35 | 3435 | 0.594 | 0.592 | 20 | 21 | 19 | 0.117 | 0.95 |
| expression | count>20,frac>=0.7 | 53 | 4480 | 0.599 | 0.6 | 20 | 28 | 19 | 0.164 | 0.95 |
| expression | switching share>=0.1,frac>=0.05 | 15 | 5811 | 0.76 | 0.665 | 1 | 1 | 1 | 0.707 | 0.95 |
| expression | switching share>=0.1,frac>=0.25 | 15 | 6236 | 0.63 | 0.571 | 1 | 2 | 1 | 0.6 | 0.95 |
| expression | count>10,frac>=0.9 | 51 | 4492 | 0.44 | 0.481 | 20 | 29 | 20 | 0.138 | 0.95 |
| expression | legacy count>10,frac>=0.7 | 12 | 4933 | 0.229 | 0.248 | 1 | 3 | 1 | 0.268 | 0.95 |
| expression | no filter | 17 | 2831 | 0.324 | 0.252 | 1 | 4 | 1 | 0.31 | 0.95 |
| minor_isoform | min_usage=0.05 | 39 | 4742 | 0.427 | 0.453 | 20 | 6 | 15 | 0.0422 | 0.95 |
| minor_isoform | min_usage=0.1 | 35 | 4633 | 0.285 | 0.336 | 20 | 7 | 15 | 0.0204 | 0.95 |
| minor_isoform | min_usage=0.01 | 14 | 6847 | 0.988 | 0.966 | 1 | 1 | 1 | 1 | 0.95 |

`n_age_sig_retained`: published FDR-significant module–age associations whose best-Jaccard refit counterpart is FDR-significant with the same sign. Agreement metrics are over genes assigned in both partitions.
