# Full refit per preprocessing setting

`switch_feature_refit.py` on **brainseq/caudate**: 18 fits of the production IsoGraph configuration, one per distinct setting of the pseudocount, expression-filter and minor-isoform axes, each compared with the committed fit.

## The noise floor

The `published` row refits the published setting itself. Its agreement with the committed partition is the ceiling for every other row; a perturbed setting that matches it has cost nothing beyond refit noise.

Published-setting refit: ARI **0.985**, NMI **0.975**, 4/4 published age-significant modules retained.

## Every setting

| axis | setting | n_modules_refit | n_genes_both | ari | nmi | n_age_sig_published | n_age_sig_refit | n_age_sig_retained | median_best_jaccard_age_sig | selected_alpha_abundance |
|---|---|---|---|---|---|---|---|---|---|---|
| published | published | 29 | 5665 | 0.985 | 0.975 | 4 | 4 | 4 | 0.994 | 0.95 |
| pseudocount | pseudocount=0.1 | 30 | 4729 | 0.721 | 0.719 | 4 | 3 | 3 | 0.452 | 0.95 |
| pseudocount | pseudocount=0.25 | 33 | 5060 | 0.746 | 0.747 | 4 | 4 | 3 | 0.516 | 0.95 |
| pseudocount | pseudocount=1 | 28 | 5284 | 0.775 | 0.763 | 4 | 4 | 3 | 0.529 | 0.95 |
| pseudocount | pseudocount=2 | 25 | 5120 | 0.656 | 0.698 | 4 | 5 | 3 | 0.16 | 0.95 |
| expression | count>5,frac>=0.5 | 26 | 2628 | 0.552 | 0.54 | 20 | 10 | 18 | 0.0543 | 0.95 |
| expression | switching share>=0.05,frac>=0.1 | 34 | 3825 | 0.493 | 0.514 | 4 | 8 | 3 | 0.163 | 0.95 |
| expression | switching share>=0.2,frac>=0.1 | 35 | 4533 | 0.476 | 0.455 | 4 | 11 | 3 | 0.115 | 0.95 |
| expression | count>10,frac>=0.5 | 35 | 3435 | 0.594 | 0.592 | 20 | 21 | 19 | 0.117 | 0.95 |
| expression | count>20,frac>=0.7 | 53 | 4480 | 0.599 | 0.6 | 20 | 28 | 19 | 0.164 | 0.95 |
| expression | switching share>=0.1,frac>=0.05 | 28 | 4614 | 0.717 | 0.696 | 4 | 4 | 4 | 0.457 | 0.95 |
| expression | switching share>=0.1,frac>=0.25 | 32 | 5051 | 0.57 | 0.597 | 4 | 7 | 3 | 0.184 | 0.95 |
| expression | count>10,frac>=0.9 | 51 | 4492 | 0.44 | 0.481 | 20 | 29 | 20 | 0.138 | 0.95 |
| expression | legacy count>10,frac>=0.7 | 45 | 3563 | 0.337 | 0.379 | 4 | 19 | 4 | 0.059 | 0.95 |
| expression | no filter | 25 | 2017 | 0.347 | 0.357 | 4 | 7 | 3 | 0.0472 | 0.95 |
| minor_isoform | min_usage=0.05 | 39 | 4742 | 0.427 | 0.453 | 20 | 6 | 15 | 0.0422 | 0.95 |
| minor_isoform | min_usage=0.1 | 35 | 4633 | 0.285 | 0.336 | 20 | 7 | 15 | 0.0204 | 0.95 |
| minor_isoform | min_usage=0.01 | 29 | 5665 | 0.985 | 0.975 | 4 | 4 | 4 | 0.994 | 0.95 |

`n_age_sig_retained`: published FDR-significant module–age associations whose best-Jaccard refit counterpart is FDR-significant with the same sign. Agreement metrics are over genes assigned in both partitions.
