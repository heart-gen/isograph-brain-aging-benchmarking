# BrainSEQ switch-QTL vs abundance-QTL

Does genetic variation preferentially regulate the **relative-isoform** axis (`S_g`) or the **total-abundance** axis (`A_g`) of the same gene? Same donors, same variants, same cis window, same covariates, one phenotype per gene per axis.

**This is same-tissue genetic anchoring, not independent replication.** BrainSEQ is the cohort the switches were discovered in. GTEx remains the external cohort.

## Read the raw counts with care

`A_g` is a directly measured log-CPM; `S_g` is PC1 of a within-gene CLR isoform composition and inherits isoform-quantification noise. Two phenotypes of unequal precision tested at the same alpha differ in yield for reasons unrelated to biology, so **the raw swQTL/eQTL counts below are confounded with measurement precision and are not the result.** They are reported because hiding them would be worse. The primary statistics are the precision-cancelling ones beneath.

| region | genes | swQTL | eQTL | switch-only | abundance-only | both |
|---|---|---|---|---|---|---|
| caudate | 12,688 | 1,155 | 5,015 | 427 | 4,287 | 728 |
| dlpfc | 12,533 | 815 | 4,098 | 313 | 3,596 | 502 |
| hippocampus | 12,488 | 567 | 2,898 | 242 | 2,573 | 325 |

## Primary statistics

`pi1` conditions on one axis being significant and asks what fraction of the OTHER axis carries signal — a ratio of ratios, with no threshold on the second axis. Effect sizes are compared only on genes where BOTH axes cleared the same bar, so detection rate largely cancels.

| region | pi1(A given S) | pi1(S given A) | n both | same lead variant | median abs slope S | median abs slope A | Wilcoxon P |
|---|---|---|---|---|---|---|---|
| caudate | 0.815 | 0.319 | 728 | 94 (0.13) | 0.619 | 0.441 | 1.45e-25 |
| dlpfc | 0.799 | 0.255 | 502 | 63 (0.13) | 0.721 | 0.560 | 1.35e-21 |
| hippocampus | 0.750 | 0.265 | 325 | 34 (0.10) | 0.704 | 0.395 | 3.11e-36 |

`yield_by_power_bin.parquet` carries the same yields within deciles of the number of variants tested, which is where most of the precision difference lives.

### Two things not to overclaim from the numbers above

**The two phenotypes are not independent.** `S_g` and `A_g` are both derived from the same transcript counts: a variant that strongly changes one isoform's abundance moves the gene total AND the within-gene composition. So some sharing is mechanically expected and `pi1(A given S)` should not be read as evidence that a switch signal *causes* an abundance signal. The informative comparison is the ASYMMETRY between the two pi1 values, and the same-lead-variant fraction.

**Effect sizes are comparable only because the transform was identical.** Both axes are rank-based inverse-normal transformed before mapping, so slopes are in SD units of the transformed phenotype and can be compared. That is a property of this pipeline, not of QTL slopes in general — comparing an untransformed composition slope against an untransformed log-CPM slope would be meaningless.

FDR: Benjamini-Hochberg at q < 0.05 on the permutation p-values, applied identically to both axes. tensorQTL's Storey q-values need rpy2 and the R `qvalue` package, which do not import in this environment.

