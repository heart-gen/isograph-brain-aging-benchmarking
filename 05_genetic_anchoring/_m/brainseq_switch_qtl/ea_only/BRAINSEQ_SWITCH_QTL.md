# BrainSEQ switch-QTL vs abundance-QTL

Does genetic variation preferentially regulate the **relative-isoform** axis (`S_g`) or the **total-abundance** axis (`A_g`) of the same gene? Same donors, same variants, same cis window, same covariates, one phenotype per gene per axis.

**This is same-tissue genetic anchoring, not independent replication.** BrainSEQ is the cohort the switches were discovered in. GTEx remains the external cohort.

## Read the raw counts with care

`A_g` is a directly measured log-CPM; `S_g` is PC1 of a within-gene CLR isoform composition and inherits isoform-quantification noise. Two phenotypes of unequal precision tested at the same alpha differ in yield for reasons unrelated to biology, so **the raw swQTL/eQTL counts below are confounded with measurement precision and are not the result.** They are reported because hiding them would be worse. The primary statistics are the precision-cancelling ones beneath.

| region | genes | swQTL | eQTL | switch-only | abundance-only | both |
|---|---|---|---|---|---|---|
| caudate | 11,133 | 966 | 4,358 | 326 | 3,718 | 640 |
| dlpfc | 11,099 | 653 | 3,610 | 272 | 3,229 | 381 |
| hippocampus | 10,605 | 523 | 2,340 | 226 | 2,043 | 297 |

## Primary statistics

`pi1` conditions on one axis being significant and asks what fraction of the OTHER axis carries signal — a ratio of ratios, with no threshold on the second axis. Effect sizes are compared only on genes where BOTH axes cleared the same bar, so detection rate largely cancels.

| region | pi1(A given S) | pi1(S given A) | n both | same lead variant | median abs slope S | median abs slope A | Wilcoxon P |
|---|---|---|---|---|---|---|---|
| caudate | 0.845 | 0.342 | 640 | 78 (0.12) | 0.631 | 0.411 | 7.18e-35 |
| dlpfc | 0.798 | 0.287 | 381 | 45 (0.12) | 0.747 | 0.542 | 1.36e-21 |
| hippocampus | 0.771 | 0.295 | 297 | 30 (0.10) | 0.685 | 0.423 | 7.06e-31 |

`yield_by_power_bin.parquet` carries the same yields within deciles of the number of variants tested, which is where most of the precision difference lives.

### Two things not to overclaim from the numbers above

**The two phenotypes are not independent.** `S_g` and `A_g` are both derived from the same transcript counts: a variant that strongly changes one isoform's abundance moves the gene total AND the within-gene composition. So some sharing is mechanically expected and `pi1(A given S)` should not be read as evidence that a switch signal *causes* an abundance signal. The informative comparison is the ASYMMETRY between the two pi1 values, and the same-lead-variant fraction.

**Effect sizes are comparable only because the transform was identical.** Both axes are rank-based inverse-normal transformed before mapping, so slopes are in SD units of the transformed phenotype and can be compared. That is a property of this pipeline, not of QTL slopes in general — comparing an untransformed composition slope against an untransformed log-CPM slope would be meaningless.

FDR: Benjamini-Hochberg at q < 0.05 on the permutation p-values, applied identically to both axes. tensorQTL's Storey q-values need rpy2 and the R `qvalue` package, which do not import in this environment.

