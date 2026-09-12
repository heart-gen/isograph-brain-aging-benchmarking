# BrainSEQ signal-level colocalization: switch (S_g) and abundance (A_g) QTLs

Arm: `ea_only` (European-ancestry donors; both BrainSEQ QTL checks passed in every region). BrainSEQ is the discovery cohort, so this is **same-tissue genetic anchoring, not replication**. GTEx is the external cohort.

GWAS side: the stage-A SuSiE cache the GTEx layer uses (1000G EUR LD). QTL side: `susie_rss` on LD computed **in-sample** from the donors each region was mapped in, as the correlation of ALT dosages residualized on that axis's covariates (the genotypes the z-scores came from), so no GTEx-style reference-LD agreement filter is applied. Cells where either trait does not fine-map fall back to `coloc.abf` and are flagged.

## What was fit

- signal-pair posteriors: 2,140
- (cell, axis) rows: 9,378; estimator by region and axis:

| region | axis | coloc.susie | coloc.abf | PP4 >= 0.8 |
|---|---|---|---|---|
| caudate | A_g | 149 | 1,551 | 23 |
| caudate | S_g | 21 | 1,401 | 2 |
| dlpfc | A_g | 114 | 1,595 | 17 |
| dlpfc | S_g | 16 | 1,414 | 2 |
| hippocampus | A_g | 53 | 1,650 | 11 |
| hippocampus | S_g | 9 | 1,405 | 1 |

Why coloc.abf, first place the cell left the signal-level pipeline: `gwas_no_credible_set` 3,594, `no_qtl_credible_set` 3,042, `gwas_locus_over_max_snps` 2,380.

QTL-side fit outcomes: `gwas_no_credible_set` 3,702, `no_qtl_credible_set` 3,042, `gwas_locus_over_max_snps` 2,380, `coloc_susie` 362, `gwas_too_few_snps` 14.

Prior robustness of the calls at the primary prior: A_g abf `intermediate` 17, A_g abf `primary_prior` 15, A_g abf `robust` 3, A_g susie `intermediate` 8, A_g susie `primary_prior` 4, A_g susie `robust` 4, S_g abf `primary_prior` 4, S_g susie `intermediate` 1.

## Paired contrast: S_g against A_g

Cells enter only when both axes cleared the shared-SNP floor; each gene contributes its best region on both axes symmetrically.

| analysis | trait | genes | S_g coloc | A_g coloc | switch-only | abundance-only | McNemar P | Wilcoxon (conditional) P |
|---|---|---|---|---|---|---|---|---|
| aging__ad | ad | 309 | 1 | 9 | 1 | 9 | 0.021 | 0.000117 |
| aging__als | als | 145 | 0 | 4 | 0 | 4 | 0.125 | 0.004 |
| aging__lbd | lbd | 53 | 0 | 3 | 0 | 3 | 0.250 | 0.207 |
| aging__pd | pd | 134 | 0 | 2 | 0 | 2 | 0.500 | 0.002 |
| aging__scz | scz | 763 | 4 | 14 | 4 | 14 | 0.031 | 5.4e-12 |
| brainseq-sczd__scz | scz | 59 | 0 | 1 | 0 | 1 | 1.000 | 0.526 |
| POOLED | ALL | 1463 | 5 | 33 | 5 | 33 | 4.26e-06 | 3.12e-19 |

## The GTEx signal-level nominations, read in BrainSEQ

Each cell is `S_g PP4 / A_g PP4`, suffixed `s` (coloc.susie) or `a` (coloc.abf); `*` marks a region whose tissue-matched GTEx tissue carries the GTEx sQTL call. `—` is untested (no BrainSEQ phenotype, or under the shared-SNP floor).

| gene | trait | caudate | dlpfc | hippocampus | S_g call in any region |
|---|---|---|---|---|---|
| ASB3 | scz | 0.111a / 0.372a | 0.053a / 0.207a | 0.049a / 0.054a | no |
| AZI2 | pd | 0.025a / 0.014a | 0.032a / 0.017a | 0.035a / 0.021a | no |
| CDIP1 | scz | 0.033a / 0.384a* | 0.046a / 0.178a | 0.040a / 0.045a* | no |
| COPA | scz | 0.025a / 0.200a | 0.121a / 0.199a | 0.082a / 0.243a | no |
| CTSH | ad | 0.042a / 0.990a | 0.036a / 0.256a | 0.072a / 0.434a* | no |
| DOC2A | ad | 0.069a / 0.005s | 0.060a / 0.136a | 0.052a / 0.164s | no |
| DOC2A | scz | 0.139a / 0.007s* | 0.052a / 0.161a* | 0.121a / 0.051s* | no |
| FAM221A | scz | 0.722s / 0.010s | 0.082a / 0.009s | 0.106a / 0.010s | no |
| G2E3 | als | 0.100a / 0.030a | 0.045a / 0.822a | 0.059a / 0.042a | no |
| GABBR2 | scz | 0.049a / 0.023a | 0.051a / 0.065a | — / 0.048a | no |
| GGNBP2 | als | 0.045a / 0.947s | 0.047a / 0.376a | 0.063a / 0.829a | no |
| GPM6A | scz | 0.071a / 0.168a* | 0.000606s / 0.026a | 0.204a / 0.027a | no |
| HMOX2 | scz | 0.036a / 0.272a | 0.045a / 0.019a* | 0.043a / 0.131a | no |
| KLC1 | scz | 0.040a / 1.13e-06s | 0.049a / 0.003s | 0.038a / 0.001s | no |
| NCOR1 | pd | 0.397a / 0.891s | 0.088a / 0.135a | 0.060a / 0.499a | no |
| NDUFS3 | ad | 0.151a / 0.000144s | 0.066a / 0.027a | 0.048a / 0.021a | no |
| NEK4 | scz | 0.058a / 0.003a | 0.211a / 0.368a* | 0.082a / 0.076a | no |
| NT5C2 | scz | 0.032a / 0.347a* | 0.085a / 0.845s* | 0.077a / 0.530a | no |
| PGS1 | als | 0.025a / 0.006a* | 0.045a / 0.026a* | 0.035a / 0.087a* | no |
| PICALM | ad | 0.114a / 0.018a | 0.048a / 0.022a | 0.032a / 0.019a | no |
| PITPNM2 | pd | 0.239a / 0.015a | 0.032a / 0.008a | 0.046a / 0.019a | no |
| PLCB2 | scz | 0.039a / 0.059a | 0.095a / 0.062a | 0.132a / 0.071a | no |
| PPIL2 | scz | 0.811a / 0.664a* | 0.297a / 0.042a* | 0.173a / 0.095a | **yes** |
| PPIP5K1 | scz | 0.059a / 0.180a | 0.085a / 0.205a | 0.181a / 0.210a | no |
| PTPRN | als | 0.074a / 0.812s | 0.035a / 0.029a* | 0.034a / 0.388a | no |
| RNASEH2C | scz | 0.079a / 0.981s* | 0.065a / 0.983s | 0.037a / 0.798a | no |
| SCFD1 | als | 0.318a / 0.946s | 0.066a / 0.936s | 0.039a / 0.931s | no |
| SH3GL2 | pd | — / 0.001a* | — / 0.964a* | — / 0.024a | no |
| SIRPA | ad | 0.024a / 0.960a* | 0.024a / 0.961a* | 0.025a / 0.946a* | no |
| SNCA | lbd | 0.041a / 0.857a | 0.209a / 0.021a* | 0.052a / 0.031a* | no |
| SNCA | pd | 0.037a / 0.012a | 0.032a / 0.020a* | 0.021a / 0.016a* | no |
| SPAG9 | ad | 0.025a / 0.063a* | 0.047a / 0.054a | 0.023a / 0.045a* | no |
| SPI1 | ad | — / 0.025a | — / 0.035a | — / 0.026a | no |
| SYT5 | scz | 0.055a / 0.018a | 0.019a / 0.004a* | 0.021a / 0.230a | no |
| TMED4 | scz | 0.080a / 0.356s* | 0.035a / 0.190s* | 0.038a / 0.285a* | no |
| TPCN1 | ad | 0.070a / 0.215a | 0.057a / 0.104a | 0.086a / 0.067a | no |
| TPP1 | als | 0.062a / 0.035a | 0.045a / 0.111a | 0.088a / 0.022a | no |
| TTC19 | pd | 0.041a / 0.737s* | 0.044a / 0.575a | 0.031a / 0.091a | no |
| TXNDC15 | als | 0.044a / 0.028a* | 0.772a / 0.288a | 0.662a / 0.317a | no |
| UNC13A | als | 0.041a / 0.009a | 0.044a / 0.027a | 0.093a / 0.066a | no |
| WIPI2 | als | 0.019a / 0.382a | 0.039a / 0.048a | 0.034a / 0.025a | no |
| ZNF232 | ad | 0.045a / 0.046a | 0.449a / 0.167a | — / 1.71e-08a | no |

## BrainSEQ's own switch-axis nominations

5 (analysis, locus, gene) cells reach PP4_S_g >= 0.8 in at least one region; 1 are also GTEx signal-level nominations. The count is a set of locus nominations selected on S_g, not a switch-specificity estimate; the paired test above is the unbiased comparison.

| gene | trait | PP4 S_g | PP4 A_g | regions coloc | best region | GTEx nomination |
|---|---|---|---|---|---|---|
| TMEM106B | ad | 0.923 | 0.039 | 1/3 | dlpfc | no |
| NGEF | scz | 0.886 | 0.368 | 1/3 | caudate | no |
| ZNF592 | scz | 0.839 | 0.125 | 1/3 | hippocampus | no |
| PPIL2 | scz | 0.811 | 0.664 | 1/3 | caudate | yes |
| NSF | scz | 0.804 | 0.670 | 1/3 | dlpfc | no |

## Reading rules

- An S_g colocalization shares a variant with a gene's switch coordinate (PC1 of its within-gene composition). It names no intron or event; the GTEx all-introns arm and the event audit do that.
- The PP4 per gene is a maximum over three regions; quote it with the region count.
- PP4_S_g >= 0.8 with low PP4_A_g is switch-preferential colocalization under prespecified thresholds, not a demonstrated switch-mediated mechanism.
- Region donors are 169-229, so a non-colocalizing A_g or S_g can be a power result.

