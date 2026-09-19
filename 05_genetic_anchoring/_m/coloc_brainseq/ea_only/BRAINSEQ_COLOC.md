# BrainSEQ signal-level colocalization: switch (S_g) and abundance (A_g) QTLs

Arm: `ea_only` (European-ancestry donors; both BrainSEQ QTL checks passed in every region). BrainSEQ is the discovery cohort, so this is **same-tissue genetic anchoring, not replication**. GTEx is the external cohort.

GWAS side: the stage-A SuSiE cache the GTEx layer uses (1000G EUR LD). QTL side: `susie_rss` on LD computed **in-sample** from the donors each region was mapped in, as the correlation of ALT dosages residualized on that axis's covariates (the genotypes the z-scores came from), so no GTEx-style reference-LD agreement filter is applied. Cells where either trait does not fine-map fall back to `coloc.abf` and are flagged.

## What was fit

- signal-pair posteriors: 7,510
- (cell, axis) rows: 35,173; estimator by region and axis:

| region | axis | coloc.susie | coloc.abf | PP4 >= 0.8 |
|---|---|---|---|---|
| caudate | A_g | 433 | 6,066 | 70 |
| caudate | S_g | 122 | 5,103 | 13 |
| dlpfc | A_g | 347 | 6,157 | 51 |
| dlpfc | S_g | 89 | 5,066 | 10 |
| hippocampus | A_g | 228 | 6,258 | 37 |
| hippocampus | S_g | 59 | 5,245 | 7 |

Why coloc.abf, first place the cell left the signal-level pipeline: `gwas_locus_over_max_snps` 15,118, `gwas_no_credible_set` 10,663, `no_qtl_credible_set` 8,114.

QTL-side fit outcomes: `gwas_locus_over_max_snps` 15,118, `gwas_no_credible_set` 11,380, `no_qtl_credible_set` 8,059, `coloc_susie` 1,278, `too_few_shared_snps` 65, `gwas_too_few_snps` 33.

Prior robustness of the calls at the primary prior: A_g abf `intermediate` 56, A_g abf `primary_prior` 49, A_g abf `robust` 20, A_g susie `intermediate` 21, A_g susie `primary_prior` 8, A_g susie `robust` 4, S_g abf `intermediate` 8, S_g abf `primary_prior` 11, S_g abf `robust` 3, S_g susie `intermediate` 4, S_g susie `primary_prior` 3, S_g susie `robust` 1.

## Paired contrast: S_g against A_g

Cells enter only when both axes cleared the shared-SNP floor; each gene contributes its best region on both axes symmetrically.

| analysis | trait | genes | S_g coloc | A_g coloc | switch-only | abundance-only | McNemar P | Wilcoxon (conditional) P |
|---|---|---|---|---|---|---|---|---|
| aging__ad | ad | 1099 | 2 | 21 | 2 | 21 | 6.6e-05 | 3.09e-15 |
| aging__als | als | 505 | 3 | 7 | 1 | 5 | 0.219 | 2.2e-08 |
| aging__lbd | lbd | 182 | 0 | 4 | 0 | 4 | 0.125 | 0.003 |
| aging__pd | pd | 467 | 1 | 7 | 1 | 7 | 0.070 | 4.48e-09 |
| aging__scz | scz | 2691 | 13 | 48 | 9 | 44 | 1.22e-06 | 1.14e-40 |
| brainseq-sczd__scz | scz | 533 | 3 | 13 | 3 | 13 | 0.021 | 2.21e-09 |
| POOLED | ALL | 5477 | 22 | 100 | 16 | 94 | 1.29e-14 | 1.36e-78 |

## The GTEx signal-level nominations, read in BrainSEQ

Each cell is `S_g PP4 / A_g PP4`, suffixed `s` (coloc.susie) or `a` (coloc.abf); `*` marks a region whose tissue-matched GTEx tissue carries the GTEx sQTL call. `—` is untested (no BrainSEQ phenotype, or under the shared-SNP floor).

| gene | trait | caudate | dlpfc | hippocampus | S_g call in any region |
|---|---|---|---|---|---|
| ACTR1B | scz | — / 4.57e-05a* | — / 0.902a* | 0.858a / 0.569a* | **yes** |
| AKT1 | ad | 0.006a / 0.007a* | 0.007a / 0.005a | 0.006a / 0.005a | no |
| ASB3 | scz | 0.050a / 0.372a | 0.057a / 0.207a | 0.046a / 0.054a | no |
| AZI2 | pd | 0.024a / 0.014a | 0.055a / 0.017a | 0.081a / 0.021a | no |
| BCKDK | ad | 0.048a / 0.035a | 0.061a / 0.036a* | 0.157a / 0.040a | no |
| C9orf72 | als | 0.893a / 6.86e-05a* | 0.983a / 4.35e-13a* | 0.139a / 0.014a* | **yes** |
| CATSPER2 | scz | 0.061a / 0.020a | 0.091a / 0.041a | 0.062a / 0.062a* | no |
| CCDC122 | scz | 0.043a / 0.956a | 0.044a / 0.161a | 0.068a / 0.953a* | no |
| CCDC62 | pd | 0.055a / 0.050a | 0.049a / 0.041a | 0.059a / 0.032a | no |
| CCS | scz | 0.037a / 0.011a | 0.056a / 0.037a | 0.054a / 0.028a | no |
| CD46 | scz | 0.095a / 0.846a | 0.088a / 0.842a | 0.777a / 0.852a | no |
| CDHR3 | pd | 0.016a / 0.008a | 0.027a / 0.009a | 0.052a / 0.018a | no |
| CDIP1 | scz | 0.683a / 0.384a* | 0.122a / 0.178a | 0.042a / 0.045a* | no |
| COG7 | ad | 0.220a / 0.173a | 0.059a / 0.208a | 0.053a / 0.026a | no |
| COPA | scz | 0.033a / 0.200a | 0.034a / 0.199a | 0.023a / 0.243a | no |
| CRELD2 | scz | 0.924s / 0.399a* | 0.917s / 0.061a* | 0.815s / 0.702a* | **yes** |
| CTSB | pd | 0.034a / 0.891s | 0.218a / 0.870s | 0.040a / 0.672s | no |
| DDRGK1 | pd | 0.015a / 0.081a | 0.220a / 0.985a | 0.039a / 0.546a | no |
| DGKZ | scz | 0.121a / 0.078a* | 0.032a / 0.020a* | 0.052a / 0.019a* | no |
| DNAJA3 | scz | 0.045a / 0.045a | 0.041a / 0.191a | 0.043a / 0.520a | no |
| DOC2A | ad | 0.049a / 0.005s | 0.055a / 0.136a | 0.055a / 0.167s | no |
| DOC2A | scz | 0.065a / 0.007s* | 0.053a / 0.161a* | 0.060a / 0.051s* | no |
| EFHB | scz | 0.001a / 0.001a | 0.507a / 0.001a* | 0.044a / 0.001a | no |
| FAM120AOS | scz | 0.448a / 0.586a | 0.123a / 0.463a | 0.048a / 0.938a | no |
| FAM184A | scz | 0.029a / 0.099a | 0.032a / 0.012a | 0.037a / 0.017a | no |
| FANCI | scz | 0.051a / 0.099a | 0.080a / 0.039a | 0.054a / 0.038a | no |
| FGFR1 | scz | 0.047a / 0.302a* | 0.057a / 0.037a | 0.080a / 0.103a | no |
| FNBP1 | als | 0.046a / 0.018a | 0.045a / 0.021a | 0.080a / 0.023a* | no |
| FOXN2 | scz | 0.978a / 0.968a | 0.894a / 0.920a | 0.887a / 0.916a | **yes** |
| G2E3 | als | 0.034a / 0.030a | 0.037a / 0.822a | 0.031a / 0.042a | no |
| GABBR2 | scz | 0.052a / 0.023a | — / 0.065a | — / 0.048a | no |
| GALNT15 | scz | 0.030a / 0.006a | — / 0.006a | — / 0.022a | no |
| GGNBP2 | als | — / 0.947s | — / 0.376a | — / 0.829a | no |
| GLYCTK | scz | 0.068a / 0.423a | 0.133a / 0.035a | 0.189a / 0.647a | no |
| GPM6A | scz | 0.101a / 0.168a* | 0.030a / 0.026a | 0.074a / 0.027a | no |
| GPR135 | scz | 0.059a / 0.107a | 0.073a / 0.080a | 0.048a / 0.051a | no |
| HMOX2 | scz | 0.043a / 0.272a | 0.034a / 0.019a* | 0.053a / 0.131a | no |
| IDH3B | scz | 0.037a / 0.729a | 0.052a / 0.538a | 0.046a / 0.179a* | no |
| IFNAR2 | ad | 0.024a / 0.700a* | 0.035a / 0.896a | 0.032a / 0.873a* | no |
| IKBIP | scz | 0.041a / 1.19e-05s | 0.064a / 1.16e-05s* | 0.045a / 1.36e-05s | no |
| INO80E | ad | 0.263a / 0.733s | 0.445a / 0.604s | 0.089a / 0.666s | no |
| INO80E | scz | 0.289a / 0.970s | 0.569a / 0.974s | 0.163a / 0.970s* | no |
| INTS8 | ad | 0.149a / 0.064a | 0.059a / 0.035a | 0.084a / 0.016a* | no |
| IRF3 | scz | 0.114a / 0.055a* | 0.085a / 0.576a* | 0.049a / 0.172a* | no |
| ITGB1BP1 | ad | 0.055a / 0.018a | 0.245a / 0.264a* | 0.115a / 0.575a | no |
| KLC1 | scz | 0.506a / 4.46e-08a | 0.029a / 7.96e-05a | 0.590a / 0.00022a | no |
| L3HYPDH | scz | 0.095a / 0.001s | 0.030a / 0.003s | 0.083a / 0.024s | no |
| LPCAT4 | scz | 0.021a / 9.21e-06a | 0.027a / 0.033a | 0.039a / 0.014a | no |
| MAD1L1 | scz | 0.037a / 0.030a | 0.053a / 0.074a* | 0.082a / 0.054a | no |
| MAP2K5 | scz | 0.039a / 0.007a | 0.034a / 0.007a | 0.042a / 0.044a | no |
| MAP7D1 | scz | 0.100a / 0.231a | 0.045a / 0.127a | 0.066a / 0.025a* | no |
| MED19 | scz | — / 0.038a | — / 0.038a | 0.051a / 0.017a* | no |
| MRPS33 | scz | 0.046a / 0.058a | 0.054a / 0.088a | 0.047a / 0.028a | no |
| NCOR1 | pd | 0.294a / 0.892s | 0.284a / 0.135a | 0.063a / 0.499a | no |
| NDUFAF7 | scz | 0.054a / 0.950a | 0.066a / 0.340a | 0.048a / 0.188a | no |
| NDUFS3 | ad | 0.036a / 0.000192a | 0.130a / 0.027a | 0.047a / 0.021a | no |
| NEK4 | scz | 0.069a / 0.003a | 0.119a / 0.368a* | 0.112a / 0.076a | no |
| NME4 | als | 0.133a / 0.028a | 0.135a / 0.028a | 0.045a / 0.030a | no |
| NMRAL1 | scz | 0.088a / 0.000719s | 0.386a / 0.000749s* | 0.043a / 0.687a* | no |
| NSMAF | als | 0.019a / 0.021a | 0.023a / 0.441a | 0.023a / 0.374a | no |
| NT5C2 | scz | 0.003s / 0.347a* | 0.104a / 0.843s* | 0.079a / 0.530a | no |
| NUCB2 | scz | 0.044a / 0.000718a | 0.044a / 0.034a | 0.051a / 0.121a | no |
| NUP50 | scz | — / 0.687a* | — / 0.768a* | 0.035a / 0.372a* | no |
| PAK6 | scz | 0.050a / 0.935a | 0.059a / 0.0006a | 0.043a / 0.021a | no |
| PAM16 | scz | 0.025a / 0.772a | 0.043a / 0.056a | 0.176a / 0.039a | no |
| PBRM1 | scz | 0.044a / 0.043a | 0.109a / 0.086a | 0.105a / 0.023a | no |
| PCBP3 | scz | 0.073a / 0.043a | 0.040a / 0.044a | 0.199a / 0.027a | no |
| PICALM | ad | 0.024a / 0.018a | 0.062a / 0.022a | 0.066a / 0.019a | no |
| PILRB | ad | 0.004s / 2.78e-14s | 2e-12s / 2.85e-14s | 3.11e-09s / 3.22e-14s | no |
| POLG | scz | 0.057a / 0.146a | 0.043a / 0.081a | 0.057a / 0.038a | no |
| PPIL2 | scz | 0.858a / 0.664a* | 0.048a / 0.042a* | 0.027a / 0.095a | **yes** |
| PPIP5K1 | scz | 0.131a / 0.180a | 0.079a / 0.205a | 0.090a / 0.210a | no |
| PRDM2 | als | 0.030a / 0.027a | 0.032a / 0.003a | 0.050a / 0.058a* | no |
| PRMT7 | scz | 0.215a / 0.956a | 0.060a / 0.985a | 0.055a / 0.976a | no |
| PSMD6 | scz | 0.055a / 0.650a | 0.041a / 0.486a | 0.039a / 0.027a* | no |
| PTPRN | als | 0.044a / 0.820s | 0.038a / 0.029a* | 0.075a / 0.388a | no |
| RAD51C | ad | 0.095a / 0.015a | 0.047a / 0.000284s | 0.143a / 0.000281s | no |
| RAI1 | scz | — / 0.012a | — / 0.024a | — / 0.019a | no |
| RBM6 | scz | 0.274a / 0.813a | 0.069a / 0.837a | 0.077a / 0.810a | no |
| RCBTB1 | scz | 0.049a / 0.977s | 0.033a / 0.975s | 0.031a / 0.973s* | no |
| RCSD1 | als | 0.046a / 0.007a | 0.029a / 0.010a | 0.029a / 0.022a | no |
| REEP2 | scz | 0.045a / 0.030a | 0.072a / 0.020a | 0.106a / 0.158a | no |
| RERE | scz | 0.031a / 0.611a | 0.178a / 0.510a | 0.238a / 0.037a | no |
| RPS6KL1 | als | 0.236a / 0.031a | 0.496a / 0.340a | 0.043a / 0.021a | no |
| SCFD1 | als | 0.627a / 0.946s | 0.053a / 0.937s | 0.041a / 0.931s | no |
| SERPINB1 | ad | — / 0.012a | 0.152a / 0.300a | 0.039a / 0.055a | no |
| SETD6 | scz | 0.206a / 0.057a | 0.193a / 0.945s | 0.070a / 0.707a | no |
| SH3GL2 | pd | — / 0.001a* | — / 0.964a* | — / 0.024a | no |
| SIRPA | ad | 0.026a / 0.960a* | — / 0.961a* | — / 0.946a* | no |
| SLC39A13 | ad | 0.218a / 0.036a | 0.146a / 0.032a | 0.252a / 0.195a* | no |
| SNAP91 | scz | 0.035a / 0.963s* | 0.034a / 0.966s | 0.032a / 0.959a | no |
| SNCA | lbd | 0.034a / 0.857a | 0.036a / 0.021a* | 0.032a / 0.031a* | no |
| SPAG9 | ad | 0.096a / 0.063a* | 0.015a / 0.054a | 0.029a / 0.045a* | no |
| SPI1 | ad | 0.044a / 0.025a | 0.052a / 0.035a | 0.135a / 0.026a | no |
| SYT5 | scz | 0.027a / 0.018a | 0.038a / 0.004a* | 0.020a / 0.230a | no |
| TAOK2 | scz | 0.034a / 0.933a | 0.005s / 0.024a | 0.094a / 0.064a | no |
| TEAD4 | scz | 0.058a / 0.029a | 0.042a / 0.056a | — | no |
| THAP3 | scz | 0.243a / 0.006s | 0.051a / 0.051a | 0.046a / 0.032a | no |
| TMED4 | scz | 0.025a / 0.228a* | 0.037a / 0.178a* | 0.033a / 0.285a* | no |
| TMEM175 | als | 0.057a / 0.028a | 0.054a / 0.055a | 0.051a / 0.051a | no |
| TNFSF13 | als | 0.049a / 0.873a | 0.061a / 0.160a | 0.033a / 0.162a | no |
| TPCN1 | ad | 0.070a / 0.215a | 0.057a / 0.104a | 0.066a / 0.067a | no |
| TPP1 | als | 0.040a / 0.035a | — / 0.111a | — / 0.022a | no |
| TSPAN31 | scz | 0.062a / 0.047a | 0.052a / 0.081a* | 0.059a / 0.039a | no |
| TTC19 | pd | 0.020a / 0.732s* | 0.047a / 0.575a | 0.039a / 0.091a | no |
| TUBGCP4 | scz | 0.044a / 0.689a | 0.039a / 0.041a | 0.055a / 0.132a | no |
| TXNDC15 | als | — / 0.028a* | — / 0.288a | 0.123a / 0.317a | no |
| UNC13A | als | — / 0.009a | — / 0.027a | — / 0.066a | no |
| VWA5B2 | ad | 0.034a / 0.078a | 0.048a / 0.078a | 0.061a / 0.093a | no |
| WHAMM | als | — / 0.037a | — / 0.044a | — / 0.174a* | no |
| WIPI2 | als | 0.088a / 0.382a | 0.058a / 0.048a | 0.036a / 0.025a | no |
| YPEL1 | scz | — / 0.095a* | — / 0.601a* | — / 0.619a | no |
| YPEL3 | ad | 0.161a / 0.028a* | 0.089a / 0.637a* | 0.043a / 0.358a | no |
| YPEL3 | scz | 0.077a / 0.037a | 0.067a / 0.053a | 0.045a / 0.095a | no |
| YWHAB | scz | 0.837a / 0.845a* | — / 0.856a* | — / 0.660a* | **yes** |
| ZDHHC12 | scz | — / 0.097a* | — / 0.059a* | — / 0.067a* | no |
| ZFYVE21 | scz | 0.057a / 0.024a* | 0.024a / 0.780a | 0.028a / 0.023a | no |
| ZNF232 | ad | 0.041a / 0.046a | 0.052a / 0.167a | 0.069a / 1.71e-08a | no |
| ZSWIM7 | pd | 0.284a / 0.775s* | 0.053a / 0.767s* | 0.078a / 0.873s* | no |

## BrainSEQ's own switch-axis nominations

22 (analysis, locus, gene) cells reach PP4_S_g >= 0.8 in at least one region; 7 are also GTEx signal-level nominations. The count is a set of locus nominations selected on S_g, not a switch-specificity estimate; the paired test above is the unbiased comparison.

| gene | trait | PP4 S_g | PP4 A_g | regions coloc | best region | GTEx nomination |
|---|---|---|---|---|---|---|
| SARM1 | als | 0.999 | 0.979 | 2/3 | hippocampus | no |
| CNOT7 | scz | 0.993 | 0.948 | 3/3 | dlpfc | no |
| C9orf72 | als | 0.983 | 0.014 | 2/3 | dlpfc | yes |
| FOXN2 | scz | 0.978 | 0.968 | 3/3 | caudate | yes |
| STX4 | pd | 0.942 | 0.332 | 1/3 | caudate | no |
| MYO19 | als | 0.930 | 0.953 | 1/3 | dlpfc | no |
| CRELD2 | scz | 0.924 | 0.702 | 3/3 | caudate | yes |
| ATE1 | scz | 0.912 | 0.716 | 1/3 | caudate | no |
| ATE1 | scz | 0.912 | 0.716 | 1/3 | caudate | no |
| PPDPF | scz | 0.890 | 0.051 | 1/3 | caudate | no |
| ZNF592 | scz | 0.888 | 0.125 | 1/3 | hippocampus | no |
| PLXNB2 | scz | 0.885 | 0.034 | 1/3 | dlpfc | no |
| MYO19 | scz | 0.880 | 0.938 | 1/3 | dlpfc | no |
| SCAMP4 | ad | 0.878 | 0.025 | 1/3 | caudate | no |
| ACTR1B | scz | 0.858 | 0.569 | 1/1 | hippocampus | yes |
| ACTR1B | scz | 0.858 | 0.569 | 1/1 | hippocampus | yes |
| PPIL2 | scz | 0.858 | 0.664 | 1/3 | caudate | yes |
| TMEM106B | ad | 0.846 | 0.039 | 1/3 | dlpfc | no |
| YWHAB | scz | 0.837 | 0.845 | 1/1 | caudate | yes |
| ARL14EP | scz | 0.834 | 0.794 | 1/3 | caudate | no |
| ARL14EP | scz | 0.834 | 0.794 | 1/3 | caudate | no |
| VPS29 | scz | 0.822 | 0.198 | 1/3 | dlpfc | no |

## Reading rules

- An S_g colocalization shares a variant with a gene's switch coordinate (PC1 of its within-gene composition). It names no intron or event; the GTEx all-introns arm and the event audit do that.
- The PP4 per gene is a maximum over three regions; quote it with the region count.
- PP4_S_g >= 0.8 with low PP4_A_g is switch-preferential colocalization under prespecified thresholds, not a demonstrated switch-mediated mechanism.
- Region donors are 169-229, so a non-colocalizing A_g or S_g can be a power result.

