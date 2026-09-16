# BrainSEQ signal-level colocalization: switch (S_g) and abundance (A_g) QTLs

Arm: `ea_only` (European-ancestry donors; both BrainSEQ QTL checks passed in every region). BrainSEQ is the discovery cohort, so this is **same-tissue genetic anchoring, not replication**. GTEx is the external cohort.

GWAS side: the stage-A SuSiE cache the GTEx layer uses (1000G EUR LD). QTL side: `susie_rss` on LD computed **in-sample** from the donors each region was mapped in, as the correlation of ALT dosages residualized on that axis's covariates (the genotypes the z-scores came from), so no GTEx-style reference-LD agreement filter is applied. Cells where either trait does not fine-map fall back to `coloc.abf` and are flagged.

## What was fit

- signal-pair posteriors: 3,005
- (cell, axis) rows: 10,749; estimator by region and axis:

| region | axis | coloc.susie | coloc.abf | PP4 >= 0.8 |
|---|---|---|---|---|
| caudate | A_g | 173 | 1,894 | 26 |
| caudate | S_g | 49 | 1,459 | 4 |
| dlpfc | A_g | 150 | 1,932 | 15 |
| dlpfc | S_g | 21 | 1,448 | 0 |
| hippocampus | A_g | 70 | 2,012 | 9 |
| hippocampus | S_g | 21 | 1,520 | 3 |

Why coloc.abf, first place the cell left the signal-level pipeline: `gwas_no_credible_set` 4,142, `no_qtl_credible_set` 3,482, `gwas_locus_over_max_snps` 2,641.

QTL-side fit outcomes: `gwas_no_credible_set` 4,290, `no_qtl_credible_set` 3,482, `gwas_locus_over_max_snps` 2,641, `coloc_susie` 484, `gwas_too_few_snps` 18.

Prior robustness of the calls at the primary prior: A_g abf `intermediate` 21, A_g abf `primary_prior` 9, A_g abf `robust` 6, A_g susie `intermediate` 10, A_g susie `primary_prior` 4, S_g abf `primary_prior` 6, S_g susie `intermediate` 1.

## Paired contrast: S_g against A_g

Cells enter only when both axes cleared the shared-SNP floor; each gene contributes its best region on both axes symmetrically.

| analysis | trait | genes | S_g coloc | A_g coloc | switch-only | abundance-only | McNemar P | Wilcoxon (conditional) P |
|---|---|---|---|---|---|---|---|---|
| aging__ad | ad | 323 | 0 | 9 | 0 | 9 | 0.004 | 0.000473 |
| aging__als | als | 136 | 0 | 1 | 0 | 1 | 1.000 | 8.66e-07 |
| aging__lbd | lbd | 57 | 0 | 1 | 0 | 1 | 1.000 | 0.008 |
| aging__pd | pd | 129 | 1 | 4 | 1 | 4 | 0.375 | 0.000117 |
| aging__scz | scz | 809 | 4 | 13 | 3 | 12 | 0.035 | 6.3e-17 |
| brainseq-sczd__scz | scz | 141 | 2 | 3 | 2 | 3 | 1.000 | 0.000816 |
| POOLED | ALL | 1595 | 7 | 31 | 6 | 30 | 6.96e-05 | 6.74e-31 |

## The GTEx signal-level nominations, read in BrainSEQ

Each cell is `S_g PP4 / A_g PP4`, suffixed `s` (coloc.susie) or `a` (coloc.abf); `*` marks a region whose tissue-matched GTEx tissue carries the GTEx sQTL call. `—` is untested (no BrainSEQ phenotype, or under the shared-SNP floor).

| gene | trait | caudate | dlpfc | hippocampus | S_g call in any region |
|---|---|---|---|---|---|
| ACTR1B | scz | — / 4.57e-05a* | — / 0.902a* | 0.858a / 0.569a* | **yes** |
| CDIP1 | scz | 0.683a / 0.384a* | 0.122a / 0.178a | 0.042a / 0.045a* | no |
| COPA | scz | 0.033a / 0.200a | 0.034a / 0.199a | 0.023a / 0.243a | no |
| CTSB | pd | 0.034a / 0.891s | 0.218a / 0.870s | 0.040a / 0.672s | no |
| DDRGK1 | pd | 0.015a / 0.081a | 0.220a / 0.985a | 0.039a / 0.546a | no |
| DGKZ | scz | 0.121a / 0.078a* | 0.032a / 0.020a* | 0.052a / 0.019a* | no |
| FGFR1 | scz | 0.047a / 0.302a* | 0.057a / 0.037a | 0.080a / 0.103a | no |
| GABBR2 | scz | 0.052a / 0.023a | — / 0.065a | — / 0.048a | no |
| GGNBP2 | als | — / 0.947s | — / 0.376a | — / 0.829a | no |
| GLYCTK | scz | 0.068a / 0.423a | 0.133a / 0.035a | 0.189a / 0.647a | no |
| IRF3 | scz | 0.114a / 0.055a* | 0.085a / 0.576a* | 0.049a / 0.172a* | no |
| KLC1 | scz | 0.506a / 8.14e-08s | 0.029a / 0.003s | 0.721s / 0.001s | no |
| LPCAT4 | scz | 0.021a / 0.002s | 0.027a / 0.033a | 0.039a / 0.014a | no |
| MED19 | scz | — / 0.038a | — / 0.038a | 0.051a / 0.017a* | no |
| MRPS33 | scz | 0.046a / 0.058a | 0.054a / 0.088a | 0.047a / 0.028a | no |
| NCOR1 | pd | 0.294a / 0.891s | 0.284a / 0.135a | 0.063a / 0.499a | no |
| NDUFAF7 | scz | 0.054a / 0.948s | 0.066a / 0.340a | 0.048a / 0.188a | no |
| NDUFS3 | ad | 0.036a / 7.03e-05s | 0.130a / 0.027a | 0.047a / 0.021a | no |
| NEK4 | scz | 0.069a / 0.003a | 0.119a / 0.368a* | 0.112a / 0.076a | no |
| NUP50 | scz | — / 0.687a* | — / 0.768a* | 0.035a / 0.372a* | no |
| POLG | scz | 0.057a / 0.146a | 0.043a / 0.081a | 0.057a / 0.038a | no |
| PPIL2 | scz | 0.858a / 0.664a* | 0.048a / 0.042a* | 0.027a / 0.095a | **yes** |
| PRDM2 | als | 0.031a / 0.027a | 0.032a / 0.003a | 0.050a / 0.058a* | no |
| PSMD6 | scz | 0.055a / 0.611s | 0.041a / 0.486a | 0.039a / 0.027a* | no |
| PTPRN | als | 0.044a / 0.834s | 0.038a / 0.029a* | 0.075a / 0.388a | no |
| RAI1 | scz | — / 0.012a | — / 0.024a | — / 0.019a | no |
| RASA1 | scz | 0.031a / 0.283a | 0.029a / 0.112a | 0.031a / 0.041a | no |
| RERE | scz | 0.031a / 0.612a | 0.178a / 0.510a | 0.239a / 0.037a | no |
| SH3GL2 | pd | — / 0.001a* | — / 0.964a* | — / 0.024a | no |
| SIRPA | ad | 0.026a / 0.960a* | — / 0.961a* | — / 0.946a* | no |
| SNAP91 | scz | 0.035a / 0.963s* | 0.034a / 0.966s | 0.032a / 0.959a | no |
| SYT5 | scz | 0.027a / 0.018a | 0.038a / 0.004a* | 0.020a / 0.230a | no |
| TAOK2 | scz | 0.034a / 0.933a | 0.005s / 0.024a | 0.094a / 0.064a | no |
| TPCN1 | ad | 0.070a / 0.215a | 0.057a / 0.104a | 0.066a / 0.067a | no |
| TTC19 | pd | 0.020a / 0.737s* | 0.047a / 0.575a | 0.039a / 0.091a | no |
| TXNDC15 | als | — / 0.028a* | — / 0.288a | 0.123a / 0.317a | no |
| UNC13A | als | — / 0.009a | — / 0.027a | — / 0.066a | no |
| WIPI2 | als | 0.088a / 0.382a | 0.058a / 0.048a | 0.036a / 0.025a | no |
| YPEL1 | scz | — / 0.095a* | — / 0.601a* | — / 0.619a | no |
| YWHAB | scz | 0.837a / 0.845a* | — / 0.856a* | — / 0.660a* | **yes** |
| ZNF232 | ad | 0.041a / 0.046a | 0.052a / 0.167a | 0.069a / 1.71e-08a | no |

## BrainSEQ's own switch-axis nominations

7 (analysis, locus, gene) cells reach PP4_S_g >= 0.8 in at least one region; 4 are also GTEx signal-level nominations. The count is a set of locus nominations selected on S_g, not a switch-specificity estimate; the paired test above is the unbiased comparison.

| gene | trait | PP4 S_g | PP4 A_g | regions coloc | best region | GTEx nomination |
|---|---|---|---|---|---|---|
| STX4 | pd | 0.942 | 0.332 | 1/3 | caudate | no |
| ZNF592 | scz | 0.888 | 0.125 | 1/3 | hippocampus | no |
| ACTR1B | scz | 0.858 | 0.569 | 1/1 | hippocampus | yes |
| ACTR1B | scz | 0.858 | 0.569 | 1/1 | hippocampus | yes |
| PPIL2 | scz | 0.858 | 0.664 | 1/3 | caudate | yes |
| YWHAB | scz | 0.837 | 0.845 | 1/1 | caudate | yes |
| ARL14EP | scz | 0.834 | 0.794 | 1/3 | caudate | no |

## Reading rules

- An S_g colocalization shares a variant with a gene's switch coordinate (PC1 of its within-gene composition). It names no intron or event; the GTEx all-introns arm and the event audit do that.
- The PP4 per gene is a maximum over three regions; quote it with the region count.
- PP4_S_g >= 0.8 with low PP4_A_g is switch-preferential colocalization under prespecified thresholds, not a demonstrated switch-mediated mechanism.
- Region donors are 169-229, so a non-colocalizing A_g or S_g can be a power result.

