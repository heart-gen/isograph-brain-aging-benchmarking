# Allelic switch test anchored on the GWAS lead — dlpfc

The primary arm fits at the switch-QTL lead and transfers the sign to the disease allele through signed LD, which fails whenever the two variants are not in strong LD — most rows. This arm fits the same within-donor contrast **at the GWAS lead itself**, so the risk haplotype is read from the donor's phased genotype and the orientation is exact. It buys validity with power: only donors heterozygous at that variant are informative.

## Phase-frame check (run before anything is fitted)

The fragment counts are phased by phASER against the Quest VCF; the anchor genotypes come from the BrainSEQ TOPMed panel. Both must be the same frame or every re-anchored sign is arbitrary, so the fitted leads were extracted from the panel and compared donor by donor.

- variants compared: **10**, donor x variant pairs: **2068**
- genotype agreement (unphased): **1.0000**
- heterozygous pairs carrying frame information: **550**
- **phase-frame agreement: 0.9982**
- gate: the arm runs only at >= 0.95

## Result

- pair x trait rows: **97** over 13 genes
- fitted: **67**; median donors informative on both haplotypes: **25.0**
- rows whose GWAS lead IS the fitted lead (nothing to re-anchor): **0**
- significant at q < 0.05: **6** (4 with an interpretable effect size, 2 at the estimator bound)
- risk allele raises T1 in **39** of the fitted rows (median |beta| among unbounded fits 0.321)

A fit **at the bound** is quasi-separation: one haplotype carries one isoform only. Its sign is informative and its magnitude is not, so those rows are counted apart rather than read as effects.

| status | rows |
| --- | ---: |
| fitted | 67 |
| too_few_paired_het_donors | 16 |
| anchor_not_in_panel | 14 |

## Fitted rows

| gene | trait | anchor | risk | beta (risk hap) | p | q | donors paired |
| --- | --- | --- | :---: | ---: | ---: | ---: | ---: |
| NT5C2 | scz | rs11191580 | T | 1.308 | 0.000384 | 0.0257 | 29 |
| NT5C2 | scz | rs11191580 | T | at bound | 0.00221 | 0.0354 | 5 |
| NT5C2 | scz | rs11191580 | T | at bound | 0.00279 | 0.0354 | 5 |
| NT5C2 | scz | rs11191580 | T | 1.745 | 0.00281 | 0.0354 | 14 |
| NEK4 | scz | rs2710323 | T | 1.683 | 0.00317 | 0.0354 | 92 |
| NEK4 | scz | rs2710323 | T | -1.691 | 0.00317 | 0.0354 | 92 |
| SCFD1 | als | rs229243 | A | -0.208 | 0.00759 | 0.0727 | 86 |
| NEK4 | scz | rs2710323 | T | -1.707 | 0.0111 | 0.0826 | 91 |
| NEK4 | scz | rs2710323 | T | -1.707 | 0.0111 | 0.0826 | 91 |
| POLG | scz | rs4702 | G | 9.578 | 0.0326 | 0.199 | 18 |
| POLG | scz | rs4702 | G | 9.578 | 0.0326 | 0.199 | 18 |
| SCFD1 | als | rs229243 | A | -0.157 | 0.0429 | 0.239 | 86 |
| SCFD1 | als | rs229243 | A | 0.161 | 0.0477 | 0.246 | 86 |
| SCFD1 | als | rs229243 | A | -0.153 | 0.056 | 0.259 | 86 |
| NT5C2 | scz | rs11191580 | T | -1.030 | 0.061 | 0.259 | 14 |
| SCFD1 | als | rs229243 | A | 0.159 | 0.0622 | 0.259 | 86 |
| POLG | scz | rs4702 | G | -1.754 | 0.0656 | 0.259 | 19 |
| SCFD1 | als | rs229243 | A | 0.146 | 0.0853 | 0.294 | 86 |
| PRDM2 | als | rs2744680 | G | 0.392 | 0.0917 | 0.294 | 48 |
| NEK4 | scz | rs2710323 | T | 0.643 | 0.0919 | 0.294 | 90 |
| PRDM2 | als | rs2744680 | G | 0.391 | 0.092 | 0.294 | 48 |
| NT5C2 | scz | rs11191580 | T | 1.069 | 0.101 | 0.296 | 14 |
| SCFD1 | als | rs229243 | A | 0.142 | 0.103 | 0.296 | 86 |
| SCFD1 | als | rs229243 | A | 0.141 | 0.106 | 0.296 | 86 |
| NT5C2 | scz | rs11191580 | T | 0.615 | 0.13 | 0.349 | 27 |
| EDEM3 | scz | rs78444298 | G | -0.671 | 0.136 | 0.351 | 8 |
| NT5C2 | scz | rs11191580 | T | -0.794 | 0.157 | 0.364 | 6 |
| EDEM3 | scz | rs78444298 | G | -0.729 | 0.157 | 0.364 | 9 |
| EDEM3 | scz | rs78444298 | G | -0.729 | 0.157 | 0.364 | 9 |
| PRDM2 | als | rs2744680 | G | 0.362 | 0.193 | 0.431 | 10 |
| TMED4 | scz | rs6656 | T | -5.073 | 0.237 | 0.512 | 82 |
| TMED4 | scz | rs6656 | T | 1.211 | 0.258 | 0.524 | 85 |
| TMED4 | scz | rs6656 | T | 1.211 | 0.258 | 0.524 | 85 |
| NT5C2 | scz | rs11191580 | T | -0.587 | 0.277 | 0.547 | 6 |
| EDEM3 | scz | rs78444298 | G | -0.503 | 0.318 | 0.593 | 7 |
| EDEM3 | scz | rs78444298 | G | -0.503 | 0.318 | 0.593 | 7 |
| NT5C2 | scz | rs11191580 | T | 0.509 | 0.348 | 0.607 | 6 |
| PRDM2 | als | rs2744680 | G | -0.250 | 0.353 | 0.607 | 48 |
| PRDM2 | als | rs2744680 | G | 0.250 | 0.354 | 0.607 | 48 |
| EDEM3 | scz | rs78444298 | G | -0.287 | 0.389 | 0.607 | 9 |
| EDEM3 | scz | rs78444298 | G | 0.286 | 0.389 | 0.607 | 9 |
| EDEM3 | scz | rs78444298 | G | 0.286 | 0.389 | 0.607 | 9 |
| EDEM3 | scz | rs78444298 | G | 0.286 | 0.389 | 0.607 | 9 |
| PRDM2 | als | rs2744680 | G | -0.258 | 0.556 | 0.806 | 48 |
| PRDM2 | als | rs2744680 | G | 0.258 | 0.558 | 0.806 | 48 |
| POLG | scz | rs4702 | G | 0.692 | 0.572 | 0.806 | 19 |
| NT5C2 | scz | rs11191580 | T | -0.332 | 0.574 | 0.806 | 25 |
| PRDM2 | als | rs2744680 | G | -0.125 | 0.578 | 0.806 | 47 |
| EDEM3 | scz | rs78444298 | G | 0.244 | 0.607 | 0.814 | 8 |
| EDEM3 | scz | rs78444298 | G | 0.244 | 0.607 | 0.814 | 8 |
| PRDM2 | als | rs2744680 | G | 0.123 | 0.651 | 0.851 | 46 |
| POLG | scz | rs4702 | G | 0.329 | 0.673 | 0.851 | 19 |
| POLG | scz | rs4702 | G | 0.329 | 0.673 | 0.851 | 19 |
| TMED4 | scz | rs6656 | T | -0.321 | 0.691 | 0.857 | 58 |
| DGKZ | scz | rs12285419 | A | 0.416 | 0.773 | 0.939 | 14 |
| NEK4 | scz | rs2710323 | T | 0.302 | 0.823 | 0.939 | 90 |
| EDEM3 | scz | rs78444298 | G | -0.071 | 0.851 | 0.939 | 9 |
| EDEM3 | scz | rs78444298 | G | 0.071 | 0.851 | 0.939 | 9 |
| NEK4 | scz | rs2710323 | T | 0.220 | 0.87 | 0.939 | 94 |
| NEK4 | scz | rs2710323 | T | 0.231 | 0.873 | 0.939 | 94 |
| PRDM2 | als | rs2744680 | G | -0.041 | 0.879 | 0.939 | 46 |
| PRDM2 | als | rs2744680 | G | -0.041 | 0.879 | 0.939 | 46 |
| NT5C2 | scz | rs11191580 | T | -0.081 | 0.884 | 0.939 | 25 |
| PRDM2 | als | rs2744680 | G | 0.021 | 0.897 | 0.939 | 47 |
| NT5C2 | scz | rs11191580 | T | 0.107 | 0.912 | 0.94 | 13 |
| NT5C2 | scz | rs11191580 | T | -0.038 | 0.944 | 0.958 | 25 |
| POLG | scz | rs4702 | G | 0.005 | 0.99 | 0.99 | 19 |

## How to read a null here

A null in the primary arm was ambiguous: the sign could not be transferred for most rows, so 'no disease linkage' and 'no answer' looked alike. A null **here** is not ambiguous in that way — the orientation is exact by construction — but it is a power statement. Read `n_paired_het` before reading the p-value: the GWAS lead is chosen for disease signal, not for heterozygosity in 200 EA donors, and a row with a handful of informative donors has not tested anything.
