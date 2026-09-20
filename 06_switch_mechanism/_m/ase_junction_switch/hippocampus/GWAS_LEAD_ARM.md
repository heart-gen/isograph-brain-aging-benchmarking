# Allelic switch test anchored on the GWAS lead — hippocampus

The primary arm fits at the switch-QTL lead and transfers the sign to the disease allele through signed LD, which fails whenever the two variants are not in strong LD — most rows. This arm fits the same within-donor contrast **at the GWAS lead itself**, so the risk haplotype is read from the donor's phased genotype and the orientation is exact. It buys validity with power: only donors heterozygous at that variant are informative.

## Phase-frame check (run before anything is fitted)

The fragment counts are phased by phASER against the Quest VCF; the anchor genotypes come from the BrainSEQ TOPMed panel. Both must be the same frame or every re-anchored sign is arbitrary, so the fitted leads were extracted from the panel and compared donor by donor.

- variants compared: **22**, donor x variant pairs: **4727**
- genotype agreement (unphased): **1.0000**
- heterozygous pairs carrying frame information: **1477**
- **phase-frame agreement: 1.0000**
- gate: the arm runs only at >= 0.95

## Result

- pair x trait rows: **209** over 30 genes
- fitted: **145**; median donors informative on both haplotypes: **29.0**
- rows whose GWAS lead IS the fitted lead (nothing to re-anchor): **0**
- significant at q < 0.05: **5** (5 with an interpretable effect size, 0 at the estimator bound)
- risk allele raises T1 in **75** of the fitted rows (median |beta| among unbounded fits 0.439)

A fit **at the bound** is quasi-separation: one haplotype carries one isoform only. Its sign is informative and its magnitude is not, so those rows are counted apart rather than read as effects.

| status | rows |
| --- | ---: |
| fitted | 145 |
| anchor_not_in_panel | 27 |
| too_few_paired_het_donors | 23 |
| one_isoform_in_het_donors | 13 |
| no_ea_donors | 1 |

## Fitted rows

| gene | trait | anchor | risk | beta (risk hap) | p | q | donors paired |
| --- | --- | --- | :---: | ---: | ---: | ---: | ---: |
| KLC1 | scz | rs10873538 | G | 0.868 | 6.55e-06 | 0.000784 | 62 |
| KLC1 | scz | rs10873538 | G | 0.961 | 1.27e-05 | 0.000784 | 61 |
| KLC1 | scz | rs10873538 | G | -1.016 | 1.62e-05 | 0.000784 | 61 |
| KLC1 | scz | rs10873538 | G | 0.808 | 0.000243 | 0.00881 | 62 |
| PICALM | ad | rs10792832 | G | 0.416 | 0.00053 | 0.0154 | 46 |
| ZFYVE21 | scz | rs10873538 | G | -0.422 | 0.00738 | 0.178 | 72 |
| EFHB | scz | rs12489270 | C | at bound | 0.0106 | 0.197 | 24 |
| EFHB | scz | rs12489270 | C | at bound | 0.0119 | 0.197 | 23 |
| ZFYVE21 | scz | rs10873538 | G | 0.397 | 0.0122 | 0.197 | 71 |
| PICALM | ad | rs10792832 | G | 0.185 | 0.0148 | 0.215 | 46 |
| PICALM | ad | rs10792832 | G | 0.181 | 0.0169 | 0.221 | 46 |
| PICALM | ad | rs10792832 | G | 0.177 | 0.0183 | 0.221 | 50 |
| SNCA | lbd | rs7680557 | A | 2.205 | 0.0258 | 0.26 | 6 |
| SNCA | lbd | rs7680557 | A | 2.205 | 0.0258 | 0.26 | 6 |
| MYO19 | als | rs11650008 | C | -2.481 | 0.0295 | 0.26 | 23 |
| RPS6KL1 | als | rs3742783 | C | 0.657 | 0.0329 | 0.26 | 54 |
| EFHB | scz | rs12489270 | C | -9.604 | 0.0332 | 0.26 | 12 |
| EFHB | scz | rs12489270 | C | -9.604 | 0.0332 | 0.26 | 12 |
| RPS6KL1 | als | rs3742783 | C | 0.651 | 0.035 | 0.26 | 52 |
| SYT5 | scz | rs2365727 | G | 1.163 | 0.0359 | 0.26 | 9 |
| RPS6KL1 | als | rs3742783 | C | 0.247 | 0.0478 | 0.33 | 65 |
| MYO19 | als | rs11650008 | C | -1.661 | 0.0635 | 0.418 | 26 |
| MYO19 | als | rs11650008 | C | -1.735 | 0.0663 | 0.418 | 32 |
| RPS6KL1 | als | rs3742783 | C | 0.668 | 0.08 | 0.483 | 53 |
| MYO19 | als | rs11650008 | C | -1.110 | 0.0955 | 0.536 | 37 |
| TTC19 | pd | rs4566208 | A | at bound | 0.104 | 0.536 | 15 |
| SNAP91 | scz | rs217336 | C | -0.595 | 0.109 | 0.536 | 29 |
| SNAP91 | scz | rs217336 | C | 0.594 | 0.109 | 0.536 | 29 |
| SNAP91 | scz | rs217336 | C | -0.592 | 0.12 | 0.536 | 28 |
| SNAP91 | scz | rs217336 | C | 0.585 | 0.12 | 0.536 | 28 |
| SNAP91 | scz | rs217336 | C | 0.585 | 0.12 | 0.536 | 28 |
| EFHB | scz | rs12489270 | C | -1.585 | 0.126 | 0.536 | 41 |
| SNAP91 | scz | rs217336 | C | 0.572 | 0.128 | 0.536 | 29 |
| EDEM3 | scz | rs78444298 | G | 7.767 | 0.134 | 0.536 | 6 |
| EDEM3 | scz | rs78444298 | G | -7.530 | 0.134 | 0.536 | 6 |
| PICALM | ad | rs10792832 | G | -0.204 | 0.136 | 0.536 | 46 |
| EFHB | scz | rs12489270 | C | -1.504 | 0.139 | 0.536 | 50 |
| EFHB | scz | rs12489270 | C | -1.499 | 0.14 | 0.536 | 50 |
| PICALM | ad | rs10792832 | G | -0.194 | 0.157 | 0.555 | 46 |
| EFHB | scz | rs12489270 | C | -1.426 | 0.16 | 0.555 | 45 |
| EFHB | scz | rs12489270 | C | -1.426 | 0.16 | 0.555 | 45 |
| RERE | scz | rs3795310 | C | 1.308 | 0.161 | 0.555 | 31 |
| RCBTB1 | scz | rs7322886 | T | -1.566 | 0.165 | 0.555 | 53 |
| RERE | scz | rs3795310 | C | 0.925 | 0.176 | 0.567 | 31 |
| ATE1 | scz | rs7902292 | T | 5.461 | 0.176 | 0.567 | 5 |
| SYT5 | scz | rs2365727 | G | 0.640 | 0.208 | 0.653 | 9 |
| EDEM3 | scz | rs78444298 | G | -0.789 | 0.241 | 0.653 | 6 |
| EDEM3 | scz | rs78444298 | G | 0.797 | 0.241 | 0.653 | 6 |
| MAP7D1 | scz | rs11210892 | G | at bound | 0.249 | 0.653 | 45 |
| MAP7D1 | scz | rs11210892 | G | at bound | 0.249 | 0.653 | 45 |
| MAP7D1 | scz | rs11210892 | G | at bound | 0.249 | 0.653 | 45 |
| MAP7D1 | scz | rs11210892 | G | at bound | 0.249 | 0.653 | 45 |
| ATE1 | scz | rs7902292 | T | at bound | 0.251 | 0.653 | 6 |
| ATE1 | scz | rs7902292 | T | at bound | 0.253 | 0.653 | 6 |
| ARL14EP | scz | rs605765 | T | -3.325 | 0.254 | 0.653 | 28 |
| ATE1 | scz | rs7902292 | T | at bound | 0.258 | 0.653 | 7 |
| RPS6KL1 | als | rs3742783 | C | 0.738 | 0.261 | 0.653 | 64 |
| RCBTB1 | scz | rs7322886 | T | 0.729 | 0.267 | 0.653 | 74 |
| SCAMP4 | ad | rs12151021 | A | -6.180 | 0.269 | 0.653 | 9 |
| SCAMP4 | ad | rs12151021 | A | -1.365 | 0.271 | 0.653 | 10 |
| SCAMP4 | ad | rs12151021 | A | -3.498 | 0.284 | 0.653 | 13 |
| RCBTB1 | scz | rs7322886 | T | 0.780 | 0.286 | 0.653 | 74 |
| SCAMP4 | ad | rs12151021 | A | -3.572 | 0.288 | 0.653 | 13 |
| RBM6 | scz | rs2710323 | T | -5.946 | 0.288 | 0.653 | 61 |
| KLC1 | scz | rs10873538 | G | -1.275 | 0.296 | 0.661 | 7 |
| SYT5 | scz | rs2365727 | G | 0.398 | 0.306 | 0.672 | 10 |
| PICALM | ad | rs10792832 | G | -0.122 | 0.327 | 0.697 | 50 |
| RPS6KL1 | als | rs3742783 | C | 0.699 | 0.331 | 0.697 | 55 |
| SYT5 | scz | rs2365727 | G | 0.391 | 0.332 | 0.697 | 10 |
| EDEM3 | scz | rs78444298 | G | 0.440 | 0.37 | 0.717 | 7 |
| EDEM3 | scz | rs78444298 | G | -0.439 | 0.371 | 0.717 | 7 |
| EDEM3 | scz | rs78444298 | G | -0.439 | 0.371 | 0.717 | 7 |
| EDEM3 | scz | rs78444298 | G | -0.439 | 0.371 | 0.717 | 7 |
| RBM6 | scz | rs2710323 | T | -5.947 | 0.379 | 0.717 | 61 |
| RBM6 | scz | rs2710323 | T | -5.947 | 0.379 | 0.717 | 61 |
| RBM6 | scz | rs2710323 | T | -5.941 | 0.384 | 0.717 | 61 |
| RBM6 | scz | rs2710323 | T | -5.865 | 0.388 | 0.717 | 64 |
| RBM6 | scz | rs2710323 | T | -5.864 | 0.388 | 0.717 | 64 |
| ATE1 | scz | rs7902292 | T | 0.651 | 0.395 | 0.717 | 7 |
| RERE | scz | rs3795310 | C | 1.080 | 0.408 | 0.717 | 31 |
| TTC19 | pd | rs4566208 | A | 0.509 | 0.413 | 0.717 | 23 |
| SNAP91 | scz | rs217336 | C | 0.186 | 0.413 | 0.717 | 55 |
| SNAP91 | scz | rs217336 | C | 0.186 | 0.414 | 0.717 | 56 |
| MRPS33 | scz | rs6662 | T | -0.446 | 0.415 | 0.717 | 31 |
| PICALM | ad | rs10792832 | G | -0.063 | 0.439 | 0.738 | 61 |
| SNAP91 | scz | rs217336 | C | 0.265 | 0.445 | 0.738 | 55 |
| SNAP91 | scz | rs217336 | C | 0.265 | 0.445 | 0.738 | 55 |
| ATE1 | scz | rs7902292 | T | 0.520 | 0.452 | 0.738 | 6 |
| SNAP91 | scz | rs217336 | C | 0.285 | 0.453 | 0.738 | 56 |
| PICALM | ad | rs10792832 | G | -0.059 | 0.463 | 0.746 | 61 |
| PICALM | ad | rs10792832 | G | -0.058 | 0.473 | 0.747 | 61 |
| MRPS33 | scz | rs6662 | T | -0.393 | 0.474 | 0.747 | 29 |
| RPS6KL1 | als | rs3742783 | C | 0.838 | 0.484 | 0.747 | 56 |
| PICALM | ad | rs10792832 | G | -0.056 | 0.485 | 0.747 | 61 |
| EDEM3 | scz | rs78444298 | G | 0.271 | 0.5 | 0.763 | 7 |
| SCAMP4 | ad | rs12151021 | A | 0.657 | 0.51 | 0.771 | 5 |
| RPS6KL1 | als | rs3742783 | C | 0.735 | 0.545 | 0.803 | 64 |
| PRMT7 | scz | rs6499508 | A | -0.809 | 0.548 | 0.803 | 7 |
| MRPS33 | scz | rs6662 | T | 0.323 | 0.549 | 0.803 | 30 |
| PRMT7 | scz | rs6499508 | A | -0.363 | 0.582 | 0.842 | 10 |
| PRMT7 | scz | rs6499508 | A | -0.319 | 0.587 | 0.842 | 11 |
| SCAMP4 | ad | rs12151021 | A | 0.535 | 0.602 | 0.856 | 5 |
| MAP7D1 | scz | rs11210892 | G | -0.759 | 0.609 | 0.857 | 40 |
| MRPS33 | scz | rs6662 | T | -0.361 | 0.619 | 0.863 | 29 |
| MRPS33 | scz | rs6662 | T | -0.346 | 0.64 | 0.884 | 27 |
| KLC1 | scz | rs10873538 | G | 0.125 | 0.678 | 0.892 | 53 |
| DDRGK1 | pd | rs2295547 | A | 0.090 | 0.684 | 0.892 | 47 |
| RPS6KL1 | als | rs3742783 | C | 0.295 | 0.687 | 0.892 | 50 |
| ATE1 | scz | rs7902292 | T | 0.143 | 0.69 | 0.892 | 8 |
| MRPS33 | scz | rs6662 | T | -0.278 | 0.69 | 0.892 | 30 |
| SYT5 | scz | rs2365727 | G | -0.366 | 0.694 | 0.892 | 10 |
| SCAMP4 | ad | rs12151021 | A | 0.390 | 0.695 | 0.892 | 13 |
| MYO19 | als | rs11650008 | C | -0.271 | 0.701 | 0.892 | 18 |
| SYT5 | scz | rs2365727 | G | -0.086 | 0.702 | 0.892 | 10 |
| MRPS33 | scz | rs6662 | T | 0.253 | 0.719 | 0.896 | 27 |
| RBM6 | scz | rs2710323 | T | 0.244 | 0.722 | 0.896 | 64 |
| DDRGK1 | pd | rs2295547 | A | -0.081 | 0.727 | 0.896 | 36 |
| TTC19 | pd | rs4566208 | A | -0.127 | 0.736 | 0.896 | 32 |
| TTC19 | pd | rs4566208 | A | -0.115 | 0.757 | 0.896 | 33 |
| ATE1 | scz | rs7902292 | T | 0.254 | 0.758 | 0.896 | 7 |
| PICALM | ad | rs10792832 | G | -0.089 | 0.758 | 0.896 | 57 |
| EDEM3 | scz | rs78444298 | G | 0.142 | 0.76 | 0.896 | 6 |
| EDEM3 | scz | rs78444298 | G | 0.142 | 0.76 | 0.896 | 6 |
| PRMT7 | scz | rs6499508 | A | -0.327 | 0.799 | 0.924 | 7 |
| TTC19 | pd | rs4566208 | A | 0.201 | 0.811 | 0.924 | 26 |
| RBM6 | scz | rs2710323 | T | 0.132 | 0.812 | 0.924 | 64 |
| MYO19 | als | rs11650008 | C | 0.166 | 0.818 | 0.924 | 17 |
| MYO19 | als | rs11650008 | C | -0.166 | 0.818 | 0.924 | 17 |
| KLC1 | scz | rs10873538 | G | 0.180 | 0.826 | 0.924 | 51 |
| PRMT7 | scz | rs6499508 | A | 0.124 | 0.828 | 0.924 | 11 |
| ATE1 | scz | rs7902292 | T | 0.069 | 0.842 | 0.925 | 8 |
| PRMT7 | scz | rs6499508 | A | 0.135 | 0.842 | 0.925 | 10 |
| KLC1 | scz | rs10873538 | G | 0.060 | 0.855 | 0.926 | 53 |
| ATE1 | scz | rs7902292 | T | 0.059 | 0.863 | 0.926 | 8 |
| KLC1 | scz | rs10873538 | G | -0.132 | 0.867 | 0.926 | 51 |
| TTC19 | pd | rs4566208 | A | 0.731 | 0.878 | 0.926 | 16 |
| EDEM3 | scz | rs78444298 | G | 0.060 | 0.881 | 0.926 | 8 |
| EDEM3 | scz | rs78444298 | G | 0.060 | 0.881 | 0.926 | 8 |
| SCAMP4 | ad | rs12151021 | A | -0.100 | 0.922 | 0.954 | 13 |
| PRMT7 | scz | rs6499508 | A | -0.062 | 0.923 | 0.954 | 10 |
| KLC1 | scz | rs10873538 | G | -0.036 | 0.928 | 0.954 | 6 |
| ATE1 | scz | rs7902292 | T | -0.007 | 0.973 | 0.99 | 7 |
| SYT5 | scz | rs2365727 | G | 0.010 | 0.976 | 0.99 | 10 |
| SYT5 | scz | rs2365727 | G | 0.004 | 0.989 | 0.991 | 10 |
| PRMT7 | scz | rs6499508 | A | 0.008 | 0.991 | 0.991 | 10 |

## How to read a null here

A null in the primary arm was ambiguous: the sign could not be transferred for most rows, so 'no disease linkage' and 'no answer' looked alike. A null **here** is not ambiguous in that way — the orientation is exact by construction — but it is a power statement. Read `n_paired_het` before reading the p-value: the GWAS lead is chosen for disease signal, not for heterozygosity in 200 EA donors, and a row with a handful of informative donors has not tested anything.
