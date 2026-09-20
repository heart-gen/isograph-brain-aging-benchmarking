# Allelic switch test anchored on the GWAS lead — caudate

The primary arm fits at the switch-QTL lead and transfers the sign to the disease allele through signed LD, which fails whenever the two variants are not in strong LD — most rows. This arm fits the same within-donor contrast **at the GWAS lead itself**, so the risk haplotype is read from the donor's phased genotype and the orientation is exact. It buys validity with power: only donors heterozygous at that variant are informative.

## Phase-frame check (run before anything is fitted)

The fragment counts are phased by phASER against the Quest VCF; the anchor genotypes come from the BrainSEQ TOPMed panel. Both must be the same frame or every re-anchored sign is arbitrary, so the fitted leads were extracted from the panel and compared donor by donor.

- variants compared: **27**, donor x variant pairs: **6416**
- genotype agreement (unphased): **1.0000**
- heterozygous pairs carrying frame information: **2694**
- **phase-frame agreement: 0.9996**
- gate: the arm runs only at >= 0.95

## Result

- pair x trait rows: **250** over 35 genes
- fitted: **159**; median donors informative on both haplotypes: **39.0**
- rows whose GWAS lead IS the fitted lead (nothing to re-anchor): **0**
- significant at q < 0.05: **1** (0 with an interpretable effect size, 1 at the estimator bound)
- risk allele raises T1 in **93** of the fitted rows (median |beta| among unbounded fits 0.270)

A fit **at the bound** is quasi-separation: one haplotype carries one isoform only. Its sign is informative and its magnitude is not, so those rows are counted apart rather than read as effects.

| status | rows |
| --- | ---: |
| fitted | 159 |
| anchor_not_in_panel | 69 |
| one_isoform_in_het_donors | 19 |
| too_few_paired_het_donors | 3 |

## Fitted rows

| gene | trait | anchor | risk | beta (risk hap) | p | q | donors paired |
| --- | --- | --- | :---: | ---: | ---: | ---: | ---: |
| RPS6KL1 | als | rs3742783 | C | at bound | 2.16e-05 | 0.00344 | 74 |
| MRPS33 | scz | rs6662 | T | 1.166 | 0.0119 | 0.683 | 22 |
| TTC19 | pd | rs4566208 | A | -1.060 | 0.0183 | 0.683 | 50 |
| TTC19 | pd | rs4566208 | A | -1.111 | 0.0186 | 0.683 | 53 |
| SCAMP4 | ad | rs12151021 | A | 4.142 | 0.027 | 0.683 | 17 |
| FNBP1 | als | rs7870884 | T | 1.376 | 0.0273 | 0.683 | 14 |
| SCFD1 | als | rs229243 | A | -0.185 | 0.0301 | 0.683 | 90 |
| TMED4 | scz | rs6656 | T | 1.126 | 0.0389 | 0.683 | 105 |
| TMED4 | scz | rs6656 | T | 1.126 | 0.0389 | 0.683 | 105 |
| SCAMP4 | ad | rs12151021 | A | 5.041 | 0.0472 | 0.683 | 14 |
| SCFD1 | als | rs229243 | A | -0.153 | 0.0618 | 0.683 | 90 |
| NME4 | als | rs11645554 | T | -1.454 | 0.0642 | 0.683 | 13 |
| MRPS33 | scz | rs6662 | T | 0.734 | 0.0676 | 0.683 | 28 |
| MRPS33 | scz | rs6662 | T | -0.618 | 0.0782 | 0.683 | 28 |
| PICALM | ad | rs10792832 | G | 0.131 | 0.0845 | 0.683 | 58 |
| ZFYVE21 | scz | rs10873538 | G | -0.126 | 0.0952 | 0.683 | 105 |
| COPA | scz | rs1289003 | T | 0.526 | 0.0995 | 0.683 | 65 |
| RPS6KL1 | als | rs3742783 | C | at bound | 0.11 | 0.683 | 79 |
| ZNF592 | scz | rs2456020 | C | -7.952 | 0.118 | 0.683 | 18 |
| PICALM | ad | rs10792832 | G | 0.122 | 0.12 | 0.683 | 73 |
| COPA | scz | rs1289003 | T | 0.956 | 0.12 | 0.683 | 65 |
| SCFD1 | als | rs229243 | A | -0.134 | 0.128 | 0.683 | 90 |
| PAK6 | scz | rs56205728 | A | 0.963 | 0.131 | 0.683 | 56 |
| EFHB | scz | rs12489270 | C | 1.093 | 0.139 | 0.683 | 60 |
| ZFYVE21 | scz | rs10873538 | G | 0.112 | 0.14 | 0.683 | 105 |
| EFHB | scz | rs12489270 | C | 1.086 | 0.14 | 0.683 | 59 |
| PAK6 | scz | rs56205728 | A | 1.008 | 0.151 | 0.683 | 57 |
| FNBP1 | als | rs7870884 | T | -1.304 | 0.156 | 0.683 | 11 |
| EFHB | scz | rs12489270 | C | 1.139 | 0.161 | 0.683 | 48 |
| BCKDK | ad | rs1140239 | C | 0.741 | 0.165 | 0.683 | 67 |
| RPS6KL1 | als | rs3742783 | C | 0.384 | 0.166 | 0.683 | 79 |
| EFHB | scz | rs12489270 | C | 1.033 | 0.168 | 0.683 | 50 |
| EFHB | scz | rs12489270 | C | 1.033 | 0.168 | 0.683 | 50 |
| SNAP91 | scz | rs217336 | C | 0.397 | 0.17 | 0.683 | 39 |
| SNAP91 | scz | rs217336 | C | 0.397 | 0.17 | 0.683 | 39 |
| SNAP91 | scz | rs217336 | C | 0.397 | 0.17 | 0.683 | 39 |
| SNAP91 | scz | rs217336 | C | 0.396 | 0.171 | 0.683 | 39 |
| SNAP91 | scz | rs217336 | C | 0.389 | 0.175 | 0.683 | 39 |
| SNAP91 | scz | rs217336 | C | 0.389 | 0.175 | 0.683 | 39 |
| TAOK2 | scz | rs3814883 | C | 8.268 | 0.18 | 0.683 | 139 |
| SCAMP4 | ad | rs12151021 | A | 0.951 | 0.193 | 0.683 | 7 |
| TAOK2 | scz | rs3814883 | C | -8.608 | 0.194 | 0.683 | 139 |
| FNBP1 | als | rs7870884 | T | -6.643 | 0.212 | 0.683 | 12 |
| SCFD1 | als | rs229243 | A | 0.119 | 0.217 | 0.683 | 90 |
| POLG | scz | rs4702 | G | at bound | 0.219 | 0.683 | 22 |
| POLG | scz | rs4702 | G | at bound | 0.219 | 0.683 | 22 |
| POLG | scz | rs4702 | G | at bound | 0.219 | 0.683 | 22 |
| TAOK2 | scz | rs3814883 | C | 6.380 | 0.223 | 0.683 | 139 |
| BCKDK | ad | rs1140239 | C | 8.682 | 0.234 | 0.683 | 66 |
| FNBP1 | als | rs7870884 | T | -0.728 | 0.234 | 0.683 | 21 |
| CCS | scz | rs58950470 | T | 9.694 | 0.236 | 0.683 | 22 |
| POLG | scz | rs4702 | G | 2.767 | 0.241 | 0.683 | 22 |
| FNBP1 | als | rs7870884 | T | 0.218 | 0.242 | 0.683 | 33 |
| SETD6 | scz | rs149165 | T | 0.731 | 0.243 | 0.683 | 31 |
| PAK6 | scz | rs56205728 | A | -1.050 | 0.243 | 0.683 | 56 |
| FNBP1 | als | rs7870884 | T | 0.200 | 0.243 | 0.683 | 28 |
| PAK6 | scz | rs56205728 | A | 1.042 | 0.245 | 0.683 | 56 |
| CCS | scz | rs58950470 | T | 8.134 | 0.254 | 0.69 | 22 |
| ZNF592 | scz | rs2456020 | C | -9.989 | 0.258 | 0.69 | 18 |
| FNBP1 | als | rs7870884 | T | 0.189 | 0.265 | 0.69 | 34 |
| ARL14EP | scz | rs605765 | T | -1.189 | 0.267 | 0.69 | 36 |
| FNBP1 | als | rs7870884 | T | 0.208 | 0.28 | 0.69 | 29 |
| MRPS33 | scz | rs6662 | T | 0.800 | 0.283 | 0.69 | 10 |
| AZI2 | pd | rs6808178 | T | -9.737 | 0.283 | 0.69 | 25 |
| FNBP1 | als | rs7870884 | T | 0.586 | 0.283 | 0.69 | 19 |
| FNBP1 | als | rs7870884 | T | -0.732 | 0.287 | 0.69 | 21 |
| COPA | scz | rs1289003 | T | -0.396 | 0.291 | 0.69 | 7 |
| FNBP1 | als | rs7870884 | T | 0.183 | 0.298 | 0.696 | 30 |
| AZI2 | pd | rs6808178 | T | -4.321 | 0.302 | 0.696 | 24 |
| AZI2 | pd | rs6808178 | T | -3.805 | 0.308 | 0.698 | 24 |
| AZI2 | pd | rs6808178 | T | 2.658 | 0.313 | 0.698 | 25 |
| TTC19 | pd | rs4566208 | A | -8.745 | 0.32 | 0.698 | 41 |
| SCFD1 | als | rs229243 | A | -0.090 | 0.321 | 0.698 | 90 |
| PAK6 | scz | rs56205728 | A | 0.485 | 0.328 | 0.704 | 57 |
| SCFD1 | als | rs229243 | A | 0.089 | 0.365 | 0.774 | 90 |
| BCKDK | ad | rs1140239 | C | -0.239 | 0.371 | 0.775 | 66 |
| RPS6KL1 | als | rs3742783 | C | at bound | 0.403 | 0.832 | 67 |
| MRPS33 | scz | rs6662 | T | 0.223 | 0.417 | 0.834 | 68 |
| RCSD1 | als | rs1214598 | G | -7.730 | 0.42 | 0.834 | 17 |
| PAK6 | scz | rs56205728 | A | 0.244 | 0.42 | 0.834 | 53 |
| ATE1 | scz | rs7902292 | T | -0.152 | 0.437 | 0.834 | 15 |
| SCFD1 | als | rs229243 | A | 0.076 | 0.44 | 0.834 | 90 |
| EFHB | scz | rs12489270 | C | 0.928 | 0.449 | 0.834 | 19 |
| EFHB | scz | rs12489270 | C | -0.926 | 0.449 | 0.834 | 19 |
| SCFD1 | als | rs229243 | A | 0.075 | 0.45 | 0.834 | 90 |
| RCSD1 | als | rs1214598 | G | -7.159 | 0.471 | 0.834 | 28 |
| MRPS33 | scz | rs6662 | T | -0.195 | 0.472 | 0.834 | 67 |
| RCSD1 | als | rs1214598 | G | 7.233 | 0.473 | 0.834 | 27 |
| EFHB | scz | rs12489270 | C | 0.857 | 0.473 | 0.834 | 31 |
| MRPS33 | scz | rs6662 | T | 0.190 | 0.478 | 0.834 | 67 |
| SCAMP4 | ad | rs12151021 | A | -0.507 | 0.483 | 0.834 | 17 |
| MAP7D1 | scz | rs11210892 | G | -0.619 | 0.493 | 0.834 | 54 |
| EFHB | scz | rs12489270 | C | 0.811 | 0.496 | 0.834 | 31 |
| PAK6 | scz | rs56205728 | A | -0.216 | 0.498 | 0.834 | 53 |
| MRPS33 | scz | rs6662 | T | 0.184 | 0.503 | 0.834 | 66 |
| ATE1 | scz | rs7902292 | T | -0.159 | 0.514 | 0.834 | 29 |
| SNAP91 | scz | rs217336 | C | 0.115 | 0.524 | 0.834 | 77 |
| SNAP91 | scz | rs217336 | C | 0.115 | 0.524 | 0.834 | 77 |
| SNAP91 | scz | rs217336 | C | 0.115 | 0.524 | 0.834 | 77 |
| FNBP1 | als | rs7870884 | T | -0.101 | 0.534 | 0.834 | 33 |
| BCKDK | ad | rs1140239 | C | -0.127 | 0.541 | 0.834 | 66 |
| TTC19 | pd | rs4566208 | A | -0.162 | 0.542 | 0.834 | 65 |
| MRPS33 | scz | rs6662 | T | -0.204 | 0.545 | 0.834 | 32 |
| TTC19 | pd | rs4566208 | A | -0.156 | 0.55 | 0.834 | 65 |
| ATE1 | scz | rs7902292 | T | 0.127 | 0.551 | 0.834 | 30 |
| PICALM | ad | rs10792832 | G | 0.067 | 0.576 | 0.863 | 58 |
| FNBP1 | als | rs7870884 | T | 0.099 | 0.585 | 0.863 | 29 |
| PICALM | ad | rs10792832 | G | 0.062 | 0.59 | 0.863 | 74 |
| PAK6 | scz | rs56205728 | A | 0.437 | 0.597 | 0.863 | 54 |
| PAK6 | scz | rs56205728 | A | 0.437 | 0.597 | 0.863 | 54 |
| MRPS33 | scz | rs6662 | T | 0.290 | 0.604 | 0.865 | 23 |
| RPS6KL1 | als | rs3742783 | C | 0.270 | 0.635 | 0.9 | 81 |
| ATE1 | scz | rs7902292 | T | -0.146 | 0.642 | 0.9 | 28 |
| ZFYVE21 | scz | rs10873538 | G | -0.050 | 0.672 | 0.9 | 104 |
| SCAMP4 | ad | rs12151021 | A | 0.496 | 0.679 | 0.9 | 10 |
| SNAP91 | scz | rs217336 | C | -0.047 | 0.699 | 0.9 | 77 |
| NME4 | als | rs11645554 | T | -0.501 | 0.699 | 0.9 | 8 |
| BCKDK | ad | rs1140239 | C | 0.121 | 0.7 | 0.9 | 66 |
| SCAMP4 | ad | rs12151021 | A | 0.289 | 0.701 | 0.9 | 14 |
| SNAP91 | scz | rs217336 | C | 0.040 | 0.702 | 0.9 | 77 |
| MRPS33 | scz | rs6662 | T | -0.084 | 0.703 | 0.9 | 67 |
| CNOT7 | scz | rs876983 | T | -0.196 | 0.72 | 0.9 | 29 |
| NME4 | als | rs11645554 | T | -2.250 | 0.728 | 0.9 | 6 |
| ATE1 | scz | rs7902292 | T | -0.106 | 0.733 | 0.9 | 28 |
| ATE1 | scz | rs7902292 | T | -0.106 | 0.733 | 0.9 | 28 |
| CNOT7 | scz | rs876983 | T | -0.352 | 0.733 | 0.9 | 24 |
| POLG | scz | rs4702 | G | 0.160 | 0.737 | 0.9 | 23 |
| POLG | scz | rs4702 | G | 0.160 | 0.737 | 0.9 | 23 |
| POLG | scz | rs4702 | G | 0.160 | 0.737 | 0.9 | 23 |
| CNOT7 | scz | rs876983 | T | -0.178 | 0.745 | 0.9 | 29 |
| MRPS33 | scz | rs6662 | T | 0.070 | 0.752 | 0.9 | 68 |
| FNBP1 | als | rs7870884 | T | -0.056 | 0.752 | 0.9 | 33 |
| ATE1 | scz | rs7902292 | T | 0.105 | 0.769 | 0.9 | 29 |
| ATE1 | scz | rs7902292 | T | 0.105 | 0.769 | 0.9 | 29 |
| CNOT7 | scz | rs876983 | T | -0.300 | 0.77 | 0.9 | 24 |
| MAP7D1 | scz | rs11210892 | G | -0.044 | 0.774 | 0.9 | 56 |
| PAK6 | scz | rs56205728 | A | 0.168 | 0.788 | 0.9 | 54 |
| PAK6 | scz | rs56205728 | A | 0.168 | 0.788 | 0.9 | 54 |
| NME4 | als | rs11645554 | T | 1.668 | 0.792 | 0.9 | 6 |
| FNBP1 | als | rs7870884 | T | -0.229 | 0.794 | 0.9 | 6 |
| POLG | scz | rs4702 | G | 0.115 | 0.798 | 0.9 | 23 |
| TTC19 | pd | rs4566208 | A | -0.197 | 0.821 | 0.908 | 21 |
| PICALM | ad | rs10792832 | G | -0.059 | 0.821 | 0.908 | 56 |
| PICALM | ad | rs10792832 | G | 0.028 | 0.823 | 0.908 | 84 |
| BCKDK | ad | rs1140239 | C | 0.070 | 0.828 | 0.908 | 66 |
| ATE1 | scz | rs7902292 | T | -0.064 | 0.848 | 0.923 | 25 |
| COPA | scz | rs1289003 | T | 0.038 | 0.893 | 0.955 | 7 |
| ARL14EP | scz | rs605765 | T | -0.576 | 0.898 | 0.955 | 36 |
| MAP7D1 | scz | rs11210892 | G | -0.009 | 0.898 | 0.955 | 56 |
| ARL14EP | scz | rs605765 | T | -0.514 | 0.909 | 0.955 | 36 |
| COPA | scz | rs1289003 | T | 0.040 | 0.913 | 0.955 | 65 |
| SCAMP4 | ad | rs12151021 | A | 0.071 | 0.913 | 0.955 | 17 |
| MRPS33 | scz | rs6662 | T | -0.017 | 0.937 | 0.961 | 69 |
| MAP7D1 | scz | rs11210892 | G | -0.008 | 0.937 | 0.961 | 56 |
| MAP7D1 | scz | rs11210892 | G | -0.008 | 0.937 | 0.961 | 56 |
| PICALM | ad | rs10792832 | G | -0.010 | 0.943 | 0.961 | 58 |
| RPS6KL1 | als | rs3742783 | C | -0.013 | 0.976 | 0.983 | 79 |
| RPS6KL1 | als | rs3742783 | C | 0.013 | 0.976 | 0.983 | 79 |
| MAP7D1 | scz | rs11210892 | G | -0.135 | 0.994 | 0.994 | 12 |

## How to read a null here

A null in the primary arm was ambiguous: the sign could not be transferred for most rows, so 'no disease linkage' and 'no answer' looked alike. A null **here** is not ambiguous in that way — the orientation is exact by construction — but it is a power statement. Read `n_paired_het` before reading the p-value: the GWAS lead is chosen for disease signal, not for heterozygosity in 200 EA donors, and a row with a handful of informative donors has not tested anything.
