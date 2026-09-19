# Allele-specific switch test — caudate

Within a donor heterozygous at the switch-QTL lead, the two haplotypes share a nucleus, a cell-type mixture and an environment, so a cis effect of the lead on isoform choice shows up as a difference between them. Each fragment crosses an isoform-specific junction (which isoform) and carries a phased heterozygous site (which haplotype); the haplotype carrying the lead's ALT allele is read off the donor's phased genotype at the lead, in the same phase frame as phASER's `PW`.

## Model

- units: donor x haplotype; `y` = T1 fragments, `n` = T1 + T2 fragments
- `logit p = mu + u_donor + beta * ALT`, `u_donor ~ N(0, sigma^2)`, `y ~ BetaBinomial(n, p, rho)`; likelihood-ratio test of `beta = 0`
- `beta` = within-donor log odds of T1 on the ALT haplotype vs the REF haplotype; the random donor intercept absorbs each donor's baseline isoform ratio, and with it the donor's cell composition
- fitted when >= 10 lead-heterozygous donors have fragments on both haplotypes; BH over the fitted gate-family pairs
- built-in null: the same model on donors HOMOZYGOUS at the lead (hap 1 as a pseudo-ALT)
- sensitivity: a model-free stratified score test with empirical variance, and the between-donor Spearman of lead dosage vs the donor's T1 fraction from the same reads

## Result

| family | pairs | fitted | genes fitted | p < 0.05 | **q < 0.05** | genes q < 0.05 | sign = between-donor (all / q < 0.05) |
|---|---|---|---|---|---|---|---|
| gate_family | 3,865 | 3,361 | 481 | 898 | **590** | 185 | 0.65 / 0.87 |
| coloc_nominated | 250 | 189 | 29 | 32 | **5** | 2 | 0.60 / 1.00 |

q-values are computed within the gate family only; a coloc-nominated pair outside the gate family is fitted but carries no q. Pairs of one gene share donors and fragments, so pair-level counts are not independent; the gene counts are the conservative read.

## Calibration

- homozygous-at-lead null: 3,190 fits, fraction p < 0.05 = 0.059, lambda_GC = 1.07
- lead-heterozygous test: lambda_GC = 2.98 (inflation is expected here where the leads are real switch-QTLs)
- GLMM beta vs score-test Z, Spearman: 0.86
- median lead-to-het-site distance: 18 kb
- gate family: 35 fits did not converge; 132 hit the |beta| = 10 bound (19 of them at q < 0.05). A bound fit is quasi-separation (one haplotype carries one isoform only): its sign is informative, its magnitude is not.

## IsoGraph modules

Every switch pair belongs to a gene of an age-selected IsoGraph co-switching module. The within-donor test asks whether that gene's switch is under cis-genetic control that cell composition cannot produce; it does not explain why genes co-switch.

| module | age trait | age effect | pairs fitted | genes fitted | pairs q < 0.05 | genes q < 0.05 |
|---|---|---|---|---|---|---|
| M000 | Age_linear | +0.26 | 549 | 84 | 88 | 29 |
| M002 | Age_linear | +0.15 | 687 | 92 | 111 | 34 |
| M003 | Age_spline | +2.94 | 631 | 100 | 81 | 30 |
| M004 | Age_linear | +0.15 | 545 | 75 | 103 | 30 |
| M006 | Age_spline | +0.84 | 585 | 80 | 146 | 40 |
| M008 | Age_spline | +1.67 | 290 | 38 | 40 | 17 |
| M011 | Age_spline | +1.42 | 74 | 12 | 21 | 5 |

**Does the lead move the isoforms along the module's switch axis?** Within a gene, ALT is one fixed allele, so beta should track each pair's module polarity (r(T1) - r(T2) against the module score) if the variant and the module move the same axis. Over 417 gate-family genes with >= 3 fitted pairs, median |Spearman(beta, polarity)| = 0.40, against 0.32 with polarity permuted within gene (95th percentile 0.35; permutation p = 0.000999). Among genes with >= 2 significant pairs, 44 of 137 have every significant pair on the same side of the module axis.

`beta_along_module` in `allelic_test.parquet` is beta signed by the pair's module polarity; across genes its sign is arbitrary until oriented to a risk allele (`risk_along_module`: > 0 when the risk allele shifts the isoforms the way the module score rises).

## Caveats

- Cell composition cancels within a donor only if allelic effects do not differ by cell type.
- Some reference mapping bias survives WASP. It shifts ALT vs REF fragment totals; it moves the T1 fraction within a haplotype only if it differs between the two isoforms' junction reads.
- A statistical phase switch between the lead and a fragment's heterozygous site swaps ALT and REF for that fragment and pulls beta toward zero; the lead-to-site distance is recorded per pair.
- The lead is the all_samples switch-QTL lead, not necessarily the causal variant; a causal variant in imperfect LD with it attenuates beta.
- `beta` is oriented to the QTL ALT allele. Orienting it to the GWAS risk allele needs the rsID -> risk-allele table (`--risk-alleles`, built on Bridges-2 with `coloc_direction._gwas_risk`).

## Top gate-family pairs

| gene | T1 | T2 | lead | het donors (paired) | beta | p | q | score p | between rho | hom-null p |
|---|---|---|---|---|---|---|---|---|---|---|
| ENSG00000143337.19 | ENST00000606911.7 | ENST00000528443.6 | rs10913906 | 176 (169) | -7.69 | 1.2e-116 | 2.0e-113 | 9.0e-31 | -0.77 | 0.44 |
| ENSG00000143337.19 | ENST00000435319.8 | ENST00000528443.6 | rs10913906 | 176 (169) | -7.69 | 1.2e-116 | 2.0e-113 | 9.0e-31 | -0.77 | 0.44 |
| ENSG00000130770.18 | ENST00000335514.10 | ENST00000497986.5 | rs8559 | 195 (192) | +6.19 | 6.1e-95 | 6.8e-92 | 1.4e-27 | +0.77 | 0.30 |
| ENSG00000143337.19 | ENST00000606911.7 | ENST00000271583.7 | rs10913906 | 179 (173) | -5.16 | 1.3e-88 | 8.7e-86 | 2.2e-30 | -0.78 | 0.26 |
| ENSG00000143337.19 | ENST00000435319.8 | ENST00000271583.7 | rs10913906 | 179 (173) | -5.16 | 1.3e-88 | 8.7e-86 | 2.2e-30 | -0.78 | 0.26 |
| ENSG00000291131.2 | ENST00000795213.1 | ENST00000795217.1 | rs17272903 | 171 (138) | -4.94 | 1.3e-76 | 7.1e-74 | 7.0e-21 | -0.77 | 0.79 |
| ENSG00000143337.19 | ENST00000528443.6 | ENST00000474875.5 | rs10913906 | 177 (174) | +2.59 | 1.4e-72 | 6.8e-70 | 8.9e-27 | +0.76 | 0.79 |
| ENSG00000291131.2 | ENST00000795217.1 | ENST00000795219.1 | rs17272903 | 186 (142) | +3.36 | 8.3e-63 | 3.5e-60 | 7.3e-21 | +0.77 | 0.41 |
| ENSG00000291131.2 | ENST00000795213.1 | ENST00000597550.6 | rs17272903 | 168 (129) | -4.77 | 7.5e-61 | 2.8e-58 | 8.7e-19 | -0.52 | 1.00 |
| ENSG00000147408.16 | ENST00000692225.2 | ENST00000397998.7 | rs10888163 | 191 (178) | -2.96 | 4.2e-49 | 1.4e-46 | 5.7e-21 | -0.87 | 0.81 |
| ENSG00000147408.16 | ENST00000695896.1 | ENST00000397998.7 | rs10888163 | 190 (175) | -2.96 | 1.6e-48 | 4.4e-46 | 2.3e-20 | -0.87 | 0.88 |
| ENSG00000147408.16 | ENST00000522854.5 | ENST00000397998.7 | rs10888163 | 190 (175) | -2.96 | 1.6e-48 | 4.4e-46 | 2.3e-20 | -0.87 | 0.82 |
| ENSG00000291131.2 | ENST00000597550.6 | ENST00000795219.1 | rs17272903 | 186 (136) | +2.88 | 3.5e-48 | 9.1e-46 | 9.7e-19 | +0.56 | 0.42 |
| ENSG00000188186.11 | ENST00000468582.5 | ENST00000460732.5 | rs12878 | 234 (233) | -0.54 | 2.6e-44 | 6.3e-42 | 2.1e-26 | -0.34 | 0.39 |
| ENSG00000188186.11 | ENST00000468582.5 | ENST00000441173.1 | rs12878 | 234 (233) | -0.53 | 8.4e-44 | 1.9e-41 | 2.8e-26 | -0.35 | 0.38 |
| ENSG00000268362.6 | ENST00000668060.1 | ENST00000654372.1 | rs62117569 | 136 (135) | +2.61 | 2.0e-42 | 4.3e-40 | 7.8e-18 | +0.50 | — |
| ENSG00000268362.6 | ENST00000668060.1 | ENST00000594934.6 | rs62117569 | 136 (135) | +2.56 | 8.0e-42 | 1.6e-39 | 9.8e-18 | +0.50 | — |
| ENSG00000268362.6 | ENST00000659757.1 | ENST00000594934.6 | rs62117569 | 136 (135) | +2.52 | 5.4e-41 | 1.0e-38 | 1.2e-17 | +0.49 | 0.98 |
| ENSG00000196557.14 | ENST00000348261.11 | ENST00000562079.6 | rs75197420 | 63 (63) | -2.16 | 3.0e-38 | 5.3e-36 | 9.1e-12 | -0.49 | 0.15 |
| ENSG00000196557.14 | ENST00000562079.6 | ENST00000565831.7 | rs75197420 | 63 (63) | +2.16 | 5.5e-38 | 9.2e-36 | 9.0e-12 | +0.52 | 0.11 |
| ENSG00000106591.4 | ENST00000223324.3 | ENST00000413995.1 | rs631132 | 131 (127) | -2.47 | 1.5e-37 | 2.4e-35 | 1.0e-13 | -0.42 | 0.47 |
| ENSG00000268362.6 | ENST00000659757.1 | ENST00000654372.1 | rs62117569 | 136 (135) | +2.24 | 2.7e-37 | 4.1e-35 | 2.2e-17 | +0.49 | 0.38 |
| ENSG00000196557.14 | ENST00000348261.11 | ENST00000639478.1 | rs75197420 | 63 (63) | -1.92 | 4.3e-36 | 6.3e-34 | 1.3e-11 | -0.44 | 0.08 |
| ENSG00000291131.2 | ENST00000597550.6 | ENST00000795217.1 | rs17272903 | 122 (34) | -6.51 | 4.6e-36 | 6.5e-34 | 9.6e-07 | -0.85 | 0.92 |
| ENSG00000196557.14 | ENST00000639478.1 | ENST00000565831.7 | rs75197420 | 63 (63) | +1.90 | 7.0e-36 | 9.4e-34 | 1.3e-11 | +0.34 | 0.04 |

