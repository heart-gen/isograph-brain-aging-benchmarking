# Allele-specific switch test — hippocampus

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
| gate_family | 2,213 | 1,892 | 291 | 448 | **247** | 92 | 0.63 / 0.86 |
| coloc_nominated | 202 | 114 | 17 | 10 | **6** | 2 | 0.60 / 1.00 |

q-values are computed within the gate family only; a coloc-nominated pair outside the gate family is fitted but carries no q. Pairs of one gene share donors and fragments, so pair-level counts are not independent; the gene counts are the conservative read.

## Calibration

- homozygous-at-lead null: 1,802 fits, fraction p < 0.05 = 0.056, lambda_GC = 1.03
- lead-heterozygous test: lambda_GC = 2.66 (inflation is expected here where the leads are real switch-QTLs)
- GLMM beta vs score-test Z, Spearman: 0.87
- median lead-to-het-site distance: 22 kb
- gate family: 48 fits did not converge; 61 hit the |beta| = 10 bound (5 of them at q < 0.05). A bound fit is quasi-separation (one haplotype carries one isoform only): its sign is informative, its magnitude is not.

## IsoGraph modules

Every switch pair belongs to a gene of an age-selected IsoGraph co-switching module. The within-donor test asks whether that gene's switch is under cis-genetic control that cell composition cannot produce; it does not explain why genes co-switch.

| module | age trait | age effect | pairs fitted | genes fitted | pairs q < 0.05 | genes q < 0.05 |
|---|---|---|---|---|---|---|
| M000 | Age_linear | +0.27 | 193 | 28 | 36 | 12 |
| M001 | Age_linear | +0.44 | 295 | 44 | 39 | 16 |
| M002 | Age_linear | +0.35 | 195 | 28 | 17 | 6 |
| M003 | Age_linear | +0.41 | 90 | 13 | 2 | 2 |
| M004 | Age_linear | +0.42 | 451 | 77 | 57 | 22 |
| M005 | Age_spline | -3.65 | 80 | 12 | 13 | 3 |
| M006 | Age_linear | +0.17 | 33 | 7 | 2 | 1 |
| M007 | Age_linear | +0.35 | 203 | 31 | 13 | 6 |
| M008 | Age_linear | +0.39 | 155 | 22 | 31 | 9 |
| M009 | Age_linear | +0.30 | 99 | 12 | 9 | 5 |
| M010 | Age_linear | +0.28 | 11 | 3 | 8 | 3 |
| M011 | Age_linear | +0.39 | 61 | 8 | 19 | 6 |
| M012 | Age_linear | +0.28 | 16 | 4 | 0 | 0 |
| M013 | Age_linear | +0.25 | 3 | 1 | 0 | 0 |
| M015 | Age_linear | +0.45 | 7 | 1 | 1 | 1 |

**Does the lead move the isoforms along the module's switch axis?** Within a gene, ALT is one fixed allele, so beta should track each pair's module polarity (r(T1) - r(T2) against the module score) if the variant and the module move the same axis. Over 239 gate-family genes with >= 3 fitted pairs, median |Spearman(beta, polarity)| = 0.42, against 0.33 with polarity permuted within gene (95th percentile 0.37; permutation p = 0.000999). Among genes with >= 2 significant pairs, 25 of 64 have every significant pair on the same side of the module axis.

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
| ENSG00000133943.21 | ENST00000523816.5 | ENST00000412671.6 | rs4900071 | 134 (123) | +10.00 (bound) | 7.9e-80 | 1.5e-76 | 4.9e-23 | +0.82 | 0.94 |
| ENSG00000133943.21 | ENST00000428926.6 | ENST00000412671.6 | rs4900071 | 134 (123) | +8.39 | 1.7e-79 | 1.6e-76 | 3.5e-23 | +0.84 | 0.12 |
| ENSG00000133943.21 | ENST00000256324.15 | ENST00000428926.6 | rs4900071 | 134 (126) | -6.35 | 2.6e-74 | 1.7e-71 | 5.0e-23 | -0.79 | 0.84 |
| ENSG00000133943.21 | ENST00000256324.15 | ENST00000523816.5 | rs4900071 | 134 (126) | -5.79 | 1.5e-71 | 7.0e-69 | 1.3e-22 | -0.69 | 0.01 |
| ENSG00000133943.21 | ENST00000523816.5 | ENST00000523461.5 | rs4900071 | 134 (128) | +4.79 | 9.5e-71 | 3.6e-68 | 1.2e-22 | +0.79 | 0.03 |
| ENSG00000133943.21 | ENST00000428926.6 | ENST00000523461.5 | rs4900071 | 134 (128) | +4.81 | 1.3e-70 | 4.1e-68 | 1.3e-22 | +0.79 | 0.02 |
| ENSG00000130770.18 | ENST00000335514.10 | ENST00000497986.5 | rs8559 | 179 (165) | +5.61 | 5.5e-62 | 1.5e-59 | 1.2e-20 | +0.73 | 0.16 |
| ENSG00000133943.21 | ENST00000412671.6 | ENST00000520328.5 | rs4900071 | 134 (129) | -2.80 | 3.9e-46 | 9.1e-44 | 2.2e-22 | -0.77 | 0.45 |
| ENSG00000133943.21 | ENST00000256324.15 | ENST00000520328.5 | rs4900071 | 134 (131) | -2.61 | 4.4e-44 | 9.3e-42 | 4.9e-22 | -0.71 | 0.95 |
| ENSG00000103494.16 | ENST00000621565.5 | ENST00000262135.9 | rs9925121 | 74 (64) | -2.85 | 3.6e-38 | 6.8e-36 | 9.1e-13 | -0.83 | 0.67 |
| ENSG00000103494.16 | ENST00000563746.5 | ENST00000621565.5 | rs9925121 | 78 (64) | +2.71 | 1.6e-36 | 2.7e-34 | 9.2e-13 | +0.80 | 0.58 |
| ENSG00000161551.15 | ENST00000640955.1 | ENST00000638827.1 | rs2191669 | 211 (202) | -1.52 | 1.5e-34 | 2.4e-32 | 9.2e-20 | -0.53 | 0.62 |
| ENSG00000133943.21 | ENST00000520328.5 | ENST00000523461.5 | rs4900071 | 134 (130) | +1.91 | 1.0e-29 | 1.5e-27 | 7.1e-20 | +0.62 | 0.70 |
| ENSG00000106591.4 | ENST00000223324.3 | ENST00000413995.1 | rs3757560 | 120 (118) | -2.33 | 1.6e-29 | 2.2e-27 | 3.6e-14 | -0.45 | — |
| ENSG00000103494.16 | ENST00000564374.5 | ENST00000262135.9 | rs9925121 | 108 (87) | -2.09 | 7.5e-28 | 9.5e-26 | 1.7e-12 | -0.81 | 0.86 |
| ENSG00000103494.16 | ENST00000563746.5 | ENST00000564374.5 | rs9925121 | 109 (87) | +1.99 | 7.9e-27 | 9.3e-25 | 2.2e-12 | +0.78 | 0.88 |
| ENSG00000106591.4 | ENST00000223324.3 | ENST00000432845.1 | rs3757560 | 120 (118) | -2.30 | 4.3e-26 | 4.8e-24 | 3.4e-14 | -0.42 | — |
| ENSG00000251600.10 | ENST00000793869.1 | ENST00000793864.1 | rs66570349 | 90 (36) | +3.67 | 3.4e-23 | 3.5e-21 | 2.3e-05 | +0.69 | 0.94 |
| ENSG00000103494.16 | ENST00000565343.2 | ENST00000262135.9 | rs9925121 | 104 (77) | -1.88 | 4.5e-23 | 4.5e-21 | 1.6e-12 | -0.77 | 0.91 |
| ENSG00000103494.16 | ENST00000563746.5 | ENST00000565343.2 | rs9925121 | 101 (77) | +1.92 | 5.4e-23 | 5.1e-21 | 1.9e-12 | +0.78 | 0.87 |
| ENSG00000170074.21 | ENST00000697110.1 | ENST00000697109.1 | rs13169820 | 115 (70) | -2.78 | 1.9e-20 | 1.7e-18 | 1.8e-08 | -0.55 | 0.99 |
| ENSG00000170074.21 | ENST00000697112.1 | ENST00000697109.1 | rs13169820 | 115 (70) | -2.78 | 2.0e-20 | 1.7e-18 | 2.0e-08 | -0.49 | 0.97 |
| ENSG00000170074.21 | ENST00000697103.1 | ENST00000697109.1 | rs13169820 | 115 (70) | -2.78 | 2.1e-20 | 1.7e-18 | 2.3e-08 | -0.54 | 0.98 |
| ENSG00000170074.21 | ENST00000503845.6 | ENST00000697109.1 | rs13169820 | 115 (70) | -2.78 | 2.1e-20 | 1.7e-18 | 2.3e-08 | -0.54 | 0.98 |
| ENSG00000181378.14 | ENST00000341552.10 | ENST00000295729.6 | rs6736922 | 169 (95) | -2.90 | 2.8e-20 | 2.1e-18 | 3.4e-09 | -0.34 | 0.58 |

