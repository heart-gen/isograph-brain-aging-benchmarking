# Allele-specific switch test — dlpfc

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
| gate_family | 1,014 | 853 | 115 | 189 | **98** | 37 | 0.59 / 0.82 |
| coloc_nominated | 97 | 69 | 9 | 7 | **0** | 0 | 0.55 / — |

q-values are computed within the gate family only; a coloc-nominated pair outside the gate family is fitted but carries no q. Pairs of one gene share donors and fragments, so pair-level counts are not independent; the gene counts are the conservative read.

## Calibration

- homozygous-at-lead null: 851 fits, fraction p < 0.05 = 0.043, lambda_GC = 1.01
- lead-heterozygous test: lambda_GC = 2.61 (inflation is expected here where the leads are real switch-QTLs)
- GLMM beta vs score-test Z, Spearman: 0.87
- median lead-to-het-site distance: 18 kb
- gate family: 13 fits did not converge; 21 hit the |beta| = 10 bound (6 of them at q < 0.05). A bound fit is quasi-separation (one haplotype carries one isoform only): its sign is informative, its magnitude is not.

## IsoGraph modules

Every switch pair belongs to a gene of an age-selected IsoGraph co-switching module. The within-donor test asks whether that gene's switch is under cis-genetic control that cell composition cannot produce; it does not explain why genes co-switch.

| module | age trait | age effect | pairs fitted | genes fitted | pairs q < 0.05 | genes q < 0.05 |
|---|---|---|---|---|---|---|
| M000 | Age_spline | -2.05 | 441 | 58 | 47 | 20 |
| M007 | Age_spline | -1.51 | 278 | 32 | 35 | 11 |
| M008 | Age_linear | +0.23 | 134 | 25 | 16 | 6 |

**Does the lead move the isoforms along the module's switch axis?** Within a gene, ALT is one fixed allele, so beta should track each pair's module polarity (r(T1) - r(T2) against the module score) if the variant and the module move the same axis. Over 104 gate-family genes with >= 3 fitted pairs, median |Spearman(beta, polarity)| = 0.40, against 0.31 with polarity permuted within gene (95th percentile 0.37; permutation p = 0.00599). Among genes with >= 2 significant pairs, 7 of 22 have every significant pair on the same side of the module axis.

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
| ENSG00000206530.11 | ENST00000393845.9 | ENST00000461734.1 | rs2936812 | 194 (159) | -10.00 (bound) | 3.0e-28 | 2.6e-25 | 1.5e-13 | -0.51 | 0.69 |
| ENSG00000157103.13 | ENST00000287766.10 | ENST00000643396.1 | rs2930154 | 71 (58) | -3.06 | 1.3e-22 | 5.5e-20 | 5.4e-08 | -0.46 | 0.08 |
| ENSG00000224271.9 | ENST00000651403.1 | ENST00000652399.1 | rs6358 | 137 (130) | +0.94 | 1.5e-20 | 3.5e-18 | 1.0e-12 | +0.49 | 0.12 |
| ENSG00000224271.9 | ENST00000651403.1 | ENST00000651662.1 | rs6358 | 137 (130) | +0.94 | 1.6e-20 | 3.5e-18 | 1.0e-12 | +0.49 | 0.12 |
| ENSG00000150753.12 | ENST00000515390.5 | ENST00000515676.5 | rs2548546 | 177 (177) | +0.39 | 2.9e-20 | 4.2e-18 | 2.5e-15 | +0.09 | 0.40 |
| ENSG00000150753.12 | ENST00000503026.5 | ENST00000515390.5 | rs2548546 | 177 (177) | -0.39 | 3.0e-20 | 4.2e-18 | 2.6e-15 | -0.09 | 0.39 |
| ENSG00000147679.12 | ENST00000521071.1 | ENST00000521974.1 | rs979867 | 149 (85) | +2.28 | 1.2e-19 | 1.4e-17 | 1.0e-09 | +0.48 | 0.28 |
| ENSG00000147679.12 | ENST00000309822.7 | ENST00000521071.1 | rs979867 | 149 (86) | -2.24 | 1.7e-19 | 1.6e-17 | 1.3e-09 | -0.51 | 0.28 |
| ENSG00000150995.21 | ENST00000649015.2 | ENST00000650294.1 | rs7613447 | 193 (187) | -0.29 | 1.9e-19 | 1.6e-17 | 6.3e-15 | -0.29 | — |
| ENSG00000150995.21 | ENST00000650294.1 | ENST00000443694.5 | rs7613447 | 193 (186) | +0.29 | 2.0e-19 | 1.6e-17 | 6.1e-15 | +0.40 | — |
| ENSG00000150995.21 | ENST00000648266.1 | ENST00000354582.12 | rs7613447 | 193 (186) | -0.29 | 2.0e-19 | 1.6e-17 | 6.1e-15 | -0.40 | — |
| ENSG00000135778.12 | ENST00000490807.5 | ENST00000494689.5 | rs12743233 | 19 (19) | -2.60 | 2.2e-17 | 1.5e-15 | 1.2e-04 | -0.56 | 0.41 |
| ENSG00000135778.12 | ENST00000494689.5 | ENST00000366627.4 | rs12743233 | 19 (18) | +2.77 | 5.3e-17 | 3.4e-15 | 1.3e-04 | +0.58 | 0.55 |
| ENSG00000135778.12 | ENST00000366628.10 | ENST00000494689.5 | rs12743233 | 19 (18) | -2.77 | 6.6e-17 | 4.0e-15 | 1.2e-04 | -0.61 | 0.96 |
| ENSG00000150995.21 | ENST00000649015.2 | ENST00000354582.12 | rs7613447 | 193 (186) | -0.25 | 7.2e-15 | 4.1e-13 | 7.1e-13 | -0.25 | 0.28 |
| ENSG00000150995.21 | ENST00000354582.12 | ENST00000443694.5 | rs7613447 | 193 (187) | +0.25 | 7.8e-15 | 4.2e-13 | 7.0e-13 | +0.10 | 0.31 |
| ENSG00000150995.21 | ENST00000648266.1 | ENST00000650294.1 | rs7613447 | 193 (187) | -0.23 | 3.7e-13 | 1.8e-11 | 8.0e-12 | -0.15 | 0.28 |
| ENSG00000157103.13 | ENST00000643396.1 | ENST00000646487.1 | rs2930154 | 71 (58) | +1.77 | 1.7e-12 | 8.0e-11 | 3.2e-06 | +0.11 | 0.13 |
| ENSG00000116688.18 | ENST00000444836.5 | ENST00000675959.1 | rs12567779 | 169 (156) | -0.83 | 4.1e-11 | 1.8e-09 | 3.5e-08 | -0.39 | 0.25 |
| ENSG00000116688.18 | ENST00000674817.1 | ENST00000675959.1 | rs12567779 | 169 (156) | -0.83 | 4.1e-11 | 1.8e-09 | 3.5e-08 | -0.39 | 0.25 |
| ENSG00000116688.18 | ENST00000444836.5 | ENST00000235329.10 | rs12567779 | 169 (156) | -0.83 | 5.0e-11 | 1.9e-09 | 4.1e-08 | -0.53 | 0.20 |
| ENSG00000116688.18 | ENST00000675872.1 | ENST00000675959.1 | rs12567779 | 169 (156) | -0.83 | 5.0e-11 | 1.9e-09 | 4.1e-08 | -0.53 | 0.20 |
| ENSG00000116688.18 | ENST00000235329.10 | ENST00000674817.1 | rs12567779 | 169 (156) | +0.83 | 5.0e-11 | 1.9e-09 | 4.1e-08 | +0.53 | 0.21 |
| ENSG00000116688.18 | ENST00000235329.10 | ENST00000675872.1 | rs12567779 | 169 (156) | +0.82 | 5.5e-11 | 2.0e-09 | 4.3e-08 | +0.51 | 0.19 |
| ENSG00000150753.12 | ENST00000503454.5 | ENST00000515676.5 | rs2548546 | 177 (177) | +0.30 | 1.8e-10 | 6.1e-09 | 2.5e-09 | -0.00 | 0.78 |

