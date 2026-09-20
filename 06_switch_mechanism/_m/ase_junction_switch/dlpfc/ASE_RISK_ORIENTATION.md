# Risk-allele orientation of the allelic switch test — dlpfc

`beta` from the allelic test is the within-donor log odds of T1 on the haplotype carrying the switch-QTL lead's ALT allele. ALT is arbitrary, so this stage re-signs it to the allele that raises disease risk. The risk allele is the trait-increasing allele at the GWAS lead of the locus where the gene colocalizes (sQTL PP4 >= 0.8); the flip follows the signed LD between that variant and the fitted lead, in the BrainSEQ genotypes the QTL mapping used, gated at |r| >= 0.8.

`risk_along_module` > 0 means the risk allele shifts the isoforms the way the module score rises.

## Coverage

- nominated pair x trait rows: **97** over 13 genes (69 fitted)
- oriented: **0** (0 genes)

| orientation status | rows |
|---|---:|
| r_gate | 72 |
| no_ld | 24 |
| palindromic_gwas_lead | 1 |

## LD between the fitted lead and the disease lead

This is what limits the arm. The switch-QTL lead the test is fitted on and the GWAS lead of the locus are usually not on the same haplotype, so the disease allele cannot be placed on the fitted lead at all.

- distinct lead x disease-lead variant pairs with LD: **10**
- median |r|: **0.080**
- at |r| >= 0.8: **1**

## Ancestry

BrainSEQ is roughly half African-American and the GWAS are European, so a pooled fit and a pooled LD estimate both mix ancestries. `beta_ea` (EA donors only) is what gets oriented, against EA-panel LD; the AA arm is fitted beside it and never pooled.

- donors: EA 208, AA 226, neither 64
- median |r| lead-to-risk-allele: EA 0.080, AA 0.076
- EA and AA betas agree in sign in 10 of 19 pairs fitted in both

## Oriented pairs at q < 0.05 (EA refit, BH within this family)

- 0 pairs over 0 genes
- risk allele raises the module score in 0, lowers it in 0

## Caveats

- The locus GWAS lead is not necessarily the causal variant; orientation inherits that assumption, and imperfect LD with the fitted lead attenuates `beta` itself.
- Every caveat of the allelic test carries over (cell-type-specific allelic effects, residual mapping bias, phase switches).
- Pairs left unoriented are not negative results: the disease allele could not be placed on a haplotype with the fitted lead.
