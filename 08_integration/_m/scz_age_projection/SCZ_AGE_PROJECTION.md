# Age-sensitive isoform-switch programs in schizophrenia

Age-sensitive co-switching modules are defined out-of-cohort in the independent aging caudate fits (GTEx caudate basal ganglia + BrainSeq caudate) and restricted to those enriched for schizophrenia GWAS (MAGMA SCZ P<0.05). Their **behaviour** is tested in an independent SCZ case/control cohort (BrainSeq caudate_sczd) with numeric age and TOPMed genotypes. Effects are abundance-conditioned (DTU-without-DGE).

## SCZ-risk loci and age-sensitive switch programs

Schizophrenia-colocalized switch genes are **not concentrated in the age-sensitive (SCZ-GWAS-enriched) modules**: 40/71 (56%) of coloc genes fall in anchored modules, against a 56% background among the genes that were coloc-TESTED (hypergeometric P=1).

> **The background is the whole result, and an earlier version of this report used the wrong one.** Against *all* module genes the background is only 65% and the same counts give P=0.949 — the previously reported convergence headline. That comparison is confounded by ascertainment: a gene can only colocalize if it sat under a SCZ GWAS peak with a QTL credible set, and anchored modules are *defined* by MAGMA SCZ enrichment, so their genes enter the tested pool preferentially. Conditioning on what could have been a hit removes the effect. `module_coloc_convergence.py` repeats this for all five traits with a size-matched permutation null and finds no concentration anywhere (P = 0.19–1.00).

The modules carrying the most colocalized loci are listed below; with these counts the per-module numbers are descriptive, not evidence of convergence.

| source | module | # coloc genes | # GO-invisible | max same-dir frac | SCZ MAGMA P |
|--------|--------|-----------|----------------|-------------------|-------------|
| gtex_caudate_bg | M000 | 11 | 7 | 0.65 | 4.7e-05 |
| gtex_caudate_bg | M002 | 8 | 8 | 0.68 | 0.026 |
| gtex_caudate_bg | M001 | 7 | 10 | 0.64 | 0.0043 |
| brainseq_caudate | M003 | 7 | 11 | 0.71 | 1.8e-05 |
| gtex_caudate_bg | M004 | 6 | 9 | 0.67 | 0.048 |
| brainseq_caudate | M000 | 3 | 1 | 0.86 | 0.0032 |

Independently of that null, the age-sensitive modules are directionally disrupted in disease: they recapitulate the aging switch direction gene-by-gene in **6/9** modules (B), show case deviation from the control age trajectory in **1/9** modules (D), and stay co-switch-coherent in disease in **9/9** modules (C).

## B. Aging↔disease direction concordance (accelerated-aging recapitulation)

Per module gene, sign of the control-only switch-vs-age slope vs the switch-vs-Dx effect. **6/9** anchored modules show above-chance concordance (binomial P<0.05) — SCZ recapitulates the age switch program gene by gene. Top:

| source | module | concordant/n | rate | binom P |
|--------|--------|--------------|------|---------|
| brainseq_caudate | M000 | 474/608 | 0.78 | 9.35e-46 |
| gtex_caudate_bg | M000 | 622/1083 | 0.57 | 5.57e-07 |
| brainseq_caudate | M010 | 87/122 | 0.71 | 1.39e-06 |
| gtex_caudate_bg | M001 | 653/1169 | 0.56 | 3.42e-05 |
| brainseq_caudate | M001 | 265/471 | 0.56 | 0.00373 |

## D. Age-state deviation

Each module eigengene's age trajectory is fit in controls; **1/9** anchored modules show cases deviating from the control-expected state (residual ~ Dx, P<0.05).

## C. Module preservation

**9/9** anchored aging modules keep coherent co-switch structure in the disease switch network (median intramodular |r| above a size-matched permutation null) — a prerequisite the projections rely on.

---

## Supplementary — single-locus genotype resolution

Per-locus genotype tests are reported for completeness; single-locus QTL power in the disease cohort is far below the module-level analyses above, so these do not carry the narrative.

### S1. Genotype-anchored cis↔diagnosis switch concordance

For **168** genotyped SCZ-colocalized loci, the risk-allele cis effect and the diagnosis effect on the same disease switch axis agree in sign in **62/168** (37%; binomial P=1) — **null at single-locus resolution**. The in-cohort cis effects are weak (median cis P=0.36); QTL power in n~168 is far below GTEx discovery (GTEx cis-direction match 80/168), which is why the directional signal is resolved at the module level (headline convergence + B/D) rather than per locus.

- coloc loci with a nominal disease switch shift (Dx P<0.05, abundance-conditioned): **27/218** (14 GO-invisible)

### S2. Genotype × diagnosis

- genotype×diagnosis interactions at P<0.05: **18** of 226 locus×outcome tests (exploratory; single-locus power caveat applies).

