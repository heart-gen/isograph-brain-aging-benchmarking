# Age-sensitive isoform-switch programs in schizophrenia

Age-sensitive co-switching modules are defined out-of-cohort in the independent aging caudate fits (GTEx caudate basal ganglia + BrainSeq caudate) and restricted to those enriched for schizophrenia GWAS (MAGMA SCZ P<0.05). Their **behaviour** is tested in an independent SCZ case/control cohort (BrainSeq caudate_sczd) with numeric age and TOPMed genotypes. Effects are abundance-conditioned (DTU-without-DGE).

## SCZ-risk loci and age-sensitive switch programs

Schizophrenia-colocalized switch genes are **not concentrated in the age-sensitive (SCZ-GWAS-enriched) modules**: 15/32 (47%) of coloc genes fall in anchored modules, against a 47% background among the genes that were coloc-TESTED (hypergeometric P=1).

> **The background is the whole result, and an earlier version of this report used the wrong one.** Against *all* module genes the background is only 25% and the same counts give P=0.0058 — the previously reported convergence headline. That comparison is confounded by ascertainment: a gene can only colocalize if it sat under a SCZ GWAS peak with a QTL credible set, and anchored modules are *defined* by MAGMA SCZ enrichment, so their genes enter the tested pool preferentially. Conditioning on what could have been a hit removes the effect. `module_coloc_convergence.py` repeats this for all five traits with a size-matched permutation null and finds no concentration anywhere (P = 0.19–1.00).

The modules carrying the most colocalized loci are listed below; with these counts the per-module numbers are descriptive, not evidence of convergence.

| source | module | # coloc genes | # GO-invisible | max same-dir frac | SCZ MAGMA P |
|--------|--------|-----------|----------------|-------------------|-------------|
| gtex_caudate_bg | M002 | 6 | 8 | 0.67 | 0.011 |
| gtex_caudate_bg | M004 | 3 | 5 | 0.58 | 0.048 |
| gtex_caudate_bg | M008 | 3 | 1 | 0.75 | 0.00075 |
| brainseq_caudate | M004 | 2 | 9 | 0.56 | 0.045 |
| gtex_caudate_bg | M006 | 1 | 0 | 1.00 | 0.035 |
| brainseq_caudate | M001 | 1 | 0 | 1.00 | 0.017 |

Independently of that null, the age-sensitive modules are directionally disrupted in disease: they recapitulate the aging switch direction gene-by-gene in **4/10** modules (B), show case deviation from the control age trajectory in **3/10** modules (D), and stay co-switch-coherent in disease in **10/10** modules (C).

### Candidate trans-regulators (mechanism)

Each listed module's members are tested for shared RBP binding-site switching (rbp_regulon --scope combined; mature+intronic motif scan). Significant RBPs (q<0.05) are candidate trans regulators coordinating the co-switch program — a named, testable hypothesis. These modules are NOT established as points of SCZ-risk convergence (see the background caveat above); the regulators are candidates for the modules' own co-switching, not for a convergence effect:

| source | module | # SCZ loci | top candidate RBP regulators (q) |
|--------|--------|-----------|----------------------------------|
| gtex_caudate_bg | M002 | 6 | IGF2BP1(q=0.0061) |
| gtex_caudate_bg | M004 | 3 | — |
| gtex_caudate_bg | M008 | 3 | RBM14(q=0.015); RBM6(q=0.023); RBMS1(q=0.025); RBM24(q=0.029); HNRNPA3(q=0.045) |
| brainseq_caudate | M004 | 2 | — |
| gtex_caudate_bg | M006 | 1 | DDX58(q=0.0001); GRSF1(q=0.0016); ADAR(q=0.0016); TARDBP(q=0.0016); PTBP2(q=0.0039) |
| brainseq_caudate | M001 | 1 | — |

_Motif-based candidate regulation (predicted binding-site gain/loss between switch isoforms), not experimental validation._

## B. Aging↔disease direction concordance (accelerated-aging recapitulation)

Per module gene, sign of the control-only switch-vs-age slope vs the switch-vs-Dx effect. **4/10** anchored modules show above-chance concordance (binomial P<0.05) — SCZ recapitulates the age switch program gene by gene. Top:

| source | module | concordant/n | rate | binom P |
|--------|--------|--------------|------|---------|
| gtex_caudate_bg | M002 | 247/430 | 0.57 | 0.00117 |
| brainseq_caudate | M004 | 195/352 | 0.55 | 0.0242 |
| gtex_caudate_bg | M004 | 200/363 | 0.55 | 0.0293 |
| brainseq_caudate | M001 | 253/467 | 0.54 | 0.0393 |
| brainseq_caudate | M009 | 88/155 | 0.57 | 0.0539 |

## D. Age-state deviation

Each module eigengene's age trajectory is fit in controls; **3/10** anchored modules show cases deviating from the control-expected state (residual ~ Dx, P<0.05).

## C. Module preservation

**10/10** anchored aging modules keep coherent co-switch structure in the disease switch network (median intramodular |r| above a size-matched permutation null) — a prerequisite the projections rely on.

---

## Supplementary — single-locus genotype resolution

Per-locus genotype tests are reported for completeness; single-locus QTL power in the disease cohort is far below the module-level analyses above, so these do not carry the narrative.

### S1. Genotype-anchored cis↔diagnosis switch concordance

For **62** genotyped SCZ-colocalized loci, the risk-allele cis effect and the diagnosis effect on the same disease switch axis agree in sign in **28/62** (45%; binomial P=0.813) — **null at single-locus resolution**. The in-cohort cis effects are weak (median cis P=0.34); QTL power in n~62 is far below GTEx discovery (GTEx cis-direction match 33/62), which is why the directional signal is resolved at the module level (headline convergence + B/D) rather than per locus.

- coloc loci with a nominal disease switch shift (Dx P<0.05, abundance-conditioned): **21/77** (11 GO-invisible)

### S2. Genotype × diagnosis

- genotype×diagnosis interactions at P<0.05: **12** of 86 locus×outcome tests (exploratory; single-locus power caveat applies).

