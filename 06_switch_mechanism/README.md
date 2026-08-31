# 06 — Switch mechanism

**Question:** are the switches real measurements rather than model artifacts, what do
they do to the protein, and do the results survive sensitivity analysis?

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–02 | `switch_consequence`, `switch_consequence_meta` | Structural consequence of the switch axis under a within-gene permutation null; cross-region meta |
| 03–04 | `validate_switch_splicing_{brainseq,gtex}` | Orthogonal PSI / junction validation of switch pairs |
| 05 | `switch_orthogonal_confirm` | Anchored + global-null confirmation against the compositional-closure baseline |
| 06 | `isa_concordance` | Concordance with satuRn / ISA differential transcript usage |
| 07 | `longread_switch_confirm` | ONT DLPFC long-read confirmation (Aguzzoli-Heberle 2024) |
| 08–09 | `download_clinical`, `clinical_consequence` | gnomAD LOEUF constraint + ClinVar pathogenic density |
| 10–11 | `scz_confound_sensitivity`, `switch_feature_sensitivity` | Medication/toxicology/smoking availability audit; five feature-construction axes |

## Key results

- Consequence is **productive remodeling, not decay**: of 9 structural classes only
  UTR-remodeled (1.27×, 10/10 regions) and CDS-remodeled (1.04×, 10/10) are enriched;
  NMD routing, biotype switch and coding-status loss are depleted. Identical in
  GO-invisible and GO-visible modules.
- Switch genes are more LoF-constrained than genome-wide (median LOEUF 0.72 vs 0.94),
  while switched exons carry *lower* ClinVar P/LP density (ratio 0.18) — expected
  alternative-exon biology.

## Two sensitivity results that must be reported as stated

- **rRNA rate is not a confounder.** It has no marginal association with diagnosis
  (d = −0.009, p = 0.93); its association appears only after conditioning on the
  published covariates (partial r = +0.108), which is the collider signature. The
  modules *most* correlated with it survive and the least correlated fall — the opposite
  of technical-artifact removal. Do not add it to the inference model on confounding
  grounds. Report as a sensitivity.
- **Quantification is the striking axis.** Per-gene switch-age effects are essentially
  uncorrelated between Salmon and RSEM on matched regions (Pearson 0.007 caudate,
  −0.002 hippocampus, sign concordance 0.499). This confounds quantifier with cohort, so
  it is an upper bound, not an isolated estimate.
- The feature-sensitivity harness holds the **module partition fixed**, so it measures
  the stability of the representation and its trait signal, not of an independently
  refit network.

**CLIs:** `isograph_benchmark/real_data/{switch_consequence,switch_consequence_meta,validate_switch_splicing,switch_orthogonal_confirm,isa_concordance,longread_switch_confirm,clinical_consequence,clinical_consequence_meta,scz_confound_sensitivity,switch_feature_sensitivity}.py`.

## Display items

S-real-5 `figSwitchConsequence`, S-real-7 `figClinicalConsequence`.
