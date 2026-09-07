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
| 12 | `junction_coloc_confirm` | Short-read BrainSEQ junction confirmation of the anchored SNCA/CTSH switches |

## Key results

- Consequence is **productive remodeling, not decay**: of 9 structural classes only
  UTR-remodeled (1.27×, 10/10 regions) and CDS-remodeled (1.04×, 10/10) are enriched;
  NMD routing, biotype switch and coding-status loss are depleted. Identical in
  GO-invisible and GO-visible modules.
- Switch genes are more LoF-constrained than genome-wide (median LOEUF 0.72 vs 0.94),
  while switched exons carry *lower* ClinVar P/LP density (ratio 0.18) — expected
  alternative-exon biology.
- **SNCA's anchored switch is confirmed in short read; CTSH's is not** (step 12). On the
  exact contrast Fig 4A draws — the anchored proximal first exon against the canonical
  distal one — minor-form usage is **0.189 in DLPFC (n = 222)** and **0.234 in caudate
  (n = 238)**, far above the pre-registered 0.05 threshold. The ONT long-read check (step
  07) had put the anchored isoform at **0.29%** and failed it; short read measures the
  junction the sQTL actually tags rather than a whole-transcript proxy, at ~20x the n, so
  the long-read result is an assay limitation rather than a refutation. **CTSH reaches only
  0.016–0.020** in hippocampus and therefore stays off any main figure.

  Two things about step 12 are worth knowing before reusing this machinery. First, the
  statistic is `min(median PSI, 1 - median PSI)`: PSI orientation is undocumented for the
  LIBD tables, and this is invariant to it. Second, the analysis n is the **aging bundle**
  (238 adult controls in hippocampus), not the 452 sample columns in the PSI file — the
  same sample definition as every other aging analysis here.

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
