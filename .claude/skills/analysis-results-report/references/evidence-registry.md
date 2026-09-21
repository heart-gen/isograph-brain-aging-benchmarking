# Evidence registry

Where each biological question's evidence lives. Paths are repository-relative; address
stage trees through `isograph_benchmark/paths.py` (`stage_out`, `region_store`,
`region_artifact_dir`) rather than literally.

Every stage report has the same shape — **§5** narrative, **§6** the numbers table with
its regeneration date, **§7** the counter-evidence, **§12** the audit trail with file
paths. When a number is needed, go stage report §6 → §12 → the result file.

---

## Question 1 — Recovery where ground truth is known *(synthetic)*

| | |
| --- | --- |
| Report | `reports/pi/01_synthetic_benchmark.md` |
| CLIs | `benchmark/run_synthetic.py`, `run_one.py`, `interpret_modules.py`, `residual_cost.py`; `stats/summarize.py`, `sweep_breakdown.py` |
| Results | `01_synthetic_benchmark/01_synthetic/_m/synthetic_results.parquet`; `03_metrics/_m/synthetic_metric_{long,summary}.parquet`, `synthetic_pairwise_tests.parquet`, `tableS_benchmark_summary.csv` |
| Figures | `01_synthetic_benchmark/03_metrics/figures/fig1_benchmark_overview.{pdf,png}`, `figS1`–`figS12` |
| Display | Fig 1D–E, S1–S12, `tableS_benchmark_summary` |
| Bound | Wins where switching dominates (7/15 scenarios), not everywhere. The archived datasets are **not regenerable** (generator-pinned after 2026-08-14); the archive is the record. |

Also here: the A1 synthetic-genetics arm — WGCNA fails module recovery on identical
features (mean 0.47 against ~1.0 for IsoGraph; one-sided paired Wilcoxon P = 3.0e-18,
n = 117 matched datasets — manuscript repo `drafts/analysis-summaries/A1-synthetic-genetics.md`,
recomputed from `synthetic_results.parquet`, `run_scenario == "genetic_anchoring"`), the cleanest features-vs-inference separator.

---

## Question 2 — A measurable layer, separable from abundance

| | |
| --- | --- |
| Reports | `reports/pi/02_module_discovery.md`, `03_module_characterization.md` |
| CLIs | `real_data/run_models.py`, `module_sizes.py`, `interpret_modules.py`, `module_enrichment.py`, `incremental_association.py`, `abundance_structure_separation.py`, `go_invisible_gate.py` |
| Results | `<store>/isograph_vae/`, `<store>/{module_interpret,module_enrichment,incremental_association}/`; `03_module_characterization/_m/incremental_effect_sizes.parquet`; `<caudate_sczd store>/go_invisible_gate.parquet` |
| Figures | `manuscript/_m/figures/figConceptOverview`, `figSeparation`, `figGoInvisible` |
| Display | Fig 1A–C, S-real-1, S-real-3, S-real-4, S-real-5, S6 |
| Baselines | `wgcna_gene` (classical abundance), `wgcna_switch_only`, `wgcna_multiplex` — the last two receive **identical** switch features, so they separate features from inference |
| Guard | Leiden module ids are re-assigned every fit; `real_data/partition_provenance.py` fingerprints each partition and the join is only valid against the fit that wrote it |

Generation: switching transcript filter since 2026-09-14, Leiden resolution **2.0** since
2026-09-16. `leiden_sweep` is `partial_coverage` (1 of 17 stores) — say so if cited.

---

## Question 3 — Trust and reproducibility

| | |
| --- | --- |
| Report | `reports/pi/04_module_trust.md` (§6 has the full numbers table) |
| CLIs | `real_data/{stability,module_trust,replication,replication_go,replication_permutation,replication_functional}.py` |
| Results | `04_module_trust/_m/stability/{partitions,stability_summary.parquet,module_trust,modules_meta,lr_validation}`, `_m/replication/`; `_m/baseline_comparison/` |
| Figures | `manuscript/_m/figures/figTrustFunnel` (Fig 2), `figCrossCohortReplication` (supplement), `figBaselineRates` |
| Randomness | 5 split-half seeds (default 13); 10,000 permutation draws per null cell |

Report together, always:

- IsoGraph chance-trusted **93/118 (79%)**; WGCNA **62/72 (86%)** — WGCNA is higher.
- Split-half age concordance **22/23** both-significant pairs (WGCNA 14/14).
- Cross-cohort **eigengene projection 33/38 and 44/55** (P = 2.1e-6, 4.4e-6) vs matched
  WGCNA null (28/52 P = 0.34; 3/10 P = 0.95) — the strongest inference-specific result.
- Cross-cohort module-**pair** count **null and retired** (1/38, perm P = 0.91).
- Curvature 0/118 — aging is linear at module scale, which licenses the linear model.

The held-out DTU added-value test (module context vs co-expression, partial r 0.056–0.102
in 3/6 analyses) also lives here; frame it per `claim-calibration.md` §4.

---

## Question 4 — Are the switches real molecular events?

| | |
| --- | --- |
| Reports | `reports/pi/06_switch_mechanism.md`, `06a_allelic_switch_arm.md` |
| CLIs | `real_data/{switch_consequence,switch_consequence_meta,validate_switch_splicing,switch_orthogonal_confirm,longread_switch_confirm,isa_concordance,clinical_consequence}.py`; `ase_junction_switch.py`, `ase_junction_allelic.py`, `ase_risk_orientation.py` |
| Results | `06_switch_mechanism/_m/{switch_consequence_meta,clinical_consequence_meta}.parquet`; `_m/{switch_validation,switch_orthogonal_confirm,longread_switch_confirm,isa_concordance}/`; `_m/ase_junction_switch/<region>/{junction_allelic_counts,allelic_test,risk_orientation}.parquet` |
| Figures | `manuscript/_m/figures/figOrthogonalConfirm`, `figSwitchConsequence`, `figIsaConcordance`, `figClinicalConsequence` |
| External | ONT DLPFC BA9/46 Bambu quants (Aguzzoli-Heberle 2024 NBT, Zenodo 8180677, n = 12); GENCODE v47 |
| Nulls | 2,000 matched-null draws |

Headline: long-read switch-like **0.675 vs a 0.243 matched null** over 206 pairs,
**26/30** genes individually. Allelic arm: isoform choice under within-donor cis control
in **185/92/37** genes (caudate/hippocampus/DLPFC) against a flat homozygous-lead null —
the project's only composition-proof genetic evidence. It is a **gene-level** statement: the
cis signal is not concentrated in the aging modules (pooled ρ = −0.183, permutation p = 0.47;
06a §5c). Disease direction is a **set-level negative** measured two ways — through LD from the
QTL lead (481/556 unorientable, median |r| 0.03–0.08) and at the GWAS lead itself (371/556
testable, 214/371 positive, chance). Never orient at the QTL lead. Two single-transcript
effects are reported against that null: PICALM-219 × AD (q = 0.015, a Result) and KLC1-213 ×
SCZ (a caveat). The allelic recount ran on Quest (`b1042`), not Bridges.

---

## Question 5 — Cell-type composition

| | |
| --- | --- |
| Report | `reports/pi/03_module_characterization.md` |
| CLIs | `real_data/celltype_composition.py`, `characterize_composition_unique.py`; MuSiC deconvolution in R (`gtex_music_deconv.R`) |
| Results | `03_module_characterization/_m/composition_adjustment.parquet`, `COMPOSITION_ADJUSTMENT_SUMMARY.md`; `02_module_discovery/gtex/_m/composition/` |
| Figures | `manuscript/_m/figures/figCompositionRobustness` (**Fig 5**), Table S13 |

Report the survivals **and** the collapses: BrainSEQ caudate 45→14, DLPFC 22→28;
6 of 8 deconvolved GTEx regions survive (hippocampus 44→29, caudate 87→50, ACC 928→365);
**the two GTEx cortical regions collapse** (frontal cortex BA9 797→2, cortex 1,169→2); the
SCZD disease signal is largely composition-confounded (60→11, 18% retained). Aging
hippocampus contributes 1 gene. Confounder vs mediator is **unresolvable from these data**
and the report must say so.

---

## Question 6 — Genetics of the switch layer

| | |
| --- | --- |
| Reports | `reports/pi/05_genetic_anchoring.md`, `05a_signal_level_genetics.md`, `08_integration.md` |
| CLIs | `real_data/{qtl_anchoring,qtl_anchoring_meta,sqtl_concordance,module_genetic_anchoring,coloc_*,ldsc_*,coloc_modality_contrast,coloc_signal_susie,coloc_brainseq,smr_heidi,brainseq_switch_qtl,locus_event_audit,gene_deep_dive,anchored_gene_summary}.py`; `isograph_benchmark/gwas/` |
| Results | `05_genetic_anchoring/_m/{qtl_anchoring_meta,coloc,coloc_modality_contrast,coloc_signal_susie,coloc_brainseq,smr_heidi,brainseq_switch_qtl,ldsc,gwas,module_coloc_convergence}/`; `08_integration/_m/{deep_dive,anchored_gene_summary}/` |
| Figures | `manuscript/_m/figures/figQtlSpecificity` (Fig 3), `figGeneticAnchoring` (Fig 4), `figSczConvergence`, `figGwasResolution` |
| Inputs | GTEx v11 brain sQTL/eQTL + all-pairs (the sQTL filename carries an extra `cis_sqtl.`; sGenes is one row per gene); 5 GWAS incl. PGC3 EUR for SCZ; `g1000_eur` hg19 panel with `NCBI37.3.gene.loc` |
| Tools | SMR 1.4.2 at `shared/opt/SMR/build/Release/smr`; coloc in `R_env` |

Rules specific to this question:

- S-LDSC on **`coef_p`**, never `enrichment_p` (register D-2026-09-10-a). PD splicing-led
  (coef p = 0.0050, survives Bonferroni and the joint model); SCZ expression-led.
- Set-level splicing specificity **1.065** (p = 0.020) with `wgcna_multiplex` **1.084**
  (p = 0.0044) in the same breath — it is a representation effect.
- Per-gene counts run the other way and belong in the same section.
- Coloc SNP guard: 12,000 primary; AD 30,000 is the sensitivity root
  `sensitivity/max_snps_30000/`.
- Fig 4A is **PRDM2**. SNCA is a falsification example. Fig 4E's SCZ convergence panel is
  retracted; the replacement is the null permutation result (P = 0.19–1.00 in all 10 cells).
- Anchored layer: **160** colocalized genes, **30** resolve to a specific switch pair, only
  **6** reach CLPP ≥ 0.05; 76% of concordant events are cerebellar, tracking GTEx sample
  size rather than disease-relevant tissue.

**RBP regulation** (`reports/pi/07_rbp_regulation.md`; `07_rbp_regulation/_m/rbp/`,
`_m/neuronal_clip/`; `figRbpRegulon`) enters only as a secondary annotation layer, with the
HepG2/K562 eCLIP tissue mismatch stated.

---

## Cross-cutting

| Source | Use |
| --- | --- |
| `ANALYSIS_MAP.md` | analysis → CLI → wrapper → output → display item |
| `reports/pi/_evidence/inventory.md` | status board; last content-change commit |
| `reports/pi/_evidence/DISCREPANCY_REGISTER.md` | numbers a document gets wrong about its own output |
| `manuscript/FIGURE_ORDERING.md` | the one honest claim each display item carries — reuse its wording |
| `manuscript/MANUSCRIPT_PLAN.md` | findings, framing, section plan |
| `../../manuscript/isograph-brain-manuscript/{TODO.md,content/}` | what the prose still says and what is queued to change |
| `reports/pi/pi_briefing.tex` | Beamer preamble, palette and `\fig` / `\srcnote` macros to reuse |

**Run dates:** `_m/` filesystem mtimes are a checkout artifact, not run dates. Use the last
git content-change commit, or SLURM log mtimes (gitignored, written by the jobs themselves).
