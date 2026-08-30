# Real-data analyses (BrainSEQ + GTEx)

This directory applies IsoGraph and a gene-level WGCNA baseline to two postmortem
human-brain cohorts and characterizes the resulting modules against age, diagnosis,
biological function, and (GTEx) cross-region conservation. It consumes the dataset
bundles built in `inputs/` (see `inputs/README.md`) and the IsoGraph software described
in `isograph_benchmark/README.md`.

> **Status (draft).** This README documents the analysis as it currently stands. The
> GTEx IsoGraph results are complete; the GTEx WGCNA baseline is being **re-run** after
> two fixes to its R script (a soft-power selection bug and a spline-FDR bug — see
> [WGCNA fixes](#wgcna-fixes)), and the GTEx downstream stages that consume WGCNA
> (module enrichment) regenerate from it. The cross-cohort
> replication IsoGraph arm is complete; the WGCNA arm now has a BrainSEQ aging WGCNA
> baseline (new, running) and regenerates for both methods once it finishes — see
> [Cross-cohort replication](#cross-cohort-replication-brainseq-vs-gtex).
> Numbers that depend on WGCNA should be treated as provisional until that chain finishes.

## Cohorts

| Cohort | Regions | Trait | Notes |
|---|---|---|---|
| **BrainSEQ** | `caudate` (Phase 3), `hippocampus`, `dlpfc` (Phase 2) | Age | Controls only, adults (Age ≥ 18) |
| **BrainSEQ SCZD** | `caudate_sczd` | Dx | Control + schizophrenia |
| **GTEx v11** | 13 brain regions (`amygdala` … `substantia_nigra`) | AGE | Exact age (v8-preferred) |

Covariates (from `configs/real_data.yaml`): BrainSEQ adjusts for Sex, MoD, RIN, mapping
rate, mito rate, and SNP PC1–PC5; GTEx adjusts for SEX, SMRIN, SMTSISCH, SMMAPRT. Age
is modeled both linearly (eigengene–age correlation) and with a natural-spline
association (df = 4), with the spline compared against the linear fit.

## Methods

Both cohorts fit **IsoGraph** (VAE backend, switch channel; optional abundance channel)
and a **gene-level WGCNA** baseline on the same samples and bundle-filtered genes, then
associate each module's eigengene with the cohort trait via the linear and spline models
above (Benjamini–Hochberg FDR across modules). Modules are further characterized by:

- **Interpretation** — switch-driver transcript tables, transcript polarity, and
  high-vs-low contrasts per selected module (`interpret_modules`), with GENCODE v47
  structural annotation of switch pairs.
- **Incremental association** — a de-confounded test of whether IsoGraph modules add
  age/diagnosis signal *beyond* gene-level abundance, reported at both gene and module
  level (`incremental_association`).
- **Module enrichment** — per-module GO:BP enrichment plus IsoGraph network metrics,
  joined to each module's age-spline phenotype FDR, for both methods
  (`module_enrichment`).
- **Leiden resolution sweep** — cheap re-clustering of saved edges to choose module
  granularity without refitting the VAE (`sweep_leiden`).

The Python implementations live in `isograph_benchmark/real_data/` (`run_models.py`,
`incremental_association.py`, `module_enrichment.py`, `sweep_leiden.py`,
`interpret_modules.py`, `go_enrichment.py`, `characterize_composition_unique.py`); WGCNA
baselines are R scripts under each cohort's `_h/`.

## Pipelines

### GTEx (`02_module_discovery/gtex/_h/`)

| Step | Script | Output |
|---|---|---|
| 01 | `01.run_isograph.sh` | IsoGraph fit per region: `modules`, `edges`, `traits`, `age_{linear,spline}`, `feature_scores`, `module_gene_roles`, `calibration` |
| 02 | `02.wgcna_gene.{R,sh}` | WGCNA baseline: `modules`, `age_{linear,spline}` |
| 04 | `04.interpret_modules.sh` | Per-module interpretation tables + structure annotations |
| 05 | `05.incremental_association.sh` | Gene- vs module-level incremental association |
| 06 | `06.module_enrichment.sh` | GO:BP + network metrics + phenotype FDR, both methods |

All 13 regions have run through stages 01–06. Current IsoGraph module counts: 14–29 per
region (giant module 21–39%).

### BrainSEQ (`02_module_discovery/brainseq/_h/`)

| Step | Script | Purpose |
|---|---|---|
| 01 | `01.run_isograph_aging.sh`, `01.aging_models.R` | IsoGraph fit + aging models on the 3 control aging regions |
| 02 | `02.run_isograph_sczd.sh`, `02.run_wgcna_sczd.{R,sh}` | Caudate SCZD IsoGraph (Dx trait) + WGCNA baseline |
| 04 | `04.interpret_modules.sh` | Module interpretation |
| 05 | `05.sweep_leiden.sh` | Leiden resolution sweep |
| 06–07 | `06.run_isograph_with_abundance.sh`, `07.sweep_leiden_with_abundance.sh` | Abundance-channel variant + its Leiden sweep |
| 08 | `08.incremental_association.sh` | De-confounded incremental association |
| 09 | `09.characterize_composition_unique.sh` | Composition-unique module characterization (`_m/composition_unique_overlap.parquet`) |
| 10 | `10.module_enrichment.sh` | GO:BP + phenotype enrichment, both methods |

The BrainSEQ cross-region aggregations are `_m/composition_unique_overlap.parquet` and
`_m/module_interpret_summary.parquet` (overlap/interpretation summaries).

### GWAS overlap (`05_genetic_anchoring/_h/`)

Module gene sets from either cohort are tested for heritability enrichment with MAGMA:
`01.prep_module_gene_sets` → `02.run_magma` → `03.plot_magma` (config `configs/gwas_magma.yaml`).

## Cross-cohort replication (BrainSEQ vs GTEx)

BrainSEQ and GTEx are two independent cohorts (different donors, libraries, and
quantifiers) profiled with the same pipeline. For the three brain regions present in
both, `isograph_benchmark/real_data/replication.py` asks how reproducible the discovered
structure is, per method:

| Match | BrainSEQ | GTEx | Shared expressed genes |
|---|---|---|---|
| Caudate | `caudate` | `caudate_basal_ganglia` | 8,267 |
| Hippocampus | `hippocampus` | `hippocampus` | 6,408 |
| DLPFC / BA9 | `dlpfc` | `frontal_cortex_ba9` | 14,826 |

Gene IDs are versioned Ensembl IDs that match directly across cohorts (same GENCODE
annotation), so analyses run on genes expressed in both cohorts for a region. Two
complementary readouts:

- **Module preservation** — each module's best-matching module in the other cohort
  (highest Jaccard over shared genes), both directions, with significance from a
  size-preserving label-permutation null (1,000 permutations → empirical p, z-score; a
  module is *preserved* at p < 0.05).
- **Age-effect concordance** — for best-matched module pairs, whether the eigengene–age
  correlation (`age_linear` effect) agrees in sign (binomial test) and magnitude
  (Spearman) across cohorts.

Run with `python -m isograph_benchmark.real_data.replication` (driver:
`03_module_trust/_h/10.replication.sh`); outputs under `03_module_trust/_m/replication/`
(`<method>_module_match.parquet`, `replication_summary.parquet`, `replication_summary.json`).

**Results (both arms complete).** Per-module best-match Jaccard and preservation rate:

| Method | mean Jaccard | frac. preserved (vs null) |
|---|---|---|
| IsoGraph | 0.04–0.06 | 0.11–0.47 |
| WGCNA | 0.14–0.36 | 0.67–1.00 |

WGCNA's apparently far higher preservation is **largely a granularity artifact, not
stronger biology**: WGCNA produces a few coarse modules (62–90% of genes in one module),
so a giant module in one cohort trivially overlaps a giant module in the other and easily
beats the size-matched permutation null. IsoGraph's finer partitions (14–29 modules,
21–39% giant) are inherently harder to match by exact best-Jaccard. The raw preservation
metric therefore favors coarse partitions and should not be read as WGCNA finding more
reproducible structure; a granularity-matched comparison (or module-level preservation
statistics) would be needed for a fair head-to-head.

**Age-effect concordance is weak for both methods.** Eigengene–age sign concordance sits
near chance (≈0.2–0.71 across region/direction) and magnitude Spearman is inconsistent
(e.g. IsoGraph DLPFC/BA9 BrainSEQ→GTEx ρ ≈ 0.56; WGCNA values range from −0.89 on a small
matched set to +0.68). The honest read: module-level aging signal replicates only weakly
across these two independent cohorts for either method.

The BrainSEQ aging WGCNA baseline used here did not previously exist (the BrainSEQ aging
regions had only `isograph_vae`); it was added as
`02_module_discovery/brainseq/_h/01.wgcna_gene_aging.{R,sh}`, adapted from the fixed GTEx
`02.wgcna_gene` to BrainSEQ counts (log2-CPM), categorical covariates (Sex, MoD;
complete-case spline handling for NA SNP PCs), and the Age trait (caudate 21 / hippocampus
9 / dlpfc 25 modules).

### Functional (GO) consistency of replicated aging modules

The strongest form of replication is functional: for modules that are *both*
cross-cohort preserved *and* age-associated (eigengene–age linear FDR < 0.10), do the
two cohorts' matched modules enrich for the same biology?
`isograph_benchmark/real_data/replication_go.py` answers this using the **full** enriched
GO:BP term set per module — `module_enrichment.py` now persists
`<method>_module_go.parquet` (one row per module × enriched term) in addition to the
top-term summary, so the comparison uses complete term sets rather than the top-8 names.

For each preserved+aging matched pair it computes the GO-term Jaccard
(|shared terms| / |union|) and the shared term names, and runs a pooled permutation test
(1,000 permutations): is the mean GO Jaccard of the matched pairs higher than when each
source module is paired with a *random* target module from the same cohort/region? Driver:
`03_module_trust/_h/11.replication_go.sh`; outputs `<method>_go_overlap.parquet`
(per pair) and `replication_go_summary.{parquet,json}`.

This requires per-module GO for both cohorts' matched regions; the BrainSEQ aging GO
enrichment (which had never been run — only `caudate_sczd` had module enrichment) was
generated for both methods alongside a GTEx enrichment re-run that adds the full-GO-set
tables.

**Results (both arms complete).** Pooled over the preserved+aging matched pairs:

| Method | preserved+aging pairs | pairs with GO | mean GO Jaccard | median | perm. *p* (vs random pairing; null mean) |
|---|---|---|---|---|---|
| IsoGraph | 6 | 5 | 0.066 | 0.000 | 0.005 (null 0.004) |
| WGCNA | 26 | 26 | 0.201 | 0.077 | 0.001 (null 0.040) |

**Both methods' replicated aging modules are functionally consistent above chance** — the
mean GO Jaccard of the cross-cohort matched pairs significantly exceeds random
source-to-target pairing for both (IsoGraph *p* = 0.005, WGCNA *p* = 0.001). As with the
structural preservation result, WGCNA's higher *absolute* Jaccard is inflated by the same
granularity artifact: its few coarse modules carry very large GO sets (up to ~1,160 enriched
BP terms), which mechanically overlap more than IsoGraph's compact modules. The honest read
is that functional replication is the one place IsoGraph's signal holds up — its finest
preserved+aging modules still reach significant above-null GO consistency, with a concrete
hit in hippocampus (GTEx M001 ↔ BrainSEQ M000: 23 shared BP terms, GO Jaccard 0.21) and in
DLPFC/BA9 (7 shared terms, Jaccard 0.12) — even though most of its 6 pairs share no terms
(median 0). WGCNA spreads functional overlap more broadly but largely because its modules
are larger, not because they localize aging biology more precisely.

### Within-cohort split-half stability (is the cross-cohort failure a confound or method instability?)

The cross-cohort comparison conflates two things: BrainSEQ (Salmon) and GTEx (RSEM) use
different transcript quantifiers, and within-gene isoform ratios — the switch signal
IsoGraph models — are quantifier-sensitive, whereas gene abundance (WGCNA's input) is
quantifier-robust. So a weak IsoGraph cross-cohort replication could be either (a) the
Salmon-vs-RSEM confound scrambling the switch signal while the method itself is stable, or
(b) the switch network simply being hard to estimate reproducibly. To separate them we hold
the quantifier and preprocessing **fixed** (stay inside one cohort), randomly split the
samples 50/50 over 5 seeds, refit the *same* method on each half with the *same* feature
set, and measure partition agreement (ARI, NMI) over genes assigned in both halves — the
same metrics as the cross-cohort test. `isograph_benchmark/real_data/stability.py` fits
IsoGraph (drivers `03_module_trust/_h/01.stability_isograph.sh`, WGCNA reference
`stability_wgcna.R` / `02`, aggregator `03`); outputs `03_module_trust/_m/stability/stability_summary.{parquet,json}`
and `stability_pairs.parquet`.

**Results (mean over 5 seeds; cross-cohort full-data fits shown for contrast).**

| Method | Region | Within-cohort ARI | Within-cohort NMI | Cross-cohort ARI | Cross-cohort NMI |
|---|---|---|---|---|---|
| IsoGraph | Caudate (BrainSEQ) | 0.224 | 0.183 | 0.025 | 0.034 |
| IsoGraph | Caudate (GTEx) | 0.278 | 0.169 | — | — |
| IsoGraph | Hippocampus (BrainSEQ) | 0.087 | 0.113 | 0.009 | 0.027 |
| IsoGraph | Hippocampus (GTEx) | 0.237 | 0.152 | — | — |
| IsoGraph | DLPFC (BrainSEQ) | 0.006 | 0.011 | 0.007 | 0.011 |
| IsoGraph | BA9 (GTEx) | 0.180 | 0.162 | — | — |
| WGCNA | Caudate (BrainSEQ) | 0.403 | 0.532 | 0.057 | 0.235 |
| WGCNA | Caudate (GTEx) | 0.711 | 0.555 | — | — |
| WGCNA | Hippocampus (BrainSEQ) | 0.291 | 0.343 | 0.306 | 0.231 |
| WGCNA | Hippocampus (GTEx) | 0.426 | 0.440 | — | — |
| WGCNA | DLPFC (BrainSEQ) | 0.307 | 0.473 | 0.056 | 0.245 |
| WGCNA | BA9 (GTEx) | 0.685 | 0.536 | — | — |

The test cleanly separates the two explanations, and **both are present**:

- **The quantifier confound is real.** For caudate and the GTEx regions, IsoGraph's
  within-cohort ARI (0.18–0.28) is **5–10× its cross-cohort ARI** (≈0.01–0.025). Holding the
  quantifier fixed substantially restores reproducibility, so part of the cross-cohort
  collapse is genuinely the Salmon-vs-RSEM mismatch corrupting within-gene isoform ratios,
  not the method failing outright.
- **But the switch signal is intrinsically lower-reproducibility than abundance.**
  IsoGraph's *within-cohort* ceiling (NMI ≈ 0.11–0.18) is still **below WGCNA's
  *cross-cohort* NMI (≈0.23)** — i.e. WGCNA replicates better across two *independent*
  cohorts than IsoGraph replicates across split-halves of the *same* cohort. WGCNA's
  within-cohort stability (ARI 0.29–0.71, NMI 0.34–0.55) is 2–3× IsoGraph's throughout.
- **Strong region heterogeneity.** BrainSEQ DLPFC (ARI 0.006) and BrainSEQ hippocampus
  (ARI 0.087, only ~4,900 genes assigned in both halves) are near-zero even within-cohort —
  for the largest/sparsest regions the switch network is effectively unstable at this sample
  size regardless of quantifier. The well-behaved regions (caudate, the GTEx trio) carry
  whatever reproducible switch structure exists.

A structural caveat compounds this: IsoGraph only assigns multi-isoform (switch-capable)
genes, and *which* genes qualify shifts with the split (BrainSEQ hippocampus shares only
~4,900 of ~17k genes between halves), so its partition agreement is computed over a smaller,
less stable gene set, whereas WGCNA partitions all ~18–19k genes. Net read: IsoGraph's weak
cross-cohort replication is a layered limitation — a real quantifier confound on top of an
intrinsically noisier, sparser switch signal — and the honest framing is the one the GO
analysis already supports: IsoGraph is a complementary DTU-without-DGE layer, not a method
that recovers more reproducible modules than WGCNA.

**Caveat — this measures reproducibility, not accuracy, and is not a head-to-head verdict
against WGCNA.** There is no ground truth on real data, so ARI/NMI here compare the *two
split-halves to each other* (self-consistency), never to a true labeling. That is a
different quantity from the synthetic benchmark, where ARI is module *recovery* against the
planted modules and IsoGraph scores well because the switch signal is injected at clean,
high effect size. The synthetic-vs-real gap therefore reflects real-data SNR, not a defect:
the real switch signal is subtle (n≈110/half), confounded, and gated to multi-isoform genes.
Three reasons WGCNA's higher numbers do **not** mean it is "better" at IsoGraph's task:

- **Different signal SNR.** WGCNA clusters gene **abundance** (log-CPM) — smooth, high-SNR,
  partitions near-identically on any half regardless of biological meaning. IsoGraph clusters
  within-gene **isoform-switch coordinates** (PC1 of CLR composition) — low-SNR by
  construction. Lower split-half agreement is the *expected* cost of measuring a harder
  signal, not evidence the method is broken.
- **Different gene sets.** The two ARI/NMI values are computed over different `n_common`
  (~8–9k switch-capable genes vs ~18–19k) and different partition granularities, so they are
  not strictly comparable point-for-point.
- **Different biology.** WGCNA's stability reflects abundance co-expression and says nothing
  about isoform usage; IsoGraph is the only layer measuring switch (DTU-without-DGE)
  structure. "WGCNA is more reproducible" ≠ "WGCNA captures what IsoGraph captures."

The proper use of this test is as a *relative* instrument: A/B-ing changes to IsoGraph
against its **own** baseline (e.g. covariate-free isoform-estimability edge downweighting
improves within-cohort ARI in 5/6 regions and NMI in 6/6; consensus Leiden did not), not as
a cross-method ranking.

## WGCNA fixes

Two bugs were fixed in `02_module_discovery/_h/09.wgcna_gene_gtex.R`:

1. **Soft-power selection.** The baseline initially collapsed into a single giant module
   in every region (76–100% of genes in one module; one module in four regions).
   `pickSoftThreshold` was called without `networkType = "signed"`, so it picked the
   power for an *unsigned* topology (power = 1) while `blockwiseModules` built a
   **signed** network — at power 1 the network is nearly complete and TOM merges
   everything into one module. The fix selects the power for the signed topology and
   floors it at WGCNA's recommended minimum (12) for signed networks at these sample
   sizes. After the fix WGCNA finds 5–16 modules per region (giant module 62–90%; still
   larger than IsoGraph's 21–39%, a property of these data under signed correlation
   rather than a collapse).

2. **Spline-FDR length bug (latent, exposed by fix 1).** `spline_age_assoc` built the
   per-age-percentile BH-FDR column by indexing a per-label list with the full label
   vector, which produced a wrong-length vector and crashed (`replacement has N rows,
   data has M`) for any region with more than one module. It never triggered before
   because WGCNA always collapsed to one module. The fix uses `ave(pvalue, age_label,
   FUN = fdr_bh)`. Without it, `age_spline.parquet` is stale/missing and module
   enrichment (which joins the age-spline FDR) is invalid.

The WGCNA baseline plus the two downstream stages that consume it (region-shared
analysis, module enrichment) are re-running with both fixes; IsoGraph results and the
IsoGraph-only stages (interpretation, incremental association) are unaffected.

The BrainSEQ SCZD WGCNA script (`02.run_wgcna_sczd.R`) shares the same latent soft-power
pattern but produced a healthy baseline (23 modules, 24% giant) and has not been re-run.

## Outputs

Per region: `02_module_discovery/<cohort>/<region>/_m/isograph_vae/` (and
`isograph_vae_with_abundance/` for BrainSEQ), `…/wgcna_gene/`, `…/module_enrichment/`,
and the IsoGraph `module_interpret/` and `incremental_association/` subdirectories.
Cohort-level summaries: `02_module_discovery/gtex/_m/` (region Jaccard, age summaries, interpret
summary) and `02_module_discovery/brainseq/_m/` (composition-unique overlap, interpret summary).
Model outputs under `_m/` are regenerable from the bundles in `inputs/bundles/`; bulky
intermediates (`_o/`, `_m/tmp/`, logs) are git-ignored.

## Reproducibility

Run from the repo root with `ISOGRAPH_BENCHMARK_ROOT` set; IsoGraph and enrichment steps
use the `isograph` env (`/ocean/projects/bio260021p/shared/opt/envs/isograph`), WGCNA and
R analyses use the R env (`/ocean/projects/bio250020p/shared/opt/env/R_env`). The SLURM
scripts under each `_h/` are submitted in numeric order; downstream stages depend on the
model fits, so re-running a baseline (e.g. WGCNA) requires re-running the stages that
consume it. See `configs/real_data.yaml` for cohort filters, covariates, and aging-model
settings.
