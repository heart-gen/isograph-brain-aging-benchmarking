# IsoGraph module-trust analysis plan (real data)

**Status:** design doc for review (2026-06-18). No implementation yet.

## 1. Goal and reframing

The deliverable is **not** "make IsoGraph's global split-half ARI match WGCNA." Global
ARI/NMI on hard partitions is only a *relative dial* for A/B-ing method changes against
IsoGraph's own baseline (estimability helped; consensus did not). It is the wrong altitude
for the science and structurally unfair to a method that fragments/relabels modules over a
low-SNR, multi-isoform-only signal (see the caveat in `02_module_discovery/README.md`).

The real questions are per-module:

1. **Q1 — Which IsoGraph modules are stable enough to trust?**
2. **Q2 — Are their driver transcripts reproducible?**
3. **Q3 — Do their aging associations replicate (within- and cross-cohort)?**
4. **Q4 — Do they reveal biology WGCNA cannot (DTU-without-DGE)?**

These form a **funnel**: each stage filters, and only survivors are validated downstream.
The headline stops being "IsoGraph's ARI is low" and becomes "*K trustworthy modules
replicate their aging biology across cohorts, of which X% are invisible to WGCNA.*"

```
all modules
   │  Q1 stability gate (vs permuted-label null)
   ▼
trusted modules ──► Q2 driver-transcript reproducibility
   │  Q3 aging-association replication (within + cross cohort)
   ▼
age-replicating trusted modules
   │  Q4 WGCNA-complementarity (no gene-level age DE) + GO
   ▼
headline biology (the complementary DTU layer)
```

## 2. Data sources

- **Split-half fits** (`03_module_trust/_m/stability/partitions/`, 5 seeds × 2 halves × 6 regions,
  per method incl. `isograph`, `isograph_reliability`, `isograph_tin`): the resample
  ensemble for Q1 (stability) and Q2 (driver reproducibility, within-cohort).
- **Production full-data fits** (`02_module_discovery/{brainseq,gtex}/<region>/_m/isograph_vae/modules.parquet`
  and `wgcna_gene/modules.parquet`): the reference modules whose trust we report, and the
  basis for Q3 cross-cohort replication (`REGION_PAIRS`) and Q4 WGCNA comparison.
- **DGE / GO**: existing de-confounded gene-level age test (the DTU-without-DGE gene sets)
  and the GO-consistency stage for Q4.

The trust funnel is anchored on the **production** modules; the split-half ensemble supplies
the stability evidence used to gate them. (A production module is "trusted" if its gene set /
eigengene / drivers recur across the split-half ensemble of the *same* region+cohort.)

## 3. Required infrastructure change (step 0)

The split-half harness currently persists only `[gene_id, module_id]` (`_write_partition`).
Q2–Q3 need each fit to also emit, per module:

- **module eigengene** (per-sample) — already computed by `compute_trait_associations`
  (base.py:94) and returned as `eigengene_table`; just discarded by the harness.
- **module–Age association** (effect, p) — already in the trait-association table.
- **top driver transcripts** — the per-transcript PC1 **loadings** are computed in
  `gene_switch_coordinates` (switch.py:35, `vh[0]`) but dropped in `channels.py`. Surface
  them (e.g. return a `gene_id, transcript_id, loading` table) so each switch gene's
  dominant isoform(s) are recoverable; a module's driver transcripts = top-|loading|
  transcripts of its hub genes.
- **hub / kME** — per-gene intramodule membership (correlation of the gene's switch score to
  the module eigengene) to rank drivers; analogous to WGCNA kME.

This is **cheap and re-fit-free** — all of it is a by-product of the existing fit; we persist
richer partition artifacts instead of the 2-column table. Proposed: write a per-fit
`partitions/<method>__<cohort>__<region>__seed<k>__<half>.parquet` (genes) **plus** a sidecar
`modules_meta/<...>.parquet` (per-module eigengene corr inputs, Age effect, driver transcripts).

## 4. Q1 — module stability gate

For each production module *m* (gene set G_m), quantify recurrence across the split-half
ensemble of the same region+cohort:

- **Best-match Jaccard**: for each half-fit, J(G_m, G') for the best-overlapping module G'
  in that fit; summarize as median over the ensemble. (Per-module analog of ARI — a few
  brittle modules no longer sink the good ones.)
- **Eigengene reproducibility**: correlate the production module eigengene with the
  best-matching half-fit eigengene (sign-aware); summarize median |r|.
- **Co-assignment density**: mean over gene pairs in G_m of the fraction of half-fits in
  which the pair stays co-clustered (module-restricted consensus).
- **(optional) Preservation Z** (Langfelder Zsummary density+connectivity) for a
  WGCNA-comparable statistic.

**Trust definition (chance-calibrated, not arbitrary):** build a null by permuting module
labels (or drawing random gene sets of size |G_m|) and recomputing the same statistics; a
module is **trusted** if its recurrence exceeds the null at FDR < 0.05 (e.g. permutation
p on best-match Jaccard / co-assignment density). Report both the calibrated set and a
descriptive cut (e.g. median Jaccard ≥ 0.5) for transparency.

**Output:** `module_stability.parquet` (method, cohort, region, module_id, n_genes,
median_jaccard, eigengene_r, coassign_density, preservation_Z, perm_p, trusted) and the
headline *"K of M modules trusted above chance"* per region.

## 5. Q2 — driver-transcript reproducibility

For each **trusted** module, take its top-k driver transcripts (top |loading| among hub
genes) and measure recurrence across the split-half ensemble:

- **Driver Jaccard / rank concordance** of the top-k transcript set vs the best-matching
  half-fit module's drivers.
- Report at the **transcript** level (IsoGraph-specific; WGCNA has no transcript drivers) —
  this tests that the module's *mechanistic identity* (which isoforms switch) is stable, not
  just its gene membership.

**Output:** `module_drivers.parquet` (module_id, driver transcripts, loadings, driver
reproducibility) + summary: fraction of trusted modules with reproducible drivers.

## 6. Q3 — aging-association replication

- **Within-cohort:** sign/magnitude concordance of the module eigengene–Age effect across
  split-halves (does a trusted module age in the same direction on independent halves?).
- **Cross-cohort (primary):** match each trusted BrainSEQ module to a GTEx module in the same
  region (`REGION_PAIRS`) by driver-transcript / gene-set overlap, then test concordance of
  the Age effect (sign + correlation of effect sizes). Surviving the Salmon↔RSEM quantifier
  gap = genuine biological replication and cleanly separates "confound scrambled the
  partition" from "the aging biology is real."

**Output:** `module_aging_replication.parquet` (module_id, region, within_concordant,
cross_match_module, age_effect_bs, age_effect_gtex, concordant, sign_match) + the count of
trusted modules whose aging signal replicates cross-cohort.

## 7. Q4 — biology WGCNA cannot see

For trusted, age-replicating modules:

- **DTU-without-DGE:** show their genes are (a) **not** age-differentially-expressed
  (de-confounded gene-level test) and (b) **not** captured by any age-associated WGCNA
  module (gene-set overlap with WGCNA modules). The IsoGraph-only age signal is the
  complementary layer (cf. the existing 34 SCZD / 43 caudate DTU-without-DGE genes).
- **Function:** GO enrichment (reuse the GO-consistency stage) to show the survivors carry
  coherent, interpretable biology.
- **Structural annotation:** For switch modules review isograph explain-module and
  annotate-structure to understand splicing-level biology.

**Output:** `module_complementarity.parquet` (module_id, frac_genes_no_age_DGE,
overlap_with_wgcna_age_modules, top_GO_terms) + the headline figure.

## 8. Fairness guardrails

- Trust thresholds are **chance-calibrated** (permutation null), never hand-set to flatter.
- Cross-cohort replication is reported on **independent** quantifiers (the confound is the
  test, not hidden).
- Method levers (estimability, TIN) are judged by **how many modules clear Q1 + replicate at
  Q3**, not by global ARI — this is the correct objective function going forward.
- Stability scores are exposed per module (we never present a consensus partition as "the
  answer"; that would inflate reproducibility by construction).

## 9. Open decisions (to settle before building)

1. Trust statistic: Use co-assignment density + permutation p as primary (robust to module
   count/granularity).
2. k for "top driver transcripts" (set to 5) and hub definition (kME threshold).
3. Cross-cohort module matching rule (driver overlap).
4. Scope: caudate-first single-region pilot (matches the TIN/estimability pilot) then fan to
   all 6 regions.

## 10. Sequencing and effort

1. **Step 0** — enrich persisted artifacts (eigengene, Age effect, driver loadings/kME).
   Small IsoGraph + harness change; re-run the cheap split-half fits to regenerate sidecars.
2. **Q1 + Q2** — stability gate + driver reproducibility from the ensemble. New analysis
   module (`module_trust.py`) + outputs. No heavy compute.
3. **Q3** — cross-cohort aging replication from production fits. No new fits.
4. **Q4** — WGCNA-complementarity + GO. Reuses DGE/GO infra.
