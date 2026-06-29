# IsoGraph module trust funnel

Modular analysis summary for Manubot integration. Generated from
`real_data/stability/_m/module_trust/` at repo commit `c5df65c`. Every numeric claim
below is reproduced from the result parquet files named in each section; do not edit the
numbers by hand — regenerate from the tables.

## Purpose

Establish *per-module* trustworthiness for IsoGraph isoform-switch modules instead of a
single global partition-similarity score. Global split-half ARI is structurally unfair to
a fine-grained, low-SNR, multi-isoform-only partition and is therefore used only as a
relative A/B dial, not as the scientific readout (see `MODULE_TRUST_PLAN.md` §1). The
funnel asks four sequential questions and only passes survivors downstream:

1. **Q1 — stability:** which production modules recur above a chance-calibrated null?
2. **Q2 — drivers:** are the switch-axis driver loadings of trusted modules reproducible?
3. **Q3 — replication:** do their aging associations replicate cross-cohort (BrainSEQ↔GTEx)?
4. **Q4 — complementarity:** do survivors carry isoform-switch biology distinct from a
   gene-abundance WGCNA baseline?

## Inputs

- **Split-half ensemble** — `real_data/stability/_m/partitions/` (5 seeds × 2 halves ×
  region × method) and the per-fit sidecars `modules_meta/` (eigengene, Age effect, driver
  loadings/kME). Supplies the resampling evidence for Q1 and Q2.
- **Production full-data fits** — `real_data/{brainseq,gtex}/<region>/_m/isograph_vae/`
  (full-multiplex primary) and `wgcna_gene/` modules. These are the reference modules whose
  trust is reported and the basis for Q3/Q4.
- **De-confounded gene-level age test** — the DTU-without-DGE gene sets feeding Q4.
- Regions analysed: BrainSEQ caudate, hippocampus, DLPFC and the paired GTEx regions
  caudate_basal_ganglia, hippocampus, frontal_cortex_ba9.

## Methods text

IsoGraph isoform-switch gene modules were fit on the full sample of each region
(full-multiplex variant, Leiden resolution 5.0) and, independently, on five random
split-half resamples per region. For each production module we quantified recurrence across
the split-half ensemble of the same region and cohort using module-restricted co-assignment
density (the mean fraction of within-module gene pairs that remained co-clustered across
half-fits) together with the best-match Jaccard against the most-overlapping half-fit module.
A module was deemed *trusted* when its co-assignment density exceeded a size-matched
permutation null (1,000 label permutations) at a Benjamini–Hochberg false discovery rate
below 0.05. The identical procedure was applied to a classical gene-abundance WGCNA baseline
on the same samples. For each trusted module we then measured switch-axis driver
reproducibility as the Spearman correlation between per-transcript PC1 switch loadings of the
production module and its best-matching half-fit module over shared genes. Cross-cohort aging
replication matched each trusted BrainSEQ module to the best GTEx module in the paired region
by gene-set and top-5 driver-transcript overlap and tested concordance of the module
eigengene–Age spline effect (sign agreement and joint significance on the independent Salmon
and RSEM quantifications). Finally, survivors were characterised for complementarity to the
abundance baseline by the fraction of member genes carrying a switch-without-differential-
expression signal, their overlap with age-associated WGCNA modules, and the structural class
of their driver isoform switches (CDS, UTR, biotype, and coding-status changes). Analyses
used IsoGraph v0.1.5 under Python 3.12.13 (numpy 2.4.4, pandas 2.3.3, scipy 1.17.1) with a
fixed random seed of 13.

## Results text

**Q1 — stability.** Across the six analysed regions IsoGraph yielded 266 production modules,
of which **236 (89%)** were trusted above the permutation null at FDR < 0.05; the matched
WGCNA baseline yielded 73 modules, of which **64 (88%)** were trusted. IsoGraph therefore
delivers ~3.6× as many chance-calibrated trustworthy modules at a comparable trusted
*fraction*. Per-module best-match Jaccard is far lower for IsoGraph (median ≈ 0.02–0.04) than
WGCNA (median ≈ 0.44–0.72), the expected signature of a finer partition over a low-SNR
switch signal; the chance-calibrated co-assignment test — not raw Jaccard — is the trust
criterion, precisely because Jaccard is granularity-confounded
(`module_stability__*__{isograph,wgcna}.parquet`).

**Q2 — driver reproducibility.** Among trusted modules the switch-axis driver loadings are
highly reproducible across split halves: the median shared-gene driver-loading Spearman ρ is
**0.77–0.82** in every region (caudate 0.78, DLPFC 0.77, hippocampus 0.79; GTEx
caudate_basal_ganglia 0.82, frontal_cortex_ba9 0.81, hippocampus 0.79), with 96–100% of
module pairs showing positive ρ. The *mechanistic identity* of a module (which isoforms
switch) is stable even where hard gene-membership Jaccard is low
(`within_cohort__*__isograph.parquet`).

**Q3 — cross-cohort aging replication.** Matching trusted BrainSEQ modules to GTEx modules in
the paired region and testing the aging effect on independent quantifiers, IsoGraph modules
replicate (sign-concordant **and** jointly significant) at **3/45 caudate, 9/35 DLPFC, and
13/50 hippocampus = 25 replicating modules**; sign concordance alone holds for 69/130. The
matched WGCNA baseline replicates at **1/21, 5/25, 0/7 = 6 modules**. Surviving the
Salmon↔RSEM quantifier gap is genuine biological replication, and IsoGraph carries ~4× as
many cross-cohort-replicating aging modules as the abundance baseline
(`module_aging_replication__*__{isograph,wgcna}.parquet`). The pooled-Stouffer table is empty
by design: with three region pairs the pooled permutation cannot reach p < 0.05 (n = 3 floor),
so replication is reported per region, not pooled.

**Q4 — complementarity.** For age-significant trusted modules the structural annotation of
the top driver switches is near-universal — CDS-change, biotype-switch, and
coding-status-change driver fractions are ≈ 1.0 and UTR-change ≈ 0.6–1.0 — confirming the
drivers are bona fide isoform switches rather than abundance artefacts. The module-level
DTU-without-DGE fraction is modest (median ≈ 0 per region) and overlap with age-associated
WGCNA modules is region-dependent (median frac_in_WGCNA-age 0.39 caudate → 1.00 hippocampus),
consistent with the project-wide honest read that gene abundance dominates the bulk aging
signal and IsoGraph contributes a complementary, mechanistically-resolved switch layer rather
than a globally-superior partition (`module_complementarity__*__isograph.parquet`).

**Headline:** *Of IsoGraph's 236 chance-trusted aging modules across six brain regions, the
switch drivers reproduce across resamples (ρ ≈ 0.77–0.82) and 25 modules replicate their
aging association across independent cohorts and quantifiers — ~4× the matched
gene-abundance baseline — with driver switches that are genuine structural isoform changes.*

## Figure and table notes

- **Main figure — module trust funnel (`real_data/stability/_m/figures/figTrustFunnel.{pdf,png}`,
  built by `real_data/stability/_h/trust_funnel_figure.R`).** Single full-width figure,
  four panels left→right mirroring the funnel, no in-panel titles (interpretation in caption):
  - **(A) Q1 stability:** per-module co-assignment density vs the size-matched null, IsoGraph
    vs WGCNA, with the trusted count annotated (236/266 vs 64/73). Dot/strip over a null band;
    do **not** plot raw Jaccard as the headline (granularity-confounded).
  - **(B) Q2 drivers:** distribution of shared-gene driver-loading ρ per region (violin/box
    hybrid), reference line at ρ = 0, showing the 0.77–0.82 medians.
  - **(C) Q3 replication:** paired BrainSEQ vs GTEx age-effect scatter for matched modules,
    colour by replicates/sign-match; annotate 25 (IsoGraph) vs 6 (WGCNA) replicating.
  - **(D) Q4 complementarity:** stacked structural-switch driver fractions
    (CDS/UTR/biotype/coding-status) for the replicating modules.
  - Rationale / key message: IsoGraph produces many trustworthy aging modules whose
    mechanistic identity reproduces and replicates cross-cohort, complementary to abundance.
- **Supplementary table:** `module_stability__*`, `within_cohort__*`,
  `module_aging_replication__*`, `module_complementarity__*` concatenated per method/region —
  full per-module trust ledger (module_id, n_genes, coassign_density, perm_p, fdr, trusted,
  driver ρ, cross-cohort age effects, replicates, structural-switch fractions).

> The figure is rendered and committed; the script is data-driven (regenerate after any
> re-fit of the module_trust tables). It reuses the project `theme_pub` / `save_fig` /
> Okabe-Ito conventions from `isograph_benchmark/figures/synthetic_benchmark.R`.

## Reproducibility information

- Analysis directory: `real_data/stability/_m/module_trust/`
- Primary script: `isograph_benchmark/real_data/module_trust.py`
  (subcommands `stability`, `within`, `meta`, `replication`, `replication-pooled`,
  `complementarity`)
- SLURM drivers: `real_data/stability/_h/01.stability_isograph.sh`,
  `02.stability_wgcna.sh`, `03.stability_aggregate.sh`, `04.module_meta.sh`; cross-cohort
  replication driver `module_trust_replication.sh` (array idx 1–8).
- Inputs: `partitions/`, `modules_meta/`, production `isograph_vae/` and `wgcna_gene/`
  modules, de-confounded gene-level age test.
- Outputs: `module_stability__*`, `within_cohort__*`, `module_aging_replication{,_pooled}__*`,
  `module_complementarity__*` parquet (per cohort × region × method).
- Key parameters: permutation null n = 1,000 (stability) / 10,000 (pooled); trust FDR < 0.05
  (Benjamini–Hochberg); top-k drivers k = 5; random seed 13.
- Git commit: `c5df65c`. Last replication run logged 2026-06-28
  (`real_data/stability/_m/logs/mtrust-rep-41805157_*.log`).
- Compute environment: PSC Bridges-2 RM-shared; IsoGraph v0.1.5, Python 3.12.13,
  numpy 2.4.4, pandas 2.3.3, scipy 1.17.1
  (`/ocean/projects/bio260021p/shared/opt/envs/isograph`).
- Missing reproducibility information: per-package versions are taken from the live
  environment, not captured in a per-run lockfile; the logs record completion timestamps but
  not a full `pip freeze`.

## Limitations and integration notes

- Trust is chance-calibrated, never hand-set; raw Jaccard is reported only descriptively
  because it penalises IsoGraph's finer granularity.
- Cross-cohort replication is deliberately tested across quantifiers (Salmon↔RSEM); the
  quantifier gap is the test, not a hidden nuisance. Pooled Stouffer is unavailable at n = 3
  region pairs, so the cross-cohort claim rests on per-region counts.
- Q4 complementarity here is the *aging* read; it is partial and region-dependent (abundance
  dominates the bulk signal). The stronger complementarity evidence is the separate
  GO-invisible-gate and incremental-association analyses (SCZD/caudate DTU-without-DGE) —
  integrate this funnel with those summaries and with the matched-baseline QTL-anchoring
  contrast before drafting the final Results section. This summary supplies the per-module
  trust scaffold; those supply the disease-axis complementarity headline.
- This is a primary analysis (the per-module trust scaffold the biological claims rest on),
  not a sensitivity check.
