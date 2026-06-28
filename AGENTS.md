# AGENTS.md — Real-Data Rigor & Biology for Nature Methods

Working plan for finishing the **real-data** analysis so it is rigorous and
biologically grounded for a Nature Methods submission. This file is the to-do
spine; it is not a status report. Update it as items close.

## North-star framing (do not drift from this)

IsoGraph is a **complementary isoform-switch network method**, *not* globally
superior to WGCNA on every real-data metric. The de-confounded gene-level test
showed abundance dominates (~75–110× more unique signal); IsoGraph's defensible
value is a small, specific **DTU-without-DGE** layer (34 SCZD / 43 caudate genes)
that is structurally invisible to any DGE/WGCNA pipeline. Every real-data claim
must be honest about this. Interpret against **three WGCNA baselines**:
`wgcna_gene` (classical abundance), `wgcna_switch_only`, `wgcna_multiplex`.

## Hard constraints (every agent, every time)

- **Heavy compute on SLURM only.** Login node = light work (reads, small
  aggregation, plotting). Account `bio260021p`, memory fixed at 2000M/cpu.
- **Never `git add -A`.** Stage selectively. Commit only when explicitly asked;
  commit messages end with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
- **Do not push** without explicit confirmation (currently on hold across all
  repos). PR bodies end with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
- **Do not commit large data caches:** `modules_meta/`, `module_trust/`,
  `partitions/`, `inputs/tin/`, `inputs/go_annotations/`, `_m/` outputs, region
  output dirs.
- Config is the single source of truth. Canonical Leiden resolution = **5.0**,
  giant-cap **OFF**. Keep code clean — rationale lives in the wiki, not comments.
- **Every official test / analysis / interpretation must be reproducible.** No
  ad-hoc inline scripts for results that go in the paper: write a committed,
  parametrized CLI module under `isograph_benchmark/real_data/` + a SLURM wrapper,
  with deterministic seeds and outputs written to disk. Inline heredoc Python is
  for exploration only — promote anything you keep.
- Interpreters: python `/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python`;
  figure R `/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript`;
  classical-WGCNA + GWAS R `/ocean/projects/bio250020p/shared/opt/env/R_env`.
- Ocean FS: use `iterdir()` not `glob()`.

## Done (foundation — read, don't redo)

- Canonical res 5.0 re-fit + GWAS (eliminates giant-module artifact: 0/10 sig
  hits ≥900g vs 17/28 at 2.0; SCZ signal survives, concentrates).
- Trust funnel complete at 5.0: 120/120 partitions, `stability_summary` written.
- WGCNA rigor fixes landed + re-run: seed-13 default; corrected spline df
  (`n_obs − rank`) + full-vs-reduced F-test + BH `fdr_ftest`; matched-feature
  baselines `wgcna_switch_only/` + `wgcna_multiplex/` across brainseq 3 + gtex 13.
- Downstream re-run at 5.0: module enrichment (`--method all`), `replication_go`
  (4 methods), `module_meta` loadings (6 regions), `lr_validation` (7 regions).
- De-confounded gene-level test, incremental association, composition-unique
  characterization infrastructure committed.

---

## What is LEFT (ordered by biological priority)

### 1. BIOLOGY GATE — GO-invisible disease switch modules — DONE (PASS, complementary)

Reproducible gate: `python -m isograph_benchmark.real_data.go_invisible_gate
--analysis brainseq-sczd` (SLURM: `real_data/brainseq/_h/12.go_invisible_gate.sh`).
Outputs under `real_data/brainseq/caudate_sczd/_m/`: `go_invisible_gate.parquet`,
`GO_INVISIBLE_GATE.md`, `go_invisible_gate_background.json`.

**Result (res 5.0).** 8 disease-sig SCZD switch modules; 2 GO-enriched (M013/M015 =
clean immune pathways), **6 GO-invisible** (M010/M017/M020/M021/M023/M026). The
GO-invisible modules carry **real isoform switches** (nearly all members have
anticorrelated transcript pairs; switch strength 1.0–1.4) that are functionally
consequential (CDS/coding-status/biotype/UTR changes) **at the same rate as the
GO-visible disease modules and the pooled background** — structurally
indistinguishable. Drivers are plausible + psychiatric-relevant (GNAL, PRKCB,
FRMPD4, STXBP5, DDX3X, DGKH…) but heterogeneous within a module (shared switch
axis, not a shared GO process). **Verdict: PASS in the complementary form** — a
genuine DTU-without-DGE layer invisible to pathway enrichment because the signal is
isoform regulation, not a shared GO term. Frame on mechanism; do NOT claim "pathways
WGCNA misses". See `memory/project_go_invisible_gate.md`.

### 1b. MECHANISM — GTEx sQTL/eQTL genetic anchoring of co-switch modules — IMPLEMENTED

Reproducible CLI `python -m isograph_benchmark.real_data.qtl_anchoring
--analysis brainseq-sczd` (SLURM array: `real_data/brainseq/_h/13.qtl_anchoring.sh`,
17 analyses = SCZD + 3 brainseq aging + 13 gtex aging). Power-matched logistic
enrichment (qtl status ~ module membership + log cis-variant count + log gene length
+ log isoform count [+ log intron group size]) within each xQTL tested universe ∩
IsoGraph genes. Data: `inputs/raw/gtex_v11/xqtl/` (sGenes/eGenes per brain tissue).
Outputs `<_m>/{qtl_anchoring.parquet, QTL_ANCHORING.md, qtl_anchoring.json}`.

**Cross-tissue meta DONE** (`qtl_anchoring_meta.py`; 17 analyses; IVW fixed + DL
random effects; outputs `real_data/_m/qtl_anchoring_meta/`). Pooling **reversed** the
favorable single-tissue SCZD tail (sQTL OR 1.37) — read it honestly:
- Co-switch genes are cis-QTL **depleted for both** sQTL and eQTL (coordinated/network
  genes are constrained → fewer common cis-QTL). This shared baseline is not the result.
- **The result is the paired splicing-specificity contrast** (sQTL OR / eQTL OR within
  analysis, removing the constraint baseline): all 1.07 (p=8e-6), pheno-sig 1.13
  (p=1.5e-4), **GO-invisible 1.13 (p=3e-3, I²=0.15 — consistent across tissues)**,
  GO-visible 1.04 (**ns**). Splicing-QTL is spared ~7–13% vs expression-QTL in
  co-switch genes, significant exactly for the phenotype-associated + GO-invisible
  modules and **null for the GO-visible (immune/abundance) modules** — a clean internal
  control. Splicing-specific genetic anchoring concentrates where IsoGraph's unique
  value is (DTU-without-DGE).
- Manuscript line: frame on the **sQTL-vs-eQTL specificity contrast**, not raw sQTL
  enrichment. Scope caveat (in report): cis-sQTL anchors member-gene splicing to
  genetics, not the co-switching coordination itself. Optional secondary:
  switch-transcript→LeafCutter-intron direction concordance.

**Matched-baseline control DONE** (`qtl_anchoring.py --method {wgcna_switch_only,
wgcna_multiplex}`; SLURM `real_data/gtex/_h/08.qtl_anchoring_matched.sh`, 2 methods × 13
GTEx tissues; `qtl_anchoring_meta.py` now multi-method). Anchors the matched WGCNA
baselines (same switch features) and compares the splicing-specificity contrast on the
**same 13 GTEx tissues**. Result is a clean **method effect**: only IsoGraph shows it —
go_invisible **1.11 (p=0.011)**, pheno_sig 1.11 (p=2e-3), go_visible 1.02 (ns, internal
control); `wgcna_switch_only` and `wgcna_multiplex` are **null everywhere** (go_invisible
1.02 / 0.99, ns). Same features + classical inference loses the splicing-genetic signal
that IsoGraph's VAE+Leiden concentrates in GO-invisible modules. (NB: a single tissue —
frontal_cortex — looked specific for wgcna_switch_only at 1.29, but it does not survive
pooling; trust the meta, not one tissue.) Complements item 2: WGCNA-switch matches
IsoGraph on phenotype-sig RATE but NOT on genetic splicing-specificity.

### 2. Three-baseline comparison synthesis — DONE

Reproducible CLI `python -m isograph_benchmark.real_data.baseline_comparison`
(login-node aggregation; no SLURM). Outputs `real_data/_m/baseline_comparison/`:
`baseline_comparison.parquet` (per region×method), `baseline_comparison_pooled.parquet`
(per method), `BASELINE_COMPARISON.md`. 17 analyses, 16 with all four methods.

**Result — read on per-module RATES, not totals** (totals scale with module count;
IsoGraph runs finer: ~35 vs 8–18 modules). Pooled mean per-module rates:
- pheno-sig rate: wgcna_switch_only **0.336** > isograph 0.268 > wgcna_multiplex
  0.189 ≈ wgcna_gene 0.180. **Phenotype signal lives in the SWITCH features**, not the
  method — both switch-fed methods beat both abundance-fed ones.
- BOTH (pheno-sig AND GO) rate: isograph is **LOWEST** (0.074) — its GO-enrichment is
  low by construction. GO rate: wgcna_gene 0.885 > wgcna_multiplex 0.758 >>
  wgcna_switch_only 0.381 > isograph 0.217 (abundance/GO bias, expected).
- **IsoGraph is NOT globally superior on module-level metrics**; classical
  wgcna_switch_only matches/beats its phenotype rate. The **one clean method effect**:
  on IDENTICAL switch+abundance features, isograph 0.268 > wgcna_multiplex 0.189 — VAE+
  Leiden recovers phenotype-linked switch structure that classical multiplex WGCNA
  dilutes back toward abundance.
- Replication (isograph vs wgcna_gene only): classical WGCNA carries more cross-cohort
  GO overlap; both beat their perm null.
- **Manuscript line:** do NOT claim IsoGraph beats WGCNA on enrichment/phenotype rates.
  Its defensible value is the DTU-without-DGE **content** (biology gate + incremental
  association + sQTL/eQTL specificity), invisible to any abundance pipeline — complementary
  layer, consistent with the de-confounded gene-level result.

### 3. Real-data confounds ablation (NM gap #6) — DONE (synthetic); #30 QC landed

Synthetic data is clean; real bulk brain RNA-seq is not. The headline this
enables: **IsoGraph WITH vs WITHOUT `residualize_covariates` under confounded
data**, proving residualization buys robustness vs WGCNA. Four canonical
confounds (zero-inflation is NOT one — bulk is NB overdispersion, already modeled):

1. RNA degradation (RIN-like 3′ coverage bias) — directly distorts PSI / transcript
   ratios = IsoGraph's signal.
2. Cell-type composition variation (neuron/glia mixture; shifts with age).
3. Batch effects (flowcell / processing batch).
4. Library depth / mapping-rate variation.

**Acceptance MET (synthetic branch):** `figS8_confound_robustness` is the ablation
figure — under composition / batch / depth confounds, `isograph_vae_residual` stays
flat while raw VAE collapses and WGCNA degrades; RNA degradation is the honest limit
(RIN regression can't recover a coordinate that isn't a single recorded axis), and
`figS9_degradation_fallback` shows the abundance-channel variants (multiplex /
reliability) recover it. Both rendered from completed runs in
`benchmark/03_metrics/figures/`.

**#30 DONE (2026-06-27):** A/B QC metrics in feature residualization. `residualize.py`
gains `residualization_qc(before, after, design, feature_info)` — per-feature
`var_retained_frac` + `confound_r2_before/after`; the working residualization drives
`confound_r2_after`→0 while `var_retained_frac` reports the collateral signal cost.
This is the diagnostic observable on REAL data (where module recovery is not). Wired
through `FitArtifacts.residualization_qc` (vae.py captures it in both branches) and
written as `residualization_qc.parquet` by `run_models.py`. Production fits already
set `residualize_covariates`, so it populates on every real fit. Unit + e2e tested
(`tests/test_residualization_qc.py`). Numerically non-invasive — does not change the
residualized features. UNCOMMITTED.

**Covariate policy DECIDED (2026-06-28): residualization is a DISCOVERY knob only.**
Covariates do two jobs: protect module *discovery* (a confound axis fabricates artifact
modules — no downstream test can un-merge them; synthetic ablation proves it), and
adjust trait *inference* (best done jointly with the trait — spline for age, conditioning
on the abundance channel). These are split, not stacked. Previously `feature_scores.parquet`
held the *residualized* matrix (vae.py:660 used the in-place-residualized `switch_matrix`;
the line-665 comment wrongly said "raw"), so `incremental_association` adjusted the SAME
covariates a second time = double residualization. **Fix (uncommitted):** vae.py now keeps
`switch_matrix_raw` and persists feature_scores from it — clustering/VAE still embed the
residualized matrix, but feature_scores is RAW, so the downstream test is the single
inference adjustment. Verified e2e: feature_scores retains the covariate axis (r²≈0.80)
while the embedded matrix has confound_r²→0. `modules.parquet` is UNCHANGED (clustering
input identical, deterministic seed); only feature_scores + downstream (incremental_association,
trait_associations, module_gene_roles, characterize) change → production real-data re-fit
needed. Task #31 (silent-skip guard on `build_design_matrix`) and #32 (finalize upstream
list as technical-confounds-only, exclude trait + biological Sex/MoD; launch re-fit) follow.

- **#10 remains the only open degradation task** (module-level QC flag + edge-type
  penalties); **do NOT start it without explicit confirmation** — it touches the
  feature residualization path.

### 4. Per-module trust funnel writeup

The funnel goal is **per-module trust**, not global ARI (global ARI is
deprioritized — finer res-5.0 isograph partitions are inherently less
split-half-stable, which is expected and not the headline). Assemble the funnel:
stability → drivers (`module_meta` loadings, done) → aging replication
(`replication`) → WGCNA-complementarity.

- For each trusted module, chain: split-half co-assignment stability → reproducible
  switch-axis driver loadings (ρ≈0.77 shared-gene) → cross-cohort aging replication
  → is it complementary to (not redundant with) WGCNA.
- **Acceptance:** a short list of high-trust modules that survive all four gates,
  carried forward as the modules the biological claims rest on.

### 5. Manuscript summary tables + figures

After 1–4 land, produce publication artifacts.

- Use the `per-analsis-summarization` skill for analysis→Manubot summary tables.
- Use the `manuscript-figures` skill for the real-data figures (three-baseline
  comparison, trust funnel, GWAS res5 specificity, confound ablation).
- **Acceptance:** main + supplemental real-data figures and Table(s) drafted,
  each tied to a specific honest claim above.

---

## Pointers

- Stability / trust harness: `real_data/stability/` (+ `MODULE_TRUST_PLAN.md`,
  `SOFTWARE_ROBUSTNESS_PLAN.md`).
- Real-data analysis modules: `isograph_benchmark/real_data/` (`interpret_modules.py`,
  `module_enrichment.py`, `replication_go.py`, `incremental_association.py`,
  `characterize_composition_unique.py`, `run_matched_wgcna.py`, `run_models.py`).
- SLURM launchers live under each cohort's `_h/` (brainseq) and per-region trees.
- Memory index: `~/.claude/.../memory/MEMORY.md` — the project_* files carry the
  detailed state (resolution5_comparison, wgcna_rigor_fixes, strengths_limitations,
  module_trust_funnel, nature_methods_gaps).
