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

### 2. Three-baseline comparison synthesis

The enrichment + replication_go outputs now exist for all four methods but are
**not yet synthesized**. Assemble the head-to-head.

- Pull per-method `module_enrichment` tables (GO-enriched %, pheno-sig count,
  BOTH count, network metrics) for isograph vs the 3 WGCNA baselines, per region.
- Pull `replication_go` cross-cohort GO consistency across the 4 methods.
- **Key question the matched baselines answer:** does giving WGCNA the *same*
  switch/multiplex feature matrix close IsoGraph's gap, or is the network
  inference (VAE + Leiden) doing real work beyond the input representation? This
  isolates method from input — the cleanest possible ablation for reviewers.
- **Acceptance:** one comparison table + narrative stating exactly where IsoGraph
  wins, where it ties, where WGCNA wins — with the matched baselines settling the
  "is it the features or the method?" question.

### 3. Real-data confounds ablation (NM gap #6, still open)

Synthetic data is clean; real bulk brain RNA-seq is not. The headline this
enables: **IsoGraph WITH vs WITHOUT `residualize_covariates` under confounded
data**, proving residualization buys robustness vs WGCNA. Four canonical
confounds, in priority order (zero-inflation is NOT one — bulk is NB
overdispersion, already modeled):

1. RNA degradation (RIN-like 3′ coverage bias) — top priority; directly distorts
   PSI / transcript ratios = IsoGraph's signal. (Partial infra exists: TIN
   reliability + gene-body coverage covariate; tasks #10/#30 remain — module-level
   QC flag + edge-type penalties, A/B QC metrics in feature residualization.)
2. Cell-type composition variation (neuron/glia mixture; shifts with age).
3. Batch effects (flowcell / processing batch).
4. Library depth / mapping-rate variation.

- **Do NOT start tasks #10/#30 without explicit confirmation** — they touch the
  feature residualization path.
- **Acceptance:** an ablation figure showing IsoGraph's metric degrades far less
  than WGCNA when residualization is on under each confound; degradation-robust
  behavior demonstrated on the real (or realistically-confounded synthetic) data.

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
