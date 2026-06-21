# IsoGraph software-robustness plan: module collapse + VAE stability

**Status:** design doc + first implementation increment landed (2026-06-21, IsoGraph repo
branch `consensus-clustering`, uncommitted). Companion to `MODULE_TRUST_PLAN.md` (which
validates *modules*); this doc hardens the *software* that produces them, for Nature Methods
reproducibility-by-strangers readiness.

**Implemented so far (behind flags, legacy defaults, unit-tested):**
- VAE **S1 gradient clipping** — `VaeModelConfig.grad_clip_norm` (None=off); clips global grad
  norm before each optimizer step (`models/vae.py:_train_single_vae`).
- VAE **S2-lite divergence guard** — non-finite val-loss -> restore best + stop (always on,
  only reachable on the pathological path; smoke-tested: lr=50 diverges gracefully, lr=50+clip
  stays finite, normal lr unaffected).
- Collapse **C-determinism** — single-run Leiden is now edge-weighted + seeded
  (`models/base.py:_detect_communities`); removes the unweighted/unseeded variance source.
- Collapse **E giant-module cap** — `VaeModelConfig.max_module_frac` (None=off); any community
  over the cap is recursively re-clustered at escalating resolution
  (`_split_oversized_communities`/`_recursive_split`; unit-tested: 100-node giant -> 4x25,
  deterministic, node-conserving; correctly leaves true cliques irreducible).

**Pending:** S2 full auto-LR-backoff+restart, S3 free-bits, S4 graph-density cap, S5 promote
estimability to production, B two-axis gene-similarity centering, D adjacency soft-threshold;
**and the validation step** — wire `leiden_resolution`+`max_module_frac`(+`grad_clip_norm` for
GTEx) into the stability harness / `run_models`, then re-run the split-half A/B across all 6
regions (ARI/NMI vs estimability baseline; giant-fraction within cap; synthetic + Type-I clean).

**Cleanup task — remove differential-TIN switch-reliability source.** Differential-TIN was a
*negative* lever in the split-half A/B (regressed ARI/NMI alongside consensus-Leiden and
median-TIN; only estimability improved stability). It currently rides along inside the
estimability commit on `keeper-trunk` (`switch_reliability_source="differential_tin"` + its
supporting functions in `reliability.py`/`switch.py`). Excise the differential-TIN code path
and its helpers so only the estimability source ships; keep `estimability` as the sole
promoted reliability source (§2.2 S5). Verify nothing references `differential_tin` after
removal and that the estimability path is unaffected.

## 0. Why this is the gating software work for NM

A tools-journal paper lives or dies on a stranger reproducing it on their own data. Two
defects currently block that bar, and both are *core IsoGraph* issues, not benchmark glue:

1. **Module collapse** — a giant catch-all module (M000) recurrently absorbs hundreds–
   thousands of genes. Production caudate M000 is 0.20-effect/p=0.002 age-associated *and*
   the single largest module: a reviewer running IsoGraph will hit it immediately.
2. **VAE-fit fragility** — the default `lr=1e-3` diverges on GTEx (val ELBO→~1e8/nan); only
   the hand-tuned `lr=3e-4` is stable, and divergence is **data-geometry/init-dependent,
   not sample-size-driven** (`nucleus_accumbens` n=285 diverged while `substantia_nigra`
   n=183 was stable). A method that needs per-dataset LR babysitting is a reproducibility
   red flag.

What this work is **not**: a fix for cross-cohort module *correspondence*. The pooled
replication test (`module_trust.py replication-pooled`) found IsoGraph best-match module
gene-set Jaccard pools to median 0.034 (max 0.21) vs WGCNA's median 0.111 (22% of modules at
J≥0.25). That gap is **largely a partition-granularity artifact, not a quality defect**:
WGCNA makes a handful of giant modules (5–7 per GTEx region, ~2–3k genes each) while IsoGraph
makes many small ones (14–24); best-match Jaccard mechanically rewards coarse partitions. The
collapse fix (finer modules) will push IsoGraph's module-Jaccard *down*, not up. So
cross-cohort best-match module Jaccard is a confounded metric and must not be used as a
pass/fail gate (see §5). The collapse and stability fixes below are justified on their own
merits — the giant M000, within-cohort split-half stability, and stranger-reproducibility —
**not** as cross-cohort-correspondence improvers.

**Objective functions (do not regress):**
- **Synthetic** module recovery / switch-gene detection (`benchmark/.../synthetic_*`).
- **Within-cohort split-half** ARI/NMI + `n_common` (`stability.py`, the committed A/B dial).
  Estimability is the *only* lever that has improved this so far (mean ΔARI +0.008, ΔNMI
  +0.012, 6/6 regions); consensus-Leiden, differential-TIN, and median-TIN all regressed it.
- **Type-I:** negative-control / non-switching-background scenarios must keep ~0 false
  module recovery (collapse fixes must not shatter the null into spurious small modules).

---

## 1. Module collapse

### 1.1 Root cause (confirmed)
`_module_table` (identical in `vae.py:276`, `baseline.py:40`, `latent.py:88`, `graph.py:79`)
builds modules with **`nx.connected_components`**. Connected components merges any two
communities that share a *single* bridge gene, so on a dense reconstruction-similarity graph
(~354k positive edges, avg degree ~106) everything fuses into one giant component. This is a
graph-algorithm defect, not a signal defect.

The fix is proven but only lives in a **benchmark-side repair script**
(`repair_caudate_sczd.py`: Leiden at `resolution=10.0, seed=13` → DRD2 in a focused
~180-gene module vs the ~2,558-gene M000 at the old connected-components/res=2.0). It must be
lifted into IsoGraph core so *every* fit benefits, not just the patched case study.

| Resolution | giant size | DRD2 module | n modules ≥20 |
|---|---|---|---|
| 2.0 (≈conn-comp) | ~1931 | ~2558 | 17 |
| 5.0 | ~748 | ~690 | 25 |
| 10.0 | ~524 | ~180 | 39 |
| 20.0 | ~444 | ~69 | 35 |

### 1.2 Interventions (ranked)

**C — Replace connected-components with seeded, edge-weighted Leiden in core `_module_table`
(highest leverage; do first).**
- Add `leiden_resolution: float | None = None` and `leiden_seed: int = 13` to
  `VaeModelConfig` (`workflow/config.py`). `None` ⇒ legacy connected-components (preserves
  old run hashes); a float ⇒ Leiden.
- Implement once in a shared helper (the four `_module_table` copies are identical — collapse
  them into one `isograph.models._modules.partition(graph, cfg)` to avoid 4-way drift).
- Weight Leiden by edge similarity; fix the seed for determinism.
- **Resolution must be data-driven, not hard-coded 10.0.** Select by stability: reuse
  `sweep_leiden.py` to pick the resolution that maximizes within-cohort split-half ARI
  subject to a giant-fraction cap (largest module ≤ X% of assigned genes, X~15–20). Persist
  the chosen resolution in run metadata.

**B — Two-axis gene-similarity centering (`vae.py:_gene_similarity`, line 264).**
- Currently Pearson on the VAE reconstruction centered by gene (axis=1) only. Add axis=0
  (across-gene) centering so a gene that is globally high/low doesn't co-rank with everything
  → fewer bridge edges feeding the giant component. Cheap; compose with C.

**D — Soft-threshold / shrink the switch adjacency before clustering.**
- Denoise low-SNR edges (the other documented variance source) by soft-thresholding
  similarities before Leiden, so weak bridge edges don't survive. Tunable on the same
  stability objective; guardrail against over-pruning `n_common`.

**E — Post-hoc giant-module split (fallback).**
- If a module still exceeds the giant cap after C+B+D, recursively re-cluster *within* it.
  Keep as a safety net, not the primary mechanism.

### 1.3 Acceptance criteria
- Largest module ≤ ~15–20% of assigned genes on **all 6 trust-funnel regions** (today M000
  violates this).
- Within-cohort split-half ARI/NMI **not worse** than the estimability baseline on any region.
- Synthetic module recovery unchanged or better; negative-control false recovery stays ~0.
- DRD2 case study reproduces from a *core* fit with no benchmark repair script.

---

## 2. VAE-fit stability

### 2.1 Root cause (confirmed)
Divergence at `lr=1e-3` is data-geometry/init-dependent (not n-driven); the degenerate
reconstruction then explodes the similarity graph into a near-complete graph and OOMs the
module step (64 GB). The current mitigation is a hand-picked `lr=3e-4` for GTEx and a
different `lr=1e-3` for BrainSEQ — exactly the per-dataset tuning a stranger can't know.

### 2.2 Interventions (ranked) — goal: one default trains every dataset without divergence

**S1 — Gradient clipping + LR warmup (do first, cheapest).**
- Global-norm gradient clipping (e.g. max-norm 1–5) and a short linear warmup. Most ELBO→1e8
  blow-ups are a few early exploding steps; clipping alone often removes the divergence cliff
  and lets a single LR serve all regions.

**S2 — Divergence guard with automatic LR backoff + restart.**
- Detect non-finite / runaway val-ELBO during training; on trip, restore the last good
  checkpoint, halve LR, and resume (bounded retries). Makes the fit *self-stabilizing*: the
  user never picks an LR. Log the backoff so it's auditable.

**S3 — KL handling: warmup/annealing + free-bits.**
- Anneal the KL weight and add a free-bits floor so posterior-collapse/term-imbalance (a
  common ELBO-explosion driver) can't drive the optimizer off a cliff.

**S4 — Graph-explosion safeguard at the boundary.**
- Even with S1–S3, cap reconstruction-graph density (top-k neighbors or hard similarity floor)
  so a bad fit degrades to "poor modules," never a 64 GB OOM. Decouples the failure mode from
  the memory wall and protects the cluster.

**S5 — Promote estimability to production.**
- `switch_reliability_source="estimability"` (`switch_estimability_min_minor_usage=0.1`) is
  the only committed *positive* reproducibility lever (IsoGraph 7b6e945) but still lives only
  in the stability harness — `run_models` production fits don't use it. Promote it (behind the
  config flag, default on for real-data) so the shipped pipeline gets the +ARI/+NMI it earns.
  Re-confirm synthetic + Type-I are unaffected before flipping the default.

### 2.3 Acceptance criteria
- A **single** documented LR/optimizer config trains all 6 trust-funnel regions *and* the
  GTEx region that diverged (`nucleus_accumbens`) to BrainSEQ-range RMSE (~1.04–1.08) with no
  divergence and no manual per-region tuning.
- No fit produces a near-complete similarity graph or OOMs at the standard fit memory
  (~48 GB ceiling; see `project_isograph_fit_memory`).
- Estimability-on production fits match or beat current production on synthetic recovery and
  within-cohort stability; Type-I unchanged.

---

## 3. Sequencing

1. **S1 + S4** (clipping/warmup + graph-density cap) — removes the divergence/OOM cliff;
   unblocks running everything on one config. *(IsoGraph core.)*
2. **C + B** (Leiden in core `_module_table` + two-axis centering) — kills the giant module
   at the source. *(IsoGraph core; collapse the 4 `_module_table` copies into one helper.)*
3. **Re-run the stability A/B** across all 6 regions; confirm ARI/NMI ≥ estimability baseline,
   giant-fraction within cap, synthetic + Type-I clean.
4. **S5** (promote estimability) + **S2/S3/D/E** as needed if (3) shows residual fragility.
5. Do **not** expect (or target) higher cross-cohort module Jaccard — finer modules lower it.
   The cross-cohort replication question moves to a granularity-invariant design
   (gene/transcript-level replication or gene-pair co-assignment agreement; see §5 and the
   Q2/Q3 reframe), where IsoGraph and WGCNA are compared at matched altitude.

## 4. Risks / honest caveats
- **Leiden resolution is a new knob.** Selecting it by stability (not hard-coding 10.0) is
  essential or we trade the giant module for an over-fragmented one — and over-fragmentation
  *lowered* ARI in the consensus-Leiden A/B. The giant-fraction cap + ARI objective must be
  enforced jointly.
- **Determinism.** Seed Leiden and any restart logic; a "self-stabilizing" fit must still be
  bit-reproducible from `(config, seed)` or it fails the NM bar it's meant to clear.
- **Some cross-cohort non-correspondence is irreducible** (different quantifiers); detection
  fixes will reduce, not eliminate, it. Do not over-promise module-level transfer.
- **Scope discipline.** This is core-IsoGraph work (repo `software/IsoGraph`,
  branch `consensus-clustering`); changes ship behind config flags with legacy defaults so
  existing run hashes and the committed benchmark results are preserved until each lever is
  validated and deliberately promoted.

## 5. Cross-cohort metric correction (evidence, 2026-06-21)

Best-match cross-cohort module gene-set Jaccard, restricted to shared genes, modules ≥10
genes, pooled over the 3 homologous region pairs:

| method | median best-J | max | frac ≥0.25 | frac ≥0.10 | module counts (BS→GTEx) |
|---|---|---|---|---|---|
| IsoGraph | 0.034 | 0.213 | 0% | 5% | 15–33 → 14–24 |
| WGCNA | 0.111 | 0.780 | 22% | 54% | 9–24 → **5–7** |

**Reading:** independently-built networks *can* share cross-cohort structure (WGCNA does), so
"separate networks can't match" is false. But the IsoGraph–WGCNA gap is **mostly a granularity
confound** — best-match Jaccard rewards WGCNA's handful of giant modules. Therefore:

- **Retire best-match module Jaccard as a cross-cohort quality metric** (keep it only as a
  descriptive, granularity-noted statistic). It is unfair to finer partitions and would be
  *worsened* by the §1 collapse fix.
- **Replace with granularity-invariant cross-cohort tests:** (a) **gene/transcript-level
  replication** — does the discovery-significant aging gene/transcript *set* replicate in the
  other cohort (matching-free meta-analysis; the `replication-pooled` Stouffer machinery
  applied at gene level instead of module level); (b) **gene-pair co-assignment agreement** —
  Rand/ARI-style, over gene pairs present in both cohorts, optionally at matched resolution.
- **Honest limitation to keep in the manuscript:** switch covariation is intrinsically less
  cross-cohort-reproducible than abundance co-expression (consistent with the committed
  within-cohort diagnosis in `real_data/README.md`). IsoGraph's contribution is the
  DTU-without-DGE layer, not module-reproducibility parity with WGCNA — do not frame it as the
  latter.
