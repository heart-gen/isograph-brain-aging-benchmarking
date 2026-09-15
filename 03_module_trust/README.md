# 03 — Module trust

**Question:** are the discovered modules reproducible, or is a fine-grained partition
just noise?

Two arms that answer one question, so they share a stage: within-cohort split-half
stability, and cross-cohort (BrainSEQ ↔ GTEx) aging replication. BrainSEQ and GTEx are
independent cohorts — different donors, libraries and quantifiers — run through the same
pipeline, which is what makes the second arm a real replication.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–03 | `stability_isograph`, `stability_wgcna`, `stability_aggregate` | Split-half partitions + `stability_summary` |
| 04 | `module_meta` | Per-module meta tables, driver loadings (per split half — re-run after any split-half re-fit, then step 17) |
| 05–07 | `lr_validation`, `lr_validation_launch`, `lr_aggregate` | Learning-rate / software-robustness validation — **retired 2026-09-14**: settled the single-LR promotion; not re-run on the switching filter, outputs only at tag `legacy_expression_filter` |
| 08 | `gcap_ab` | Giant-cap ablation (`_m/stability_gcap_ab/`) — **retired 2026-09-12**: a resolution-2.0 A/B of a cap never promoted, superseded by resolution 5.0 and cited nowhere |
| 09 | `module_trust_replication` | Q3 cross-cohort aging replication of trusted modules |
| 10–11 | `replication`, `replication_go` | Cross-cohort module matching + GO consistency |
| 12–13 | `replication_permutation`, `replication_functional` | Permutation null + functional preservation |
| 14 | `age_model_curvature` | Curvature test licensing the linear age model |
| 15 | `replication_pooled` | Q3 pooled over all region pairs (Stouffer meta-Z + permutation arms) |
| 16 | `stability_resolution_sweep` | Phenotype-blind Leiden resolution curve (re-clusters saved split-half graphs) — a disclosed sensitivity; 5.0 stays on the ≥ 900-gene criterion (PI decision 2026-09-12) |
| 17 | `within_cohort` | Within-cohort Q2 driver reproducibility + Q3 sign concordance (`module_trust within`); run after step 04. It had no launcher, which is how a stale `modules_meta` went unnoticed after the 2026-09-09 re-fit |

`_m/` holds `stability/` and `replication/`. The resolution-2.0 siblings (`stability_res2`,
`replication_res2`, which nothing read or wrote), the giant-cap ablation and the LR validation
outputs were removed from the live tree on 2026-09-14 and survive only at tag
`legacy_expression_filter`.

**CLIs:** `isograph_benchmark/real_data/{stability,module_trust,replication,replication_go,replication_permutation,replication_functional,age_model_curvature}.py`.

## Caveat

Cross-cohort **module** Jaccard is granularity-confounded — a method with few giant
modules wins it mechanically — so it is retired as a quality gate in favour of
gene-level comparison. See `docs/MODULE_TRUST_PLAN.md`.

It is retired *as a gate*, not as a reported quantity, and step 15 is the cautionary
example. Its `--min-jaccard` defaulted to 0.25, which admits WGCNA's giant modules
(max 0.65 here) and excludes IsoGraph's fine ones (max 0.12) — so both pooled tables sat
at **0 rows** for months, read by anything that opened them as "nothing replicates". The
default is now 0.0 and a zero result writes no parquet at all, only a stats json naming
the gate. If you add a Jaccard threshold anywhere, check first whether it is measuring
replication or measuring module size.

**The pooled arm is null for both methods** (IsoGraph Stouffer Z = −1.15, perm p = 0.122,
sign concordance 68/130 p = 0.066, Spearman ρ = 0.108 p = 0.110; WGCNA Z = −1.16,
p = 0.079, 27/53 p = 0.086). These are the post-2026-09-09 split-half re-fit values; the
earlier Z = −0.68 / p = 0.113 and its nominal 69/130 sign test are superseded. It does not
corroborate the per-region 23-vs-3 count, which is a different estimand — see the Q3-pooled section of
`_m/stability/module_trust/MODULE_TRUST_SUMMARY.md` for why the quantifier gap separates
them.

## The split-half partitions were a generation behind (fixed 2026-09-09)

`fit-isograph` writes the split halves that the trust funnel validates production modules
against, and `_vae_config` pins them to the production resolution *precisely so* the funnel
validates the modules the paper ships. The committed partitions did not honour that: they
carried **7–26 modules per region** — the `isograph_vae_res2` scale (10–18) — while
production is 39–50 and a fresh half at the canonical 5.0 gives 29–67. Half-sample size does
not explain it; the fresh halves have *more* modules than production, not fewer.
`brainseq/dlpfc` was a **mixture**: nine of its ten files were current and one (seed0 A, 9
modules) was not, which is the `out.exists()` resume-skip leaving stale files in place.

Re-fitting all 60 halves at 5.0 moved the funnel **in the project's favour**, which is worth
stating plainly since the correction could as easily have gone the other way:

| | before | after |
|---|---|---|
| Chance-trusted modules | 236/266 (89%) | **250/266 (94%)** |
| Split-half ARI (median over regions) | 0.21 | **0.35** |
| Cross-cohort aging replications | 25/130, perm p = 0.014 | **23/130, perm p = 0.034** |
| df = 3 spline arm | 15/130, p = 0.152 | 15/130, p = 0.146 |
| Pooled arm | null | null |

The WGCNA arm was untouched and needed no correction — its halves (5–33 modules) already
bracket its production fits (5–25).

## Choosing the resolution without looking at a trait

The canonical Leiden resolution (5.0) was selected using a phenotype-aware argument — the
≥900-gene giant module produced a GWAS artifact at 2.0 — which is a standing reviewer
objection. Step 16 answers it on split-half stability instead:

    sbatch --export=ALL,STABILITY_SAVE_EDGES=1 03_module_trust/_h/01.stability_isograph.sh
    sbatch 03_module_trust/_h/16.stability_resolution_sweep.sh
    sbatch 03_module_trust/_h/03.stability_aggregate.sh

The sweep **re-clusters saved graphs and refits nothing**. A gene-gene graph does not depend
on the Leiden resolution — only the clustering step does — so one fit pass with
`--save-edges` supports the whole grid, instead of one 30–48G fit pass per resolution. Each
non-canonical resolution is tagged `isograph_resXpY` and aggregates as its own method, so
the curve falls out of `stability_summary` and the committed canonical partitions are never
touched (the sweep skips 5.0 by design).

### Result (grid run 2026-09-09): split-half agreement cannot select a resolution

| resolution | mean ARI | mean NMI | mean modules | mean genes compared |
|---|---|---|---|---|
| 0.5 | **0.412** | 0.352 | 5.3 | 4,722 |
| 1.0 | 0.363 | 0.336 | 7.7 | 4,699 |
| 2.0 | 0.318 | 0.307 | 17.3 | 4,633 |
| 3.0 | 0.263 | 0.330 | 31.9 | 4,219 |
| **5.0 (canonical)** | 0.279 | 0.407 | 46.5 | 3,409 |
| 8.0 | 0.329 | 0.472 | 49.7 | 2,523 |
| 12.0 | 0.358 | 0.518 | 50.9 | 1,977 |
| 20.0 | 0.398 | **0.564** | 47.4 | 1,416 |

**Do not read the ARI column as a criterion.** It is U-shaped with a minimum at 3–5, so a
naive reading picks 0.5 — where there are **5.3 modules**. That is the granularity confound
this project already retired the cross-cohort Jaccard for (see the Caveat above), now in a
second guise: ARI is inflated at the coarse end by giant modules and at the fine end by a
collapsing comparison set, since the genes assigned to modules of size ≥ 20 in *both* halves
fall from 4,722 to 1,416 across the grid. The two ends are not comparable to each other, let
alone to the middle.

NMI, which is far less inflated by a few giant modules, moves the other way: it rises
monotonically through 5.0 and keeps rising. So the two agreement metrics disagree in
direction, and both track granularity rather than reproducibility.

**What this licenses us to say.** The phenotype-blind stability criterion the reviewer asks
for does not exist in the assumed form — split-half agreement does not have an interior
optimum to read a resolution off. It does not vindicate 5.0 and it does not condemn it: 5.0
sits at the ARI minimum but partway up a monotone NMI rise, still assigning 3,409 genes.
Report the curve, state why neither metric selects, and keep the ≥900-gene giant-module
criterion as the stated basis — now with the sweep as a disclosed sensitivity rather than an
absent analysis.

Reproduce with `16.stability_resolution_sweep.sh` then `03.stability_aggregate.sh`; the
per-resolution rows live in `stability_summary.parquet` under method `isograph_resXpY`.

## Display items

Main **Fig 2** `figTrustFunnel`; supplementary `tableS7_module_trust_funnel`.
Design notes: `docs/MODULE_TRUST_PLAN.md`, `docs/SOFTWARE_ROBUSTNESS_PLAN.md`.

## Age model: why the linear arm is quoted

The cross-cohort aging replication count is **23/130 matched module pairs, `complement`
covariate mode, linear age model** (matching permutation p = 0.034). The df = 3 spline
arm gives 15/130 (p = 0.146) at the same covariate mode.

That gap is a power cost, not a correction, and `14.age_model_curvature.sh` is the test
that says so. It applies a Wald second-difference contrast inside the spline's own fit —
is each module's projected age trajectory distinguishable from a straight line? Across all
266 modules in all six cohort x region fits, none is (0 at nominal p < 0.05 against 13.3
expected by chance; smallest p = 0.71). The contrast is well powered because the projected
effects share nearly all their uncertainty and every contrast row sums to zero, so the
shared part cancels. The spline also finds no *different* modules: 4/266 spline-only
significant against 78/266 linear-only, and 4 is the chance rate.

So aging is linear at module resolution here, the spline spends 2 df on noise, and the
spline arm belongs in the supplement as a declared sensitivity.

The test needs the full projection covariance (`cov_age_ij`) and refuses to run on the
diagonal alone — dropping the off-diagonal terms would inflate the variance by exactly the
shared term the contrast cancels, and so would manufacture a null. The WGCNA baseline's
`age_spline.parquet` carries only per-point standard errors, so it cannot be tested here.

    sbatch 03_module_trust/_h/14.age_model_curvature.sh
    # -> _m/stability/module_trust/age_model_curvature__isograph.{parquet,json}
