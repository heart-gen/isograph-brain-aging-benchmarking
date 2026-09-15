# 04 — Module trust

**Question:** are the discovered modules reproducible, or is a fine-grained partition
just noise?

Two arms that answer one question, so they share a stage: within-cohort split-half
stability, and cross-cohort (BrainSEQ ↔ GTEx) aging replication. BrainSEQ and GTEx are
independent cohorts — different donors, libraries and quantifiers — run through the same
pipeline, which is what makes the second arm a real replication.

Stage order: this stage was `03_module_trust` until 2026-09-15. It now follows module
characterization, because the replication, permutation, pooled and complementarity steps read
stage 03's interpretation, enrichment and composition-unique genes. The three-baseline rate
comparison moved here from characterization for the same reason in reverse: it reads this
stage's replication tables. Functional preservation reads stage 06 and is now
`08_integration/_h/01c`.

## Order

Run the whole stage with `bash 04_module_trust/_h/run_stage.sh` (add `--dry-run` to print the
plan). The leading number of a wrapper is its tier; steps in one tier run in parallel.

| Step | Wrapper | Waits on | Produces |
|---|---|---|---|
| 01a–01b | `stability_isograph`, `stability_wgcna` | stage 02 | Split-half partitions; IsoGraph run with `STABILITY_SAVE_EDGES=1` so 02a can re-cluster |
| 01c | `replication` | stage 02 | Cross-cohort module matching |
| 01d | `age_model_curvature` | stage 02 | Curvature test licensing the linear age model |
| 02a | `stability_resolution_sweep` | 01a | Phenotype-blind Leiden resolution curve (re-clusters saved split-half graphs) — a disclosed sensitivity; 5.0 stays on the ≥ 900-gene criterion (PI decision 2026-09-12) |
| 02b | `module_meta` | 01a–01b | Per-module meta tables, driver loadings, one job per region × method (re-run after any split-half re-fit, then 03b) |
| 02c | `trust_gate` | 01a–01b | Q1 trust gate (`module_trust stability`): the trusted sets read by 03c–03g |
| 02d | `replication_go` | 01c, stage 03 enrichment | GO consistency of matched modules |
| 03a | `stability_aggregate` | 01a–01b, 02a | `stability_summary`, one row per method and resolution |
| 03b | `within_cohort` | 02b | Within-cohort Q2 driver reproducibility + Q3 sign concordance (`module_trust within`), per method |
| 03c | `module_trust_replication` | 02c, stage 03 interpretation | Q3 cross-cohort aging replication of trusted modules, linear and `--model spline` |
| 03d | `replication_permutation` | 02c | Permutation null, one array per `--covariates {full,complement,none}` |
| 03e | `replication_pooled` | 02c | Q3 pooled over all region pairs (Stouffer meta-Z + permutation arms) |
| 03f | `eigengene_projection` | 02c | Frozen-weight cross-cohort eigengene projection |
| 03g | `complementarity` | 02c, stage 03 interpretation + composition-unique | Q4 DTU-without-DGE / WGCNA complementarity of trusted modules |
| 03h | `baseline_comparison` | 02d, stage 03 enrichment | Three-baseline per-module rate comparison (the scope bound) |
| 04a | `replication_model_contrast` | 03c | Linear-vs-spline decomposition of the Q3 count |
| 04b | `replication_permutation_report` | 03d | `REPLICATION_PERMUTATION.md` over all three covariate modes |
| 04c | `eigengene_projection_aggregate` | 03f | Cross-pair eigengene-projection summary |

Retired (`_h/retired/`): `lr_validation`, `lr_validation_launch`, `lr_aggregate` (learning-rate /
software-robustness validation, retired 2026-09-14 — settled the single-LR promotion; outputs only
at tag `legacy_expression_filter`) and `gcap_ab` (giant-cap ablation, retired 2026-09-12 — a
resolution-2.0 A/B of a cap never promoted, superseded by resolution 5.0 and cited nowhere).

`_m/` holds `stability/`, `replication/` and `baseline_comparison/`. The resolution-2.0 siblings
(`stability_res2`, `replication_res2`, which nothing read or wrote), the giant-cap ablation and the
LR validation outputs were removed from the live tree on 2026-09-14 and survive only at tag
`legacy_expression_filter`.

**CLIs:** `isograph_benchmark/real_data/{stability,module_trust,replication,replication_go,replication_permutation,age_model_curvature,eigengene_projection,baseline_comparison}.py`.

## Caveat

Cross-cohort **module** Jaccard is granularity-confounded — a method with few giant
modules wins it mechanically — so it is retired as a quality gate in favour of
gene-level comparison. See `docs/MODULE_TRUST_PLAN.md`.

It is retired *as a gate*, not as a reported quantity, and step 03e is the cautionary
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
objection. Step 02a answers it on split-half stability instead:

    sbatch --export=ALL,STABILITY_SAVE_EDGES=1 04_module_trust/_h/01a.stability_isograph.sh
    sbatch 04_module_trust/_h/02a.stability_resolution_sweep.sh
    sbatch 04_module_trust/_h/03a.stability_aggregate.sh

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

Reproduce with `02a.stability_resolution_sweep.sh` then `03a.stability_aggregate.sh`; the
per-resolution rows live in `stability_summary.parquet` under method `isograph_resXpY`.

## Display items

Main **Fig 2** `figTrustFunnel`; supplementary `tableS7_module_trust_funnel`; S-real-1
`figBaselineRates` and tables S1/S2 from the three-baseline comparison (03h).
Design notes: `docs/MODULE_TRUST_PLAN.md`, `docs/SOFTWARE_ROBUSTNESS_PLAN.md`.

## Age model: why the linear arm is quoted

The cross-cohort aging replication count is **23/130 matched module pairs, `complement`
covariate mode, linear age model** (matching permutation p = 0.034). The df = 3 spline
arm gives 15/130 (p = 0.146) at the same covariate mode.

That gap is a power cost, not a correction, and `01d.age_model_curvature.sh` is the test
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

    sbatch 04_module_trust/_h/01d.age_model_curvature.sh
    # -> _m/stability/module_trust/age_model_curvature__isograph.{parquet,json}
