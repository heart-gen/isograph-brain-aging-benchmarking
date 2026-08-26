# AGENTS.md — Real-Data Rigor & Biology for Cell Genomics

Working plan for finishing the **real-data** analysis so it is rigorous and
biologically grounded for a **Cell Genomics** submission (retargeted 2026-07-16
from Nature Methods; rigor bar unchanged). This file is the to-do spine; it is
not a status report. Update it as items close.

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
- **Never `git add -A`.** Stage selectively. Commit only when explicitly asked.
  Do **not** add a `Co-Authored-By:` trailer to commit messages.
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

**Result (res 5.0, post covariate-decouple re-fit, regen 2026-06-29).** 4 disease-sig
SCZD switch modules (pheno_fdr ≤ 0.1), **all 4 GO-invisible** (M026/M020/M010/M023),
0 GO-enriched. (The earlier pre-refit run reported 8 modules with 2 GO-visible
M013/M015; the re-fit moved the phenotype-FDR landscape — cite the regenerated
parquet, not the old prose. The GO-visible internal control now lives only in the
cross-tissue QTL meta, not in this gate.) The GO-invisible modules carry **real
isoform switches** (nearly all members have anticorrelated transcript pairs; max
switch strength 1.10–1.39; 93–459 sig switch tx) that are functionally consequential
(CDS/coding-status/biotype/UTR) **at or above the pooled background** (CDS 0.84 /
coding-status 0.67 / biotype 0.74 / UTR 0.61) — structurally indistinguishable.
Drivers are plausible + psychiatric-relevant (DDX3X, STXBP5, SPTBN1, ARNT2, PBX1,
DGKH, SLC25A12…) but heterogeneous within a module (shared switch axis, not a shared
GO process). **Verdict: PASS in the complementary form** — a genuine DTU-without-DGE
layer invisible to pathway enrichment because the signal is isoform regulation, not a
shared GO term. Frame on mechanism; do NOT claim "pathways WGCNA misses". Manubot
summary: `real_data/brainseq/caudate_sczd/_m/GO_INVISIBLE_GATE_SUMMARY.md`. See
`memory/project_go_invisible_gate.md`.

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
  analysis, removing the constraint baseline). Numbers below are the 2026-06-29
  post-covariate-decouple re-fit inputs (`qtl_anchoring_meta_contrast.parquet`;
  regenerates bit-identically, max|diff| = 0): all 1.068 (p=1.3e-5, I²=0.29),
  pheno-sig 1.163 (p=3.6e-7, I²=0.47), **GO-invisible 1.172 (p=2.3e-5, I²=0.00,
  Q=8.5 — FE and RE identical)**, GO-visible 1.104 (p=0.022, I²=0.68, Q=31.4).
  Splicing-QTL is spared ~7–17% vs expression-QTL in co-switch genes, concentrating in
  the phenotype-associated + GO-invisible modules — where IsoGraph's unique value is
  (DTU-without-DGE).
- **GO-visible is NO LONGER a clean internal null** — do not describe it as one. On the
  re-fit inputs it is 1.104 at p=0.022. What separates it from GO-invisible is now
  *consistency*, not presence/absence: GO-invisible is homogeneous across all 10 tissues
  (I²=0.00) while GO-visible rides on heterogeneity (I²=0.68), so its nominal
  significance is carried by a few tissues rather than a consistent effect. Report it as
  the low end of a gradient.
- **The primary internal control is the matched WGCNA baselines, not GO-visible.** They
  hold the switch features fixed and vary only the inference, and they are null across
  every module set (0.970–1.017, p=0.44–0.86) — a stronger and more interpretable
  control than a module-content contrast.
- Manuscript line: frame on the **sQTL-vs-eQTL specificity contrast**, not raw sQTL
  enrichment. Scope caveat (in report): cis-sQTL anchors member-gene splicing to
  genetics, not the co-switching coordination itself. Optional secondary:
  switch-transcript→LeafCutter-intron direction concordance.

**Matched-baseline control DONE** (`qtl_anchoring.py --method {wgcna_switch_only,
wgcna_multiplex}`; SLURM `real_data/gtex/_h/08.qtl_anchoring_matched.sh`, 2 methods × 13
GTEx tissues; `qtl_anchoring_meta.py` now multi-method). Anchors the matched WGCNA
baselines (same switch features) and compares the splicing-specificity contrast on the
**same 13 GTEx tissues**. Result is a clean **method effect**: only IsoGraph shows it —
go_invisible **1.164 (p=2.0e-4)**, pheno_sig 1.146 (p=1.9e-5), go_visible 1.087 (p=0.060);
`wgcna_switch_only` and `wgcna_multiplex` are **null everywhere** (go_invisible 1.021
p=0.71 / 0.993 p=0.86; go_visible 0.959 p=0.30 / 0.993 p=0.77). This null — not
go_visible — is the primary internal control. Same features + classical inference loses
the splicing-genetic signal
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
residualized features. **COMMITTED** (IsoGraph core `c958174`).

**Covariate policy DECIDED (2026-06-28): residualization is a DISCOVERY knob only.**
Covariates do two jobs: protect module *discovery* (a confound axis fabricates artifact
modules — no downstream test can un-merge them; synthetic ablation proves it), and
adjust trait *inference* (best done jointly with the trait — spline for age, conditioning
on the abundance channel). These are split, not stacked. Previously `feature_scores.parquet`
held the *residualized* matrix (vae.py:660 used the in-place-residualized `switch_matrix`;
the line-665 comment wrongly said "raw"), so `incremental_association` adjusted the SAME
covariates a second time = double residualization. **Fix — DONE, COMMITTED** (IsoGraph core
`c958174`, "Decouple residualization; persist raw feature_scores; add QC + missing-covariate
guard"): vae.py keeps `switch_matrix_raw` and persists feature_scores from it — clustering/VAE
still embed the residualized matrix, but feature_scores is RAW, so the downstream test is the
single inference adjustment. Verified e2e: feature_scores retains the covariate axis (r²≈0.80)
while the embedded matrix has confound_r²→0. `modules.parquet` is UNCHANGED (clustering
input identical, deterministic seed); only feature_scores + downstream (incremental_association,
trait_associations, module_gene_roles, characterize) changed.

- **#31 (silent-skip guard on `build_design_matrix`) — DONE**, same commit.
- **#32 (upstream list technical-confounds-only; launch re-fit) — DONE.** Production re-fit
  regenerated (`8376556` "Regenerate canonical real-data outputs from post-decouple fit";
  res-2.0 comparison arm `fedec9c`). Verified 2026-07-16: all 13 GTEx + 4 BrainSEQ regions
  carry canonical `_m/isograph_vae/modules.parquet`; both repos clean on `main`.

- **#10 (module-level degradation QC flag + edge-type penalties) — WON'T DO for this
  submission** (decided 2026-07-16). Degradation is instead **disclosed as a characterized
  limitation** in the Discussion, citing figS8 (residualization holds for composition/batch/
  depth; degradation is the exception) and figS9 (abundance-channel fallback restores
  recovery). Do not start it without explicit confirmation — it touches the feature
  residualization path.

### 4. Per-module trust funnel writeup — DONE (delivered in manuscript)

Delivered: Results §"IsoGraph modules are reproducible and replicate aging associations"
(266 full-data modules, 236/89% trusted vs size-matched permutation null; driver ρ 0.77–0.82;
25 cross-cohort aging replications vs 6 for the matched abundance baseline), figure
`figTrustFunnel`, and supplementary `tableS7_module_trust_funnel.csv`. Original design notes
below retained for rationale.

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

### 5. Manuscript summary tables + figures — DONE (assembled 2026-07-16)

- Supplementary tables S1–S7 built by `real_data/_h/assemble_supp_tables.py`, with legends
  in `real_data/_m/supp_tables/SUPPLEMENTARY_TABLES.md`; mirrored into the manuscript at
  `content/supplementary_tables/`.
- Supplement wired in the manuscript repo as `content/06.supplement.md`: 14 supplementary
  figures — S1–S11 synthetic (`benchmark/03_metrics/figures/figS1..figS11`) and S12–S14
  real-data (`figBaselineRates`, `figGoInvisible`, `figGwasResolution`) — numbered in
  citation order via `pandoc-fignos` `tag=`, plus S1–S7 table legends. All 17 figure
  cross-refs resolve; no orphaned assets.
- **Remaining:** Cell Genomics re-target (abstract re-lead, STAR Methods + Key Resources
  Table), DOI placeholders (Zenodo ×2, benchmark repo, protocols.io, GTEx access date).

### 6. sQTL intron-direction concordance — DONE (informative null; coloc carries direction)

Approved 2026-07-16 as the single added orthogonal validation; shipped as committed CLIs
`isograph_benchmark/real_data/sqtl_concordance.py` + `sqtl_concordance_meta.py` (commit
01f3c47) with SLURM wrapper `real_data/brainseq/_h/14.sqtl_concordance.sh` (seed 13, 2000
permutation draws), parquet + `real_data/_m/sqtl_concordance_meta/SQTL_CONCORDANCE_META.md`.
**Outcome — a diagnosed null:** the within-gene rank-concordance between the switch axis and
the lead sQTL's per-transcript intron direction is at or below a full-variance permutation null
in every module set (pooled mean |rho| 0.29–0.34, p≈0.76–1.0). Cause is construction, not
biology: a single lead sQTL tags introns with a near-constant net per-transcript sign in ~2/3
of genes, so the within-gene correlation collapses to tie-breaking noise. Directionality is
therefore carried by colocalization (`sqtl_coloc`, GWAS-anchored allele direction), not this
allele-reference-free relative test; the positive anchoring evidence remains the sQTL/eQTL
specificity enrichment (`qtl_anchoring`) + coloc. Reported transparently as an underpowered
negative; does not feed a headline table but is documented for completeness.

### 7. Neuronal CLIP orthogonal validation — DONE (informative sparse null)

Stages 22–24 build switch-localized windows, quarantine the structurally confounded
NOVA1 nomination, and formally classify PTBP2/TDP-43 human contexts as descriptive or
not estimable when discordant support is sparse. Stage 25 is the independent NOVA-family
rescue: `nova_family_renomination.py` plus
`real_data/brainseq/_h/25.nova_family_renomination.sh` scans the full 17-analysis
universe with two-sided intronic opportunity control and gene-level adjusted enrichment.
The frozen result is 79,763 transcripts, 374,296 eligible pairs, and two GO-visible
BrainSEQ hippocampus regulons (M001/M002; 339 genes, 1,993 unique transcript pairs).
Label these candidates `NOVA_FAMILY`, not NOVA1/2. Stage 26 is implemented as
`nova2_ctag_clip.py` plus `real_data/brainseq/_h/26.nova2_ctag_clip.sh`: it acquires
GSE103315 and reciprocal hg38/mm10 chain resources with pinned hashes, reconstructs
strict within-gene matched windows for the frozen candidates, requires full single-block
forward mapping plus >=95% reciprocal overlap, and tests exact NOVA2 coverage separately
in Emx1 cortical, Gad2 cortical, and Pcp2 Purkinje contexts. GEO exposes one pooled
unique-tag bedGraph per context (three biological replicates pooled), not replicate-level
processed peaks or input. Therefore stage-26 effects are explicitly descriptive Tier-2
pooled-signal evidence; absence of coverage is not proof of no binding. Final result QC is
complete (SLURM 42902680, exit 0): only 1,272/42,718 matched window rows (3.0%)
passed strict reciprocal mapping. At the prespecified 100-nt width, only 9 genes had
jointly callable case/control windows and none carried pooled NOVA2 coverage in any
context; the Emx1-vs-Gad2 interaction likewise had 0 informative differences. The
50-nt Emx1 sensitivity was sparse and null (76 genes, 4 case-bound vs 2 control-bound,
OR 1.8, exact p=0.6875; matched-window OR 1.0, p=1.0). **Verdict: informative null /
not testable at adequate power**, not a NOVA2 rescue. Raw replicate reprocessing is
low priority: the nine selected SRA-lite runs total only ~474 MB, but reprocessing
cannot repair the dominant cross-species mapping bottleneck and would require a newly
pinned SRA/alignment environment. GSE103314 acquisition is complete: the six Nova2-cKO
Quantas/BED12 supplements are checksum-pinned at sample/contrast level, and the Emx1,
Gad2, and Pcp2 catalogs pass unique-event-to-coordinate QC.

Stage 27 is implemented as `nova2_perturbation.py` plus
`real_data/brainseq/_h/27.nova2_perturbation.sh`, following the frozen
`real_data/NOVA2_PERTURBATION_ANALYSIS_SPEC.md`. It deterministically reconciles
duplicate Quantas rows, retains one-to-many BED12 event coordinates, derives
strand-aware 50/100/250-nt intronic flanks, requires strict reciprocal mm10/hg38
mapping, and tests event localization in the frozen matched case/control windows with
genes as the primary unit. Counts are pinned by context and width in config. Final QC
is complete (SLURM 42928722, exit 0). At 100 nt, reciprocal mapping retained 586/11,353
Emx1 flanks (5.2%), 435/8,199 Gad2 (5.3%), and 181/2,611 Pcp2 (6.9%), yielding only
2, 3, and 0 localized human windows. Emx1 had 1 case-only versus 1 control-only gene
(327 callable genes; OR 1.0, exact p=1.0); Gad2 had 0 versus 1 (OR 0.33, p=1.0); Pcp2
had no discordance. The 50-nt sensitivity was directionally favorable but below the
prespecified 10-discordant-gene threshold (Emx1 7 vs 2, OR 3.0, p=0.180; Gad2 4 vs 2,
OR 1.8, p=0.688). No gene replicated case-specific localization across both cortical
contexts. One primary-width Emx1 case window at **SPPL2A** had both Nova2-cKO event
localization and independent pooled NOVA2 cTag coverage; this is an exact descriptive
example, not enrichment evidence. **Verdict: informative sparse null / not formally
estimable**, so this does not rescue NOVA2 and candidates remain `NOVA_FAMILY`. Close
the neuronal CLIP analysis here for this submission; integrate the null transparently
in Methods/Results or supplement, and keep controlled-access EGA reanalysis optional.

### 8. RESOLVED 2026-08-26 (option 1: re-freeze to v6) — neuronal-CLIP frozen inputs invalidated by the intronic re-scan

`configs/neuronal_clip.yaml` (version 5, frozen 2026-08-01) declares its
`candidate_source` as `rbp_regulon_intronic.parquet`, `rbp_switch_calls_intronic.parquet`
and `rbp_counts_intronic.parquet`. Those three files were regenerated on 2026-08-26
(SLURM 44484240) when the intronic scan was moved off the legacy flat 0.25 background onto
the GC-binned composition background, which the mature scope had already used since
2026-08-12. The pre-re-scan files are preserved as `*_intronic_flatbg.parquet`.

The config pins `expected_module_rbp_nominations: 17` and the new tables still yield
exactly 17 — **but they are not the same 17**. Over TARDBP/NOVA1/NOVA2/RBFOX2/PTBP2:
8 nominations shared, 9 dropped, 9 added (TARDBP 5 -> 3, NOVA2 3 -> 5); e.g.
`frontal_cortex_ba9/M008/TARDBP` is gone and `frontal_cortex_ba9/M010/NOVA2` is new.
The count guard therefore PASSES while 53% of the candidate set has silently changed —
the failure mode a pinned expectation is supposed to prevent.

**Severity: provenance, not conclusion.** Section 7's verdict is an informative sparse
null driven by the cross-species mm10/hg38 reciprocal-mapping bottleneck (~3-5% of windows
retained), not by which candidates were nominated, so the scientific finding does not move.
What is broken is the claim that the frozen suite is reproducible from its declared inputs.

**Decision taken with the user 2026-08-26: option 1, re-freeze to `version: 6`.** The
composition background is the corrected scan, and pinning the suite to the superseded
flat-background tables would have enshrined the uncorrected background inside a validation
suite with no good answer to "why does your CLIP validation use a different background than
your motif scan?". The re-run is compute-only — `inputs/raw/neuronal_clip` (1.7 GB) and
`reports/neuronal_clip/dataset_manifest.tsv` were already frozen on disk — so
`accessed_date` stays 2026-08-01 and only the candidate-derivation provenance moves.

What landed:
- `configs/neuronal_clip.yaml` -> `version: 6`, `frozen_date: 2026-08-26`, `accessed_date`
  unchanged. Re-pinned `window_stage.expected_candidate_rows` 12779 -> **11261** and
  `expected_unique_candidates` 12038 -> **10410**.
- **The count guard is replaced by a candidate-identity hash.**
  `candidate_source.expected_candidate_identity_sha256` pins the sorted
  (region, module_id, rbp, gene, transcript_id_1, transcript_id_2) digest, reusing the
  `_candidate_identity_sha256` / `_guard_candidate_identity` convention already used by the
  NOVA-family stages. Unpinned -> prints the observed hash and continues (bootstrap);
  pinned and mismatched -> raises. The count guard is kept as a cheap first check.
- 5 tests in `tests/test_neuronal_clip_fetch.py`, including the exact failure mode the count
  guard could not see: same number of nominations, different nominations.
- New committed wrapper `real_data/brainseq/_h/21b.neuronal_clip_freeze.sh` — the freeze had
  been run by hand, which was its own reproducibility gap.

Grounded diff of the re-freeze (old vs new `candidate_manifest.parquet`): 17 -> 17
nominations, **8 shared / 9 dropped / 9 added**; TARDBP 5 -> 3, NOVA2 3 -> 5;
`frontal_cortex_ba9/M008/TARDBP` gone, `frontal_cortex_ba9/M010/NOVA2` new. Candidate rows
12,779 -> 11,261 (9,287 shared). This confirms the count guard was passing on a set that had
changed by 53%.

Stages 21b/22/23 re-run and COMPLETED 2026-08-26 (jobs 44534466/44534467/44534573); the
re-run freeze passed the newly pinned identity hash silently, which is the guard working.
Stages 24-27 chained (44534675-8). Downstream pins in the NOVA-family / ctag / perturbation
stages (`expected_candidate_rows: 1993`, `expected_candidate_identity_sha256: 8970b1cf...`,
`expected_context_counts`, `consensus_sha256`) are derived from the OLD candidate set and
are expected to need re-deriving as those stages run.

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
