# Analysis TODO — Cell Genomics

Reconciled against the repository 2026-08-26. Every `[x]` names the committed CLI
and wrapper that closed it, so the claim can be checked rather than trusted.
Status vocabulary: `[x]` done, `[~]` partially done (the gap is stated), `[ ]` open.

Superseded framing note: the old "Additional analysis for Nature Communications"
section is gone — the target is Cell Genomics (retargeted 2026-07-16) and all three
of its items are now covered under Mechanistic validation below.

## Done — do not redo

* [x] **Validate major IsoGraph switches using junction counts or PSI.**
  `validate_switch_splicing.py` + `real_data/brainseq/_h/16.validate_switch_splicing.sh`.
  Caudate marginal OR ~105 vs module OR ~2.3.
* [x] **Adjust human-cohort analyses for estimated cell-type composition.**
  `celltype_composition.py` — joins the committed MuSiC BrainSEQ deconvolution and
  supplies a marker-depletion cut; consumed by `incremental_association`.
* [x] **Reconcile and validate the benchmark recovery-metric definition.**
  `benchmark/partition_metrics.py` + `backfill_metrics.py`, 15 unit tests,
  `benchmark/03_metrics/_m/PARTITION_METRICS.md`. Best-match Jaccard reproduces the
  published `module_recovery_score` to 2.22e-16 over all 13,060 runs, and the
  count-preserving Jaccard null shows WGCNA's negative-control recovery (0.500) is
  fully explained by its module-count structure (z=0.47, perm p=0.746).
* [x] **Direct genetic analysis of module eigengenes.**
  `module_genetic_anchoring.py` + `real_data/brainseq/_h/21.module_anchoring.sh` —
  per-module splicing-specificity contrast against a permutation null of random
  gene sets, closing the "pooling does not show any individual module is anchored"
  objection. Plus S-LDSC + coloc (`ldsc_annot_prep.py`, `coloc_*.py`).
* [x] **Ablations: VAE denoising, estimability weighting, multiplex integration,
  abundance-channel inclusion.** Carried as first-class methods in the synthetic
  grid — `isograph_baseline`, `isograph_spearman_leiden` (Leiden without the VAE),
  `isograph_vae`, `isograph_vae_reliability` (estimability), `isograph_vae_multiplex`,
  `isograph_vae_residual` — plus `project_tiers.py`, which derives the switch_only /
  switch_primary / multiplex edge tiers from ONE fit so the latent is held constant.
  Residualization cost on unconfounded data: `benchmark/03_metrics/_m/RESIDUAL_COST.md`.
* [x] **Sensitivity to counts versus TPM.** N/A — no method ever consumes TPM.
  BrainSEQ uses Salmon counts, GTEx RSEM expected counts, synthetic simulated counts;
  TPM is archived as a data product (`build_parquet.convert_gtex_transcript_tpm`) and
  never fed to a model (`build_bundles.py:230,249`). Manuscript methods corrected.

## Verification audit 2026-08-28 — are the open analysis tasks still needed?

Checked each of the 7 remaining analysis items against the code and the data on disk, on
the suspicion that TODO.md had gone stale. **Six are accurate as written. One is not.**

| # | task | verdict |
|---|---|---|
| 1 | long-read matched null | still open — `longread_switch_confirm.py` has no null of any kind |
| 2 | BrainSEQ<->GTEx eigengene projection | still open — `replication.py` is overlap-based, `scz_age_projection.py` projects *within* BrainSEQ (aging modules -> SCZD cohort), neither is cross-cohort eigengene projection |
| 3 | junction/PSI network + transcript-level WGCNA | still open — no match anywhere in the tree |
| 4 | phenotype-blind resolution/edge selection | still open, but **smaller than stated** — see below |
| 5 | SCZ medication/toxicology/smoking confounders | **premise is wrong — the data does not exist** — see below |
| 6 | five-part sensitivity block | still open — the knobs exist as hard-coded defaults, no harness varies them |
| 7 | SNCA orthogonal confirmation | still open, and **the machinery to do it already holds the answer** — see below |

**#4 is smaller than TODO.md implies.** `stability.py` already implements exactly the
phenotype-blind criterion the reviewer asked for — split-half refit, ARI/NMI over genes
assigned in both halves, seeds 0..k-1. It simply has no `--leiden-resolution` argument, so
it has only ever run at the canonical resolution. What is missing is a sweep, not a method.
Separately, `leiden_sweep_results.parquet` does carry two phenotype-blind columns
(`giant_fraction`, `nmi_to_prev`) alongside the phenotype-aware `n_sig_linear_fdr10` /
`n_sig_spline_fdr10` — but `nmi_to_prev` compares adjacent resolutions within one fit, which
is a smoothness diagnostic, not a resampling stability.

**#5's premise is wrong: BrainSEQ does not release these variables.** The BrainSEQ colData
carries 130 columns and every one beyond `Dx / Age / Sex / Race / PMI / MoD / RIN / Protocol
/ SNP_PC1-10` is sequencing QC (`rse-gene.bsp3.caudate-n487.gencode-v47.RData`, checked
2026-08-28). There is no medication, antipsychotic, toxicology, nicotine or smoking field —
not in the bundles (`inputs/bundles/brainseq_*/*/samples.parquet`, 28 columns), not in the
source RSEs, and not anywhere under `shared/resources/libd-data`. `clinical_consequence.py`
is unrelated despite the name: it is ClinVar/gnomAD variant constraint. The reviewer's
request was explicitly "**where available**", and the honest answer is that they are not.
The task is therefore rescoped, not dropped: audit availability as a first-class output,
run the sensitivity over every confounder that *is* available, and use molecular proxies
where a defensible one exists — flagged as proxies.

**#7's answer is already computable from existing outputs, and it is a negative.** The
SNCA sQTL that colocalizes in LBD and PD tags junction `chr4:89835692-89836127(-)` carried
by `ENST00000508895.5`. That transcript *is* detected in the ONT DLPFC long-read data
(`longread_switch_confirm/pair_confirmation.parquet`), but at a mean isoform fraction of
0.0029, its switch partner `ENST00000394989` at 0.0011, and their usage correlation is
**+0.635** — positively correlated, `switch_like = False`. The gene passes gene-level
confirmation in both regions while the specific anchored pair does not behave as a switch
in an orthogonal technology. This needs to be a committed CLI over all 12 splicing-led
genes, not an observation about one gene, and it bears directly on whether SNCA can carry
a main figure.

## Partially done — the gap is the work

* [~] **Confirm top transcript pairs and disease-linked events using long-read brain
  RNA-seq.** Tier-1 done (`longread_switch_confirm.py` +
  `real_data/brainseq/_h/29.longread_switch_confirm.sh`): ONT DLPFC BA9/46, Zenodo
  8180677 Bambu quants, confirms 60.5% of GTEx cortical switch genes.
  A matched null now exists (`switch_orthogonal_confirm.py --mode global-null`, 2026-08-28)
  and it is sobering: switch pairs are switch-like 0.6466 of the time against an
  abundance-matched null of 0.6391 from the *same genes* (p=0.022, difference +0.0075 —
  significant only because n≈18,000). Non-switch pairs sit at a mean usage correlation of
  -0.172 before any biology, because within-gene fractions sum to one. **A bare negative
  usage correlation is therefore close to vacuous as switch evidence**; read the rate
  against this null, never against zero, and prefer within-switch-universe contrasts.
  **GAP: this uses the detected-pair denominator while the published 8.4% uses all
  prespecified pairs, so the 8.4% figure itself still needs its own null before it goes
  in the paper.**
* [~] **Frozen-module projection between cohorts.** `replication.py` does bidirectional
  module preservation (best-Jaccard + size-preserving label-permutation null) and
  `scz_age_projection.py` freezes out-of-cohort aging modules and projects them onto
  the independent SCZD cohort. **GAP: no BrainSEQ->GTEx / GTEx->BrainSEQ eigengene
  projection specifically** — preservation-by-overlap is not the same test.
* [~] **Workflow-level comparators.** Done: conventional abundance WGCNA
  (`wgcna_gene`), matched-feature WGCNA (`wgcna_switch_only`, `wgcna_multiplex`,
  `run_matched_wgcna.py`), CLR switch coordinates + Leiden without the VAE
  (`isograph_spearman_leiden`). **GAP: no junction/PSI network comparator, and no
  transcript-level (as opposed to gene-level-on-switch-features) WGCNA.**
* [~] **Select Leiden resolution and edge thresholds using phenotype-blind stability
  criteria.** Resolution 5.0 is canonical and `sweep_leiden.py` +
  `05.sweep_leiden.sh` / `07.sweep_leiden_with_abundance.sh` compute the sweep; the
  giant-module size criterion (>=900 genes) is itself phenotype-blind.
  **GAP: the published justification is phenotype-AWARE** —
  `real_data/gwas/_m/GWAS_RESOLUTION_SUMMARY.md` argues 5.0 by showing it removes a
  GWAS enrichment artifact. A reviewer asked for selection on stability BEFORE any
  trait is consulted. Needed: a phenotype-blind stability-vs-resolution curve with
  the choice made on it, and the GWAS result demoted to post-hoc confirmation.
  Edge thresholds have no documented selection criterion at all.

## Open

* [ ] **Known-junction orthogonal validation of the SNCA and CTSH coloc events
  (BrainSEQ short-read, in-house).** The long-read arm failed on these two genes for two
  *different* reasons, and neither is "the switch is not real":
  **CTSH is a region mismatch** — its coloc tissue is hippocampus (CLPP 0.386) but the ONT
  data is DLPFC BA9/46, so CTSH was never tested in its own region.
  **SNCA is a platform mismatch** — its coloc tissues (cortex, frontal_cortex_BA9) DO match
  the long-read tissue, so region is not the explanation; the assay is ONT **cDNA**
  (SQK-PCS111), which truncates at the 5' end, and SNCA's event is an alternative **first
  exon** — exactly the event class that protocol under-detects. Its 0.29% usage is
  therefore not interpretable as absence.
  **The decisive test is already on disk and is far better powered.** Both colocalizing
  junctions are present as *known* junctions in the BrainSEQ junction annotation —
  SNCA `chr4:89835692-89836127(-)` and CTSH `chr15:78937423-78937686(-)` /
  `chr15:78937423-78939140(-)` — and BrainSEQ has hippocampus n=452, DLPFC n=500,
  caudate n=487, i.e. ~40x the ONT n=12, measuring the junction the sQTL actually tags
  rather than a whole-transcript proxy.
  **Work:** extend `validate_switch_splicing.py` (already built for junction/PSI
  validation) with a targeted per-junction mode — take the 12 splicing-led coloc events,
  resolve each to its annotated junction, compute per-sample PSI against the gene's
  competing junctions, and test (i) that the junction is measurably used in the coloc
  tissue and (ii) that its PSI anti-correlates with the partner junction beyond the
  compositional-closure baseline established by `switch_orthogonal_confirm --mode
  global-null` (never against zero). Run CTSH in **hippocampus** and SNCA in
  **DLPFC + caudate**. Wrapper `real_data/brainseq/_h/38.junction_coloc_confirm.sh`.
  **Decision rule, pre-registered:** if the junctions validate here, report the
  short-read junction result as the orthogonal confirmation and cite the long-read
  failure as an assay limitation (5' bias / region), not as a negative. If they do not
  validate, SNCA and CTSH stay off any main figure and the set-level result in
  `switch_orthogonal_confirm --mode anchored` (0.453 vs 0.252, p=5e-4) stands alone.
  This supersedes commissioning new long-read data as the next step; see the
  platform/region note there if long-read is revisited (hippocampus for CTSH;
  5'-complete chemistry — direct-RNA, 5'-capture, or PacBio Iso-Seq — for SNCA, since
  more ONT cDNA would reproduce the same bias; disease-specific tissue is NOT required,
  because the question is whether the anchored isoform reaches measurable usage and
  anti-correlates with its partner, which is a normal-variation question — the coloc
  already supplies the disease link).

* [x] **Test schizophrenia findings for medication, toxicology, smoking and related
  confounding where available.** DONE 2026-08-28 — `scz_confound_sensitivity.py` +
  `real_data/_h/36.scz_confound_sensitivity.sh`. The reviewer's "where available" is the
  operative clause: BrainSEQ releases **none** of medication, toxicology or smoking, so
  the CLI emits the availability audit as a re-runnable output and then tests the three
  tiers it can (measured covariates the published model omits; molecular proxies for
  smoking and antipsychotic exposure built from raw bundle counts; and what stays
  untestable). Baseline verified against `diagnosis_assoc.parquet` (max |diff| 5.7e-15).
  **The rRNA-rate result needs care, and an earlier reading of it here was wrong.**
  rRNA rate alone takes the FDR<0.05 modules from 7 to 4 and the fully adjusted model to 2,
  but that is *not* a confounding result and must not be reported as one. `r_rna_rate` has
  **no marginal association with diagnosis** (Cohen's d = -0.009, p = 0.93), so it cannot
  confound the diagnosis-outcome relationship; its association appears only after
  conditioning on the published covariates (partial r = +0.108), which is the collider-bias
  signature, not the confounder signature. Ruled out as alternative explanations: it is not
  double-correction against the VAE (`residualize_composition=False`, so `feature_scores`
  are RAW -- verified by the exact rebuild in `switch_feature_sensitivity.py`, max |diff| 0),
  it is not double-counting the mito axis (partial r with `mito_rate` given the published
  model = -0.0000), and it is not power loss (SE ratios 0.963-1.006; the effects genuinely
  shrink 4-23%). There is also no dose-response: the modules *most* correlated with rRNA
  rate survive (M010 -0.344, M020 -0.333) while the ones that fall are the *least*
  correlated (M022 -0.121, M023 -0.112) -- the opposite of what technical-artifact removal
  would produce. **Do not add `r_rna_rate` to the inference model on confounding grounds.**
  The policy-consistent action, if any, is at *discovery*: `r_rna_rate` is the one major
  technical axis absent from `BRAINSEQ_DISCOVERY_COVARIATES` (which carries RIN,
  mapping_rate, mito_rate, SNP_PC1-5). Report the sensitivity as a sensitivity.
* [x] **Remaining sensitivity analyses** — harness built and RUN 2026-08-28,
  `switch_feature_sensitivity.py` + `real_data/_h/37.switch_feature_sensitivity.sh`,
  covering all five axes (pseudocount, transcript-expression filter, minor-isoform
  threshold, identifiability by transcript number, quantification pipeline). Gated on an
  exact rebuild of the published switch channel (max |diff| 0.0). **GAP: the module
  partition is held fixed**, so this measures the stability of the representation and its
  trait signal, not of an independently refit network; a full refit per setting is the
  expensive follow-on. The quantification axis is a cross-cohort concordance and therefore
  confounds quantifier with cohort — an upper bound, not an isolated estimate.

  **Results (brainseq/caudate, gate passed at max |diff| = 0):** pseudocount is a
  non-issue (median per-gene |r| >= 0.992 over 0.1-2.0, 20/20 published-significant
  module-age associations retained). The expression filter and minor-isoform threshold
  matter more — loosening to count>5/frac>=0.5 gives median |r| 0.842 and 18/20;
  min_usage=0.10 drops 2,727 genes and retains 15/20. **No sign flips among
  published-significant modules under any setting.** Identifiability is a mild monotone
  trend (median |age r| 0.066 -> 0.102 across transcript-number strata, module membership
  0.28 -> 0.62) but the |r|>0.2 tail is flat, so it reads as a membership effect, and the
  top stratum has 69 genes. **Quantification is the striking one:** per-gene switch-age
  effects are essentially uncorrelated between Salmon and RSEM on matched regions
  (Pearson 0.007 caudate, -0.002 hippocampus, sign concordance 0.499).

## Mechanistic validation (was "Additional analysis for Nature Communications")

At least one decisive mechanistic validation was required; the RBP arm has been
carried furthest.

* [x] **Demonstrate direct genetic regulation of a co-switching module.**
  `module_genetic_anchoring.py` (see above).
* [~] **Validate a recurrent predicted RBP regulon using binding or perturbation
  data.** Binding: `rbp_binding.py` (eCLIP, length-normalised density + Haldane OR),
  `rbp_regulon.py` (covariate-adjusted binomial GLM with an estimability gate),
  neuronal-CLIP suite (`neuronal_clip_*.py`, stages 21b-27), `nova2_ctag_clip.py`,
  `nova2_perturbation.py`. Verdict is an **informative sparse null** driven by the
  cross-species mm10/hg38 reciprocal-mapping bottleneck (~3-5% of windows retained).
  Perturbation is designed but not performed:
  `real_data/_m/rbp_target_panel/WETLAB_PERTURBATION_DESIGN.md` with the assayable
  pair list from `rbp_pair_assayability.py` + `real_data/_h/34.rbp_pair_assayability.sh`
  (NONO 11 / ELAVL1 6 / KHDRBS1 8 measurable-and-responsive pairs).
  **Wet-lab CLOSED 2026-08-28 (out of scope, prioritized and handed off).**
  `real_data/_m/rbp_target_panel/WETLAB_HANDOFF.md` states the priority order the design
  doc did not: P1 = NONO + KHDRBS1 together (the adjudicating contrast; NONO alone can
  only confirm and cannot separate a regulon from a generic abundant-RBP effect), P2 =
  add ELAVL1 (broadest nomination but the thinnest assayable panel, 2/14 tier-1 genes
  responsive), P3 = secondary endpoints + same-model CLIP before any direct-regulation
  claim. The limitation is now written into the manuscript discussion
  (manuscript repo `content/04.discussion.md`, commit ba19cdd) together with the
  HepG2/K562-not-brain eCLIP caveat, rather than left implied.
  The HepG2/K562-not-brain caveat on the ENCODE eCLIP panels is written into
  `real_data/_m/rbp/RBP_REGULON_SUMMARY.md` (Results and Limitations) and into
  `rbp_binding.py`, so it survives regeneration.
* [x] **Validate the shared SNCA alternative-first-exon mechanism across independent
  data types or cohorts.** DONE 2026-08-28 (computational arm) —
  `switch_orthogonal_confirm.py` + `real_data/_h/35.switch_orthogonal_confirm.sh` scores
  the sQTL-anchored transcript pair itself, for all 12 splicing-led genes, in ONT
  long-read DLPFC against switch pairs matched on abundance decile.
  **Set-level: confirmed.** 0.453 switch-like vs a matched null of 0.252 (p=5e-4);
  0.600 vs 0.304 restricted to pairs whose anchored isoform is usably expressed.
  **SNCA specifically: not confirmed.** Its sQTL-carrying transcript `ENST00000508895`
  sits at 0.29% of SNCA's long-read output, 4 of its 5 anchored pairs are *positively*
  correlated, and it fails the abundance qualification entirely — as does CTSH, the only
  high-confidence coloc. Report the splicing-led set, and do not promote SNCA or CTSH to
  a main figure on orthogonal grounds. This reinforces Table 3's existing set-level
  caption rather than contradicting it.
  **Follow-on:** the two failures are region- and platform-explained, and the better-powered
  test is in-house short-read junctions -- see the open "Known-junction orthogonal
  validation of the SNCA and CTSH coloc events" item under **Open**.

## Manuscript mechanics — MOVED OUT 2026-08-28

Abstract re-lead, STAR Methods conversion, the Key Resources Table and the 5 identifier
placeholders now live in the manuscript repository
(`/ocean/projects/bio260021p/kbenjamin/manuscript/isograph-brain-manuscript`, `TODO.md`,
section "Cell Genomics re-target", commit 069c8e6). They are manuscript mechanics, not
analysis. Closed here.
