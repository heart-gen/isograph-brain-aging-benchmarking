# 05 — Genetic anchoring

**Question:** are the co-switch modules genetically real, and does disease risk resolve
onto specific isoform switches? This stage carried the paper's headline until 2026-09-12; it now
supplies **supporting, set-level** genetic evidence — the headline is the reproducible,
abundance-independent switch layer (`manuscript/MANUSCRIPT_PLAN.md` §10).

## Order

Run the whole stage with `bash 05_genetic_anchoring/_h/run_stage.sh` (add `--dry-run` to print
the plan). The leading number of a wrapper is its tier; steps in one tier run in parallel, except
the few `run_stage.sh` serializes because they write one shared file.

| Step | Wrapper | Waits on | Produces |
|---|---|---|---|
| 01a–01b | `qtl_anchoring`, `qtl_anchoring_matched` | stages 02–03 | Power-matched logistic sQTL/eQTL enrichment per analysis (17 analyses), IsoGraph and the matched WGCNA baselines |
| 01c | `qtl_anchoring_sensitivity` | stages 02–03 | Pre-specified sensitivity arms: constraint-adjusted (gnomAD LOEUF + missense z + log expression, both covariate sets on the identical constraint-complete subset), threshold-free continuous (rank-INT of -log10 pval_beta), and SuSiE credible-set dose. Writes to its own files; the primary arm is untouched. The matched baselines' arms are 01b with `--outcome` / `--covariate-set` |
| 01d | `sqtl_concordance` | stage 02 | Switch-transcript → LeafCutter-intron direction concordance |
| 01e | `module_anchoring` | stage 03 enrichment | Module-level genetic anchoring |
| 01f | `prep_module_gene_sets` | stage 02 | MAGMA module gene sets (`configs/gwas_magma.yaml`); also with `MAGMA_ISOGRAPH_BACKEND=isograph_vae_res2` |
| 01g | `coloc_prep` | stage 02 | Switch genes under each trait's GWAS peaks + their GTEx v11 brain QTL credible sets; per-locus GWAS z (6 analyses) |
| 01h | `ldsc_annot_prep` | stage 02 | S-LDSC sQTL/eQTL/cis switch-gene annotations (`brainseq-sczd`; the `aging` bundle) |
| 01i | `brainseq_switch_qtl` | stage 02 | BrainSEQ cis-swQTL on `S_g` and matched cis-eQTL on `A_g`, same donors/variants/covariates (GPU). **The first all_samples swQTL results (2026-09-10) are superseded:** `build_phenotypes` omitted the discovery transcript filter, so `S_g` was not the discovery switch coordinate (median \|r\| 0.31-0.35 on shared libraries; 1.000 with the filter). Re-mapped 2026-09-11; `A_g` was unaffected. `ea_only` uses within-EA genotype PCs |
| 02a | `qtl_anchoring_meta` | 01a–01c | Random-effects pooling + the paired splicing-specificity contrast; sensitivity arms under `sensitivity/<outcome>_<covariates>/` |
| 02b | `sqtl_concordance_meta` | 01d | Concordance rollup |
| 02c | `module_anchoring_meta` | 01e | Module-level anchoring rollup (Table S16) |
| 02d | `run_magma` | 01f | MAGMA gene and gene-set analyses, 8 traits |
| 02e | `locus_ld` | 01g | Per-locus 1000G EUR LD matrices + QTL variant → rsID bridge |
| 02f | `ldsc_make_annot_ldscores` | 01h | Per-chromosome annotations and LD scores |
| 02g | `brainseq_qtl_checks` | 01i | Gates before any BrainSEQ QTL result is read: `S_g` must reproduce the discovery coordinate (median \|r\| >= 0.99) and is sign-pinned to it; `A_g` eQTL positive control against tissue-matched GTEx v11 eGenes (pi1 >= 0.5 and direction concordance >= 0.90 at GTEx lead variants, both pre-specified). all_samples control PASS in all three regions (pi1 0.83-0.88, concordance 0.987-0.991) |
| 02h | `brainseq_switch_qtl_meta` | 01i | swQTL-vs-eQTL paired tables and modality contrast, both arms |
| 03a | `plot_magma` | 02d | `magma_results_combined[_res2].parquet` + figures |
| 03b | `coloc_clpp` | 02e | GWAS SuSiE + eCAVIAR CLPP against sQTL/eQTL, per analysis |
| 03c | `coloc_gwas_susie` | 02e | **Signal-level coloc, stage A:** the per-locus GWAS SuSiE cache (primary grid at 12,000 SNPs; the scoped `aging__ad` sensitivity at 30,000) |
| 03d | `ldsc_munge_h2` | 02f | Partitioned heritability per trait on baselineLD v2.2 |
| 03e | `brainseq_effect_size` | 02h | Switch vs abundance QTL effect sizes without selection asymmetry |
| 04a | `coloc_meta` | 03b | Cross-trait CLPP rollup |
| 04b | `coloc_direction` | 03b | Signed risk-allele direction + resolved isoform events (CLPP layer) |
| 04c | `module_coloc_convergence` | 03a, 03b | Module-level coloc convergence for AD/PD/LBD/ALS/SCZ: size-matched permutation null + anchored-module enrichment against the CLPP-tested pool. Generalises the SCZ-only layer in `scz_age_projection.py` |
| 04d | `coloc_modality_prep` | 03b | Per-gene sQTL-vs-eQTL coloc targets + GTEx variant bridge, per gene-pool arm |
| 04e | `ldsc_summary` | 03d | `ldsc_partitioned.parquet` |
| 05a | `coloc_modality_abf` | 04d | `coloc.abf` per brain tissue, per arm |
| 05b | `coloc_signal_susie_prep` | 04d (switch arm) | Signal-level target grid, GTEx credible sets, `work_list.tsv` |
| 05c | `deep_dive_events` | 04b | `deep_dive_events`, the per-event table stage 06 reads; the per-gene panel is `08_integration/_h/01a` |
| 06a | `coloc_modality_meta` | 05a | Paired McNemar / Wilcoxon contrast, per arm |
| 06b | `coloc_signal_susie` | 03c, 05b | **Signal-level coloc, stage B:** each gene's cis QTL re-fit with `susie_rss`, `coloc::coloc.susie` signal against signal, per (analysis, tissue); representative and all-introns sQTL arms. Estimator hierarchy is **susie > abf > CLPP**; cells SuSiE cannot speak to fall back to `coloc.abf` and are flagged, never dropped |
| 06c | `coloc_brainseq_prep` | 03c, 05b, 02g | BrainSEQ signal-level coloc targets |
| 07a | `coloc_modality_compare` | 06a | The four gene-pool arms side by side |
| 07b | `coloc_signal_susie_meta` | 06b | Cell tables and the estimator hierarchy; `--sqtl all` gives the primary nominations, `--max-snps 30000` the sensitivity arm |
| 07c | `coloc_brainseq_susie` | 06c | Array over (analysis, region); see 08c |
| 08a | `coloc_isoform_events_signal` | 07b, 04b | Isoform events for the signal-level nominations (read by stage 06) |
| 08b | `locus_event_audit` | 07b, 06a, 05b | **Does a colocalizing locus name the event it claims to?** GWAS signal → sQTL signal → intron phenotype → driver transcript, matched against `configs/known_splice_events.yaml`; abf, susie, all-introns and SNP-guard audits in turn |
| 08c | `coloc_brainseq_meta` | 07c, 07b | BrainSEQ signal-level coloc report |
| 08d | `locus_ld_robustness` | 07b | Locus LD audit for any locus carrying a biological claim: convergence, credible-set purity, boundary trims, `coloc.susie` across the p12 sweep, and LD-mismatch diagnostics (`kriging_rss` outlier drop, residual-variance refit). **Manual:** one run per locus, arguments chosen from the nominations |
| 09a | `smr_heidi_submit` | 07b (GTEx); 08c, 02g (BrainSEQ) | Runs SMR prep, then submits 10a and 11a with the array sized from `work_list.tsv` |
| 10a | `smr_heidi` | 09a | **SMR + HEIDI** on the signal-level nominations (GTEx v11 QTL, 1000G EUR LD). Orthogonal corroboration beneath coloc: `b_SMR` is not causal direction, a non-rejected HEIDI is not proof of sharing, and a HEIDI rejection does not overrule coloc. BrainSEQ source (`--qtl-source brainseq --arm ea_only`): gated on both `02g` checks; GTEx nominations plus BrainSEQ's own `S_g` nominations in all three regions, `S_g` slopes sign-pinned to the discovery axis, read against BrainSEQ's coloc for the same region and axis. **Two testing families, corrected apart** (primary confirmatory probe per gene; the gene's other introns as a secondary event-localization family), an `smr_status` axis that separates a probe never instrumented from one tested and null, `F = (b/se)²` reported per family and never used to exclude, TSS-fallback probes flagged as QC, and a SNP-attrition trace whose dominant filter is **GWAS locus coverage, not allele QC** — panel rsID overlap and allele agreement are both 100%, and BESD-SNP retention must never be compared across QTL sources because its denominator is imputation density |
| 11a | `smr_heidi_meta` | 10a | SMR assembly per arm and source: primary, `--peqtl-smr 1e-6`, `--smr-multi` (submitted by 09a) |

The per-gene deep dive (`build_deep_dive`) and the SCZ age projection read stages 06 and 07, so on
2026-09-15 they moved to `08_integration/_h/01a` and `01b`; only the per-event table (05c) stays here.

## Re-running over an existing tree

Several steps keep what is already on disk, so a re-run after an upstream change (a re-fit, a new
transcript filter) must clear these first, or it silently mixes generations:

- `_m/coloc/<analysis>/susie/` — 02e skips any `<LOCUS_ID>.unphased.vcor1.bin` present, and 01g
  re-derives the loci from the switch genes, so a reused id can pair an old LD matrix with a new SNP list.
- `_m/coloc_signal_susie/gwas_susie/` and `_m/coloc_signal_susie/sensitivity/max_snps_30000/gwas_susie/` —
  03c writes `<LOCUS_ID>.rds` only for loci with a credible set and never removes old ones; 06b and 07c
  read whatever is there.
- `_m/ldsc/<annotation>/{beds,ldscores,results}` — 02f skips existing annotation and LD-score files.
- `_m/smr_heidi/<source>/{besd,smr,sensitivity}` — 11a collects every task under `smr/`, including tasks
  from a longer earlier work list.

The GWAS-only caches (MAGMA gene analysis, munged sumstats, `coloc/_tmp`, the GTEx variant bridge) do
not depend on the modules and can stay.

## The set-level contrast, stated honestly

> **Not the headline (2026-09-12).** Every per-gene view runs the other way — GTEx signal-level
> 22 splicing-only vs 45 expression-only genes (P = 0.007), BrainSEQ in-sample coloc 5 switch-only
> vs 33 abundance-only (P = 4.3e-6) — and the matched-baseline control narrowed to p = 0.053.
> Quote the contrast as set-level support, always with those results beside it.

Co-switch genes are cis-QTL **depleted for both** sQTL and eQTL — coordinated network
genes are constrained. That shared baseline is not the result. The result is the paired
**splicing-specificity contrast** (sQTL OR / eQTL OR within analysis), which removes it:
all 1.068 (95% CI 1.037–1.100, p=1.3e-5, I²=0.29), **pheno-sig 1.111 (1.048–1.177,
p=3.6e-4, I²=0.23)**, GO-invisible 1.068 (0.993–1.148, p=0.077, n.s.), GO-visible 1.084
(1.000–1.174, p=0.050).

Numbers are the 2026-08-29 refresh (`qtl_anchoring_meta_contrast.parquet`). **Never quote
a pre-2026-08-29 anchoring number** — the earlier figures rested on a stale
`module_enrichment` join that relabelled the GO partition.

Two things not to overclaim:

- **The effect does not localise to the GO-invisible modules.** GO-invisible (1.068,
  p=0.077) and GO-visible (1.084, p=0.050) are indistinguishable, so neither is a null and
  neither is the site of the effect — it is carried by the **phenotype-associated** set.
  The DTU-without-DGE claim rests on the GO-invisible content gate (S-real-3), not on
  genetics. The primary internal control is the **matched WGCNA baselines** from stage 02,
  which reach significance in no module set (0.968–1.066, p=0.053–0.71, fixed effects).
  After the 2026-09-09 baseline re-fit the closest baseline is `wgcna_multiplex` on the
  phenotype-associated set at 1.046 (p=0.053, I²=0.75), so state the control as "no baseline
  clears 0.05" rather than "the baselines are flat"; the pre-re-fit range (0.959–1.021,
  p=0.30–0.86) is stale.
- **Per-gene CLPP posteriors are individually modest** (inclusion at ≥0.01; only 4/12
  reach 0.05; only CTSH at 0.39). The genetic-anchoring significance rests on the
  set-level contrast, not on per-locus colocalization. S-LDSC is not independent support for
  splicing: on `coef_p` the splicing annotation does not survive correction. The 12 splicing-led genes
  are the **eCAVIAR layer**: locus nominations now come from the signal-level hierarchy below
  (42 nominations, of which 6 are among the 12), so never quote the 12 as signal-level.
- Scope: cis-sQTL anchors member-gene *splicing* to genetics, not the co-switching
  coordination itself.

## Signal-level coloc and the event audit — what they changed

Two things a reader should know before quoting a per-locus number.

**Most of the grid cannot be scored at signal level.** `coloc.susie` needs credible sets on
both sides. On the primary grid only **230 of 579** GWAS loci fine-map (286 have no GWAS
credible set, 61 exceed the 12,000-SNP guard, 2 have too few SNPs), and only **2,933 of
42,130** (cell, modality) rows -- 7.0% -- are scored by `coloc.susie`. The rest fall back to
`coloc.abf` and carry `estimator = "abf"` plus a `fallback_reason`. That is a property of the
GWAS and the reference LD, not a bug. **25 of the 42 all-introns nominations are headlined by an
abf cell**: CTSH's locus has no GWAS credible set, and PICALM's was never fit because it exceeds
the SNP guard. Quote the fallback rate whenever the hierarchy is described.

**SNCA survives, and the multi-signal case is real.** At SNCA/LBD, SuSiE resolves **two**
QTL credible sets and only one colocalizes (PP4 0.974 on rs6532192; the other is
PP3 0.999). That is exactly the situation `coloc.abf`'s single-causal-variant assumption
cannot represent. PP4 stays high across the `p12` sweep (0.79 at 1e-6 to 0.997 at 1e-4), but
the 0.8 call itself holds only from p12 = 5e-6 (`prior_robustness = intermediate`); say
"near the call at the most conservative prior", not "holds across the sweep".

**Neither UNC13A nor PICALM is a recovered mechanism on the current nominations.** The
audit curates each locus' literature event to verified GRCh38 coordinates and asks whether
the colocalizing intron IS that event:

- **UNC13A** — the TDP-43 cryptic exon sits in intron 20, `chr19:17,641,557-17,642,844`
  (derived here: rs12608932 at chr19:17,641,880 in the project's own panel, placed by
  GENCODE v47). The colocalizing intron is `chr19:17,630,750-17,632,782`, ~9 kb away.
- **PICALM** — the AD-associated event is the TWAS cluster anchored on rs3851179,
  `chr11:86,026,368-86,031,468` after liftover. The colocalizing intron is
  `chr11:85,974,812-85,981,129`. Note the locus carries a *second*, different PICALM sQTL
  (`clu_7402`, lead rs540422) that does **not** carry the AD signal — picking it would
  have mis-stated the result.

In both cases the curated event **is testable in GTEx** but is not GTEx's
grouped-permutation representative intron, so only the all-introns arm puts it on trial. That
arm settles it (primary audit, 2026-09-11, 42 nominations): **0** `known_mechanism_recovered`;
**UNC13A** `context_distinct_splice_colocalization` (the cryptic-exon intron was tested in the
same two cerebellar tissues and did not colocalize); **PICALM** `disease_locus_splice_linked` on
the primary grid (abf only, locus over the SNP guard) and `context_distinct` only in the
30,000-SNP sensitivity arm; the other **40** `novel_splice_led_candidate`.

That arm writes under `_m/coloc_signal_susie/all_introns/` and is read back with
`--stage meta --sqtl all`. Both modes name their shards `<analysis>__<tissue>.parquet`, so
the directory split is what keeps an all-introns run from overwriting the representative
results it exists to be compared against.

A raised GWAS SuSiE SNP guard is split out the same way. `MAX_SNPS = 12,000` is the uniform
primary grid; a run at any other `COLOC_GWAS_MAX_SNPS` writes under
`_m/coloc_signal_susie/sensitivity/max_snps_<N>/`, is read back with
`--stage meta --max-snps N`, and is audited with `locus_event_audit --max-snps N`. The only
such arm is the scoped `aging__ad` recovery at 30,000, run to fit PICALM's locus. Loci
recovered there were empirically enriched for GWAS-reference-LD inconsistency, so a
recovered locus enters the narrative only after `08d.locus_ld_robustness.sh`. For PICALM
that audit is LD-robust — one 2-variant credible set (rs10792832/rs3851179, purity 0.994),
unchanged under 100–500 kb boundary trims, kriging outlier removal and
`estimate_residual_variance` — but **prior-sensitive**: cortex PP4 0.813 at p12 = 1e-5,
0.685 at 5e-6, 0.303 at 1e-6.

Every hierarchy cell also carries two **evidence-strength descriptors**, and the event
audit reports them per nomination at the headline tissue. They are descriptors, never
gates. `prior_robustness` is the smallest `p12` in the sweep at which the call still holds
(`robust` = 1e-6, `intermediate` = 5e-6, `primary_prior` = 1e-5); the abf layer's own p12
arms supply it for fallback cells. `fallback_reason` says why a `coloc.abf` cell was not
scored by `coloc.susie`, at the first place it left the pipeline: a locus over the SNP guard
was never tested at signal level, while `gwas_no_credible_set` was tested and the GWAS did
not fine-map, which weakens any colocalization claimed there.

**SMR + HEIDI (10a) agrees with most nominations, as expected, and flags one disagreement.** It
uses the same GWAS and GTEx summary statistics as coloc, so agreement is a consistency check, not
independent evidence. **30/42** nominations have at least one tissue where the gene's
pre-designated primary sQTL probe is SMR-significant with HEIDI not rejected (SNCA/LBD 11/11,
UNC13A 2/2, PICALM 1/1); 31/42 counting either modality, and 17/42 if the probe must be the exact
intron coloc headlined in that tissue. Of the other 12, 9 have no instrument at 5e-8, 2 are
HEIDI-rejected and 1 is instrumented but null. *(An earlier "29/42" here matched none of these
definitions; recounted 2026-09-12 from `smr_heidi/gtex/smr_results.parquet`, keyed on gene × trait.)* **SNCA/PD disagrees, and the disagreement is a HEIDI rejection rather than an
absence of SMR signal.** All 8 sQTL tissues are instrumented; 6 clear the family Bonferroni
threshold (0.05/20 = 2.5e-3) at p_SMR 1.1e-3 to 2.2e-3, and HEIDI rejects in every one of them.
The other 2 fall just short of the threshold (p_SMR 2.5e-3 and 4.2e-3) and HEIDI rejects there
too: p_HEIDI is below 1e-6 in all 8, median 6.9e-11. The eQTL arm has no instrument at 5e-8 in
any of the 8. Report it as a disagreement; it does not overrule the colocalization.
The HEIDI cut is 0.01 as prespecified; 20 of 113 not-rejected primary tissue probes would be
rejected at 0.05, and ZNF232/AD and GPM6A/SCZ would lose every tissue.

Two gates gate the top tier, and they are independent. `status: reviewed` means the
coordinates were verified here. `evidence_class` must be `functional_validation` or
`coloc_association` — **a TWAS anchor may not promote**, because TWAS does not establish a
shared causal variant and LD-driven co-regulation produces the same signal. PICALM is
therefore capped at `disease_locus_splice_linked` however well its coordinates match.

A fourth tier, `context_distinct_splice_colocalization`, sits between the top tier and
`disease_locus_splice_linked` and is reachable **only from the all-introns arm**
(`--sqtl-arm all`). It applies when the curated event was itself *tested* and did not
colocalize while a different intron of the same gene did. The distinguishing fact is
measured, not asserted: `curated_event_testability.parquet` reads the GTEx all-pairs
release to establish that a phenotype exists at the curated intron. Without it a
non-match is uninformative — the event might colocalize and we would not know — which is
the whole difference between the two tiers. Where the curated event WAS on trial and
stayed silent, the observed signal is demonstrably not a proxy for it, and that is event
specificity rather than a near miss.

`_m/{coloc,ldsc,gwas}/` per-trait subdirectories are gitignored (235 GB / 1.5 GB of LD
and SuSiE intermediates); only the lean summaries and result tables are tracked.

**CLIs:** `isograph_benchmark/real_data/{qtl_anchoring,qtl_anchoring_meta,sqtl_concordance,sqtl_concordance_meta,module_genetic_anchoring,module_coloc_convergence,coloc_*,ldsc_*,locus_event_audit,brainseq_switch_qtl,brainseq_qtl_checks,smr_heidi}.py`, `gene_deep_dive.py --part events`, `isograph_benchmark/gwas/`.

## Display items

Main **Fig 3** `figQtlSpecificity`, **Fig 4** `figGeneticAnchoring` (rendered by
`08_integration/_h/01a`), **Table 1** (`table2_qtl_specificity_contrast`); S-real-2
`figGwasResolution`; supplementary tables S3–S5. **Table 2** (`table3_splicing_led_genes`),
`figSczConvergence` / Fig 4E and S8–S12 come from `08_integration/`.
