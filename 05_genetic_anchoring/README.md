# 05 — Genetic anchoring

**Question:** are the co-switch modules genetically real, and does disease risk resolve
onto specific isoform switches? This stage carries the paper's headline.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–02 | `qtl_anchoring`, `qtl_anchoring_matched` | Power-matched logistic sQTL/eQTL enrichment per analysis (17 analyses) |
| 18 | `module_coloc_convergence` | Module-level coloc convergence for AD/PD/LBD/ALS/SCZ: size-matched permutation null + anchored-module enrichment against the CLPP-tested pool. Generalises the SCZ-only layer in `scz_age_projection.py` |
| 17 | `qtl_anchoring_sensitivity` | Pre-specified sensitivity arms: constraint-adjusted (gnomAD LOEUF + missense z + log expression, both covariate sets on the identical constraint-complete subset), threshold-free continuous (rank-INT of -log10 pval_beta), and SuSiE credible-set dose. Writes to its own files; the primary arm from 01 is untouched |
| 03 | `sqtl_concordance` | Switch-transcript → LeafCutter-intron direction concordance |
| 04 | `module_anchoring` | Module-level genetic anchoring |
| 05–07 | `prep_module_gene_sets`, `run_magma`, `plot_magma` | MAGMA module-GWAS enrichment (`configs/gwas_magma.yaml`) |
| 08–11 | `coloc_prep`, `locus_ld`, `coloc_clpp`, `coloc_direction` | GTEx v11 SuSiE × 5 GWAS colocalization (eCAVIAR CLPP) |
| 12–14 | `ldsc_annot_prep`, `ldsc_make_annot_ldscores`, `ldsc_munge_h2` | S-LDSC partitioned heritability (baselineLD v2.2) |
| 22–23 | `coloc_gwas_susie`, `coloc_signal_susie` | **Signal-level colocalization.** Caches the per-locus GWAS SuSiE once (22), then re-fits each gene's cis QTL with `susie_rss` and runs `coloc::coloc.susie` signal-against-signal (23). Estimator hierarchy is **susie > abf > CLPP**; cells SuSiE cannot speak to fall back to `coloc.abf` and are flagged, never dropped |
| 24 | `locus_event_audit` | **Does a colocalizing locus name the event it claims to?** GWAS signal → sQTL signal → intron phenotype → driver transcript, matched against `configs/known_splice_events.yaml` |
| 25 | `brainseq_switch_qtl` | BrainSEQ cis-swQTL on `S_g` and matched cis-eQTL on `A_g`, same donors/variants/covariates (GPU). **The first all_samples swQTL results (2026-09-10) are superseded:** `build_phenotypes` omitted the discovery transcript filter, so `S_g` was not the discovery switch coordinate (median \|r\| 0.31-0.35 on shared libraries; 1.000 with the filter). Re-mapped 2026-09-11; `A_g` was unaffected. `ea_only` uses within-EA genotype PCs |
| 28 | `brainseq_qtl_checks` | Gates before any BrainSEQ QTL result is read: `S_g` must reproduce the discovery coordinate (median \|r\| >= 0.99) and is sign-pinned to it; `A_g` eQTL positive control against tissue-matched GTEx v11 eGenes (pi1 >= 0.5 and direction concordance >= 0.90 at GTEx lead variants, both pre-specified). all_samples control PASS in all three regions (pi1 0.83-0.88, concordance 0.987-0.991) |
| 26 | `locus_ld_robustness` | Locus LD audit for any locus carrying a biological claim: convergence, credible-set purity, boundary trims, `coloc.susie` across the p12 sweep, and LD-mismatch diagnostics (`kriging_rss` outlier drop, residual-variance refit) |
| 27 | `smr_heidi` | **SMR + HEIDI** on the signal-level nominations (GTEx v11 QTL, 1000G EUR LD). Orthogonal corroboration beneath coloc: `b_SMR` is not causal direction, a non-rejected HEIDI is not proof of sharing, and a HEIDI rejection does not overrule coloc. BrainSEQ QTL are refused until an `ea_only` mapping exists |
| 15 | `build_deep_dive` | Per-gene deep dive: anchor → switch → consequence |
| 16 | `scz_age_projection` | Are age-sensitive switch programs disrupted in SCZ? (the *convergence* sub-test is RETRACTED 2026-08-30 — wrong background; see `module_coloc_convergence`) |

## The headline, stated honestly

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
  which are null everywhere (0.959–1.021, p=0.30–0.86).
- **Per-gene CLPP posteriors are individually modest** (inclusion at ≥0.01; only 4/12
  reach 0.05; only CTSH at 0.39). The genetic-anchoring significance rests on the
  set-level contrast and S-LDSC, not on per-locus colocalization. The 12 splicing-led genes
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
recovered locus enters the narrative only after `26.locus_ld_robustness.sh`. For PICALM
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

**SMR + HEIDI (27) agrees with most nominations, as expected, and flags one disagreement.** It
uses the same GWAS and GTEx summary statistics as coloc, so agreement is a consistency check, not
independent evidence. 29/42 nominations have at least one tissue where the colocalizing intron is
SMR-significant with HEIDI not rejected (SNCA/LBD 11/11, UNC13A 2/2, PICALM 1/1); 9 have no
instrument at 5e-8. **SNCA/PD disagrees**: 7/8 tissues are not significant after Bonferroni (min
p_SMR 1.1e-3) and HEIDI rejects in the eighth (median p_HEIDI 7e-11). Report it as a
disagreement; it does not overrule the colocalization. The HEIDI cut is 0.01 as prespecified; 20
of 113 not-rejected primary tissue probes would be rejected at 0.05, and ZNF232/AD and GPM6A/SCZ
would lose every tissue.

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

**CLIs:** `isograph_benchmark/real_data/{qtl_anchoring,qtl_anchoring_meta,sqtl_concordance,sqtl_concordance_meta,module_genetic_anchoring,coloc_*,ldsc_*,gene_deep_dive,scz_age_projection}.py`, `isograph_benchmark/gwas/`.

## Display items

Main **Fig 3** `figQtlSpecificity`, **Fig 4** `figGeneticAnchoring`, **Table 1**
(`table2_qtl_specificity_contrast`), **Table 2** (`table3_splicing_led_genes`);
S-real-2 `figGwasResolution`, `figSczConvergence`; supplementary tables S3–S5, S8–S12.
