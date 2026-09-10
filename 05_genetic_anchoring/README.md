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
| 25 | `brainseq_switch_qtl` | BrainSEQ cis-swQTL on `S_g` and matched cis-eQTL on `A_g`, same donors/variants/covariates (GPU) |
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
  set-level contrast and S-LDSC, not on per-locus colocalization. Frame the 12
  splicing-led genes as resolved candidates whose strength is cross-disease and
  GO-invisible coherence.
- Scope: cis-sQTL anchors member-gene *splicing* to genetics, not the co-switching
  coordination itself.

## Signal-level coloc and the event audit — what they changed

Two things a reader should know before quoting a per-locus number.

**Half the grid cannot be scored at signal level.** `coloc.susie` needs credible sets on
both sides. Only **234 of 528** target loci fine-map on the GWAS side, so the rest fall
back to `coloc.abf` and carry `estimator = "abf"`. That is a property of the GWAS, not a
bug: PICALM and CTSH have **zero** GWAS credible sets at their loci, so signal-level coloc
cannot speak to them at all. Quote the fallback rate whenever the hierarchy is described.

**SNCA survives, and the multi-signal case is real.** At SNCA/LBD, SuSiE resolves **two**
QTL credible sets and only one colocalizes (PP4 0.974 on rs6532192; the other is
PP3 0.999). That is exactly the situation `coloc.abf`'s single-causal-variant assumption
cannot represent, and it holds across the whole `p12` sweep (0.79 at 1e-6 to 0.997 at
1e-4).

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
grouped-permutation representative intron, so the representative arm never tested it. That
is a limitation of the nomination design, not evidence against the loci; the all-introns
arm (`COLOC_SIGNAL_SQTL=all`) is what settles it.

That arm writes under `_m/coloc_signal_susie/all_introns/` and is read back with
`--stage meta --sqtl all`. Both modes name their shards `<analysis>__<tissue>.parquet`, so
the directory split is what keeps an all-introns run from overwriting the representative
results it exists to be compared against.

Two gates gate the top tier, and they are independent. `status: reviewed` means the
coordinates were verified here. `evidence_class` must be `functional_validation` or
`coloc_association` — **a TWAS anchor may not promote**, because TWAS does not establish a
shared causal variant and LD-driven co-regulation produces the same signal. PICALM is
therefore capped at `disease_locus_splice_linked` however well its coordinates match.

`_m/{coloc,ldsc,gwas}/` per-trait subdirectories are gitignored (235 GB / 1.5 GB of LD
and SuSiE intermediates); only the lean summaries and result tables are tracked.

**CLIs:** `isograph_benchmark/real_data/{qtl_anchoring,qtl_anchoring_meta,sqtl_concordance,sqtl_concordance_meta,module_genetic_anchoring,coloc_*,ldsc_*,gene_deep_dive,scz_age_projection}.py`, `isograph_benchmark/gwas/`.

## Display items

Main **Fig 3** `figQtlSpecificity`, **Fig 4** `figGeneticAnchoring`, **Table 1**
(`table2_qtl_specificity_contrast`), **Table 2** (`table3_splicing_led_genes`);
S-real-2 `figGwasResolution`, `figSczConvergence`; supplementary tables S3–S5, S8–S12.
