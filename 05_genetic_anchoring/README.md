# 05 — Genetic anchoring

**Question:** are the co-switch modules genetically real, and does disease risk resolve
onto specific isoform switches? This stage carries the paper's headline.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–02 | `qtl_anchoring`, `qtl_anchoring_matched` | Power-matched logistic sQTL/eQTL enrichment per analysis (17 analyses) |
| 03 | `sqtl_concordance` | Switch-transcript → LeafCutter-intron direction concordance |
| 04 | `module_anchoring` | Module-level genetic anchoring |
| 05–07 | `prep_module_gene_sets`, `run_magma`, `plot_magma` | MAGMA module-GWAS enrichment (`configs/gwas_magma.yaml`) |
| 08–11 | `coloc_prep`, `locus_ld`, `coloc_clpp`, `coloc_direction` | GTEx v11 SuSiE × 5 GWAS colocalization (eCAVIAR CLPP) |
| 12–14 | `ldsc_annot_prep`, `ldsc_make_annot_ldscores`, `ldsc_munge_h2` | S-LDSC partitioned heritability (baselineLD v2.2) |
| 15 | `build_deep_dive` | Per-gene deep dive: anchor → switch → consequence |
| 16 | `scz_age_projection` | Do SCZ-risk loci converge on age-sensitive switch programs? |

## The headline, stated honestly

Co-switch genes are cis-QTL **depleted for both** sQTL and eQTL — coordinated network
genes are constrained. That shared baseline is not the result. The result is the paired
**splicing-specificity contrast** (sQTL OR / eQTL OR within analysis), which removes it:
all 1.068 (p=1.3e-5), pheno-sig 1.163 (p=3.6e-7), **GO-invisible 1.172 (p=2.3e-5,
I²=0.00)**, GO-visible 1.104 (p=0.022, I²=0.68).

Two things not to overclaim:

- **GO-visible is not a clean internal null.** It is the low, heterogeneous end of a
  gradient. The primary internal control is the **matched WGCNA baselines** from stage 02,
  which are null everywhere (0.970–1.017, p=0.44–0.86).
- **Per-gene CLPP posteriors are individually modest** (inclusion at ≥0.01; only 4/12
  reach 0.05; only CTSH at 0.39). The genetic-anchoring significance rests on the
  set-level contrast and S-LDSC, not on per-locus colocalization. Frame the 12
  splicing-led genes as resolved candidates whose strength is cross-disease and
  GO-invisible coherence.
- Scope: cis-sQTL anchors member-gene *splicing* to genetics, not the co-switching
  coordination itself.

`_m/{coloc,ldsc,gwas}/` per-trait subdirectories are gitignored (235 GB / 1.5 GB of LD
and SuSiE intermediates); only the lean summaries and result tables are tracked.

**CLIs:** `isograph_benchmark/real_data/{qtl_anchoring,qtl_anchoring_meta,sqtl_concordance,sqtl_concordance_meta,module_genetic_anchoring,coloc_*,ldsc_*,gene_deep_dive,scz_age_projection}.py`, `isograph_benchmark/gwas/`.

## Display items

Main **Fig 3** `figQtlSpecificity`, **Fig 4** `figGeneticAnchoring`, **Table 1**
(`table2_qtl_specificity_contrast`), **Table 2** (`table3_splicing_led_genes`);
S-real-2 `figGwasResolution`, `figSczConvergence`; supplementary tables S3–S5, S8–S12.
