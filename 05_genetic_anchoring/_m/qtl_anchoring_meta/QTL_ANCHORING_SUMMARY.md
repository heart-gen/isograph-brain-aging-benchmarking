# IsoGraph co-switch module genetic anchoring (sQTL/eQTL specificity)

Modular analysis summary for Manubot integration. Generated from
`05_genetic_anchoring/_m/qtl_anchoring_meta/` and the per-analysis `qtl_anchoring.parquet` tables.
Every numeric claim is reproduced from the parquet files named in each section; do not
edit the numbers by hand — regenerate from the tables.

## Purpose

Test whether IsoGraph's co-switch modules are anchored to *splicing* genetics in a way
that (a) tracks the phenotype-associated modules and (b) is a property of IsoGraph's
representation + inference rather than of the switch features alone. This is the
orthogonal genetic validation of the module-trust scaffold: if the co-switch modules
merely re-described abundance structure, their member genes would carry no
splicing-specific genetic signal beyond expression-QTL, and a matched WGCNA baseline on
the same features would behave identically.

The GO-invisible / GO-visible split is retained below as a **secondary, exploratory**
stratification. It was originally framed as the primary localisation; on the 2026-08-29
refresh the two arms are indistinguishable, so it is reported, not claimed.

## Inputs

- **xQTL universes** — GTEx v11 brain `sGenes` / `eGenes` per tissue
  (`inputs/raw/gtex_v11/xqtl/`); the tested universe per tissue defines both the QTL
  status label and the background.
- **Module gene sets** — production IsoGraph `isograph_vae/` modules and the matched
  WGCNA baselines `wgcna_switch_only/` + `wgcna_multiplex/` (same switch / switch+abundance
  features), for 17 analyses (1 SCZD + 3 BrainSEQ aging + 13 GTEx aging).
- **Module sets** — `all_modules`, `pheno_sig_modules`, `go_invisible_modules`,
  `go_visible_modules` (the last two from the GO-invisible gate), so the contrast can be
  read where DTU-without-DGE concentrates vs the immune/abundance GO-visible control.
- **Power covariates** — per gene: cis-variant count, gene length, isoform count, and
  (where available) intron-group size, so enrichment is not confounded by testable-variant
  density.

## Methods text

Within each region/cohort analysis we tested whether co-switch module membership predicts
cis-QTL status among the genes in that tissue's xQTL-tested universe, separately for
splicing-QTL (sQTL) and expression-QTL (eQTL). Enrichment was estimated by power-matched
logistic regression (QTL status ~ module membership + log cis-variant count + log gene
length + log isoform count [+ log intron-group size]), giving a module-membership log
odds-ratio per (analysis, QTL kind, module set, graph method). Log odds-ratios were pooled
across analyses by inverse-variance meta-analysis under both fixed-effect and
DerSimonian–Laird random-effects models, with Cochran's Q and I² for heterogeneity. The
primary estimand is the **splicing-specificity contrast** — the paired within-analysis
ratio of the sQTL odds-ratio to the eQTL odds-ratio — which removes the shared cis-QTL
depletion baseline of constrained network genes and isolates whether splicing genetics is
spared relative to expression genetics. The identical pipeline was applied to the matched
WGCNA baselines on the same switch features and the same 13 GTEx tissues, so that comparing
the contrast across methods isolates feature-driven from inference-driven signal. The SCZD
analysis and the per-tissue tables additionally report the raw odds-ratios. Analyses used
IsoGraph v0.1.5 under the project Python 3.12 environment with fixed seeds.

## Results text

**Raw cis-QTL enrichment — a shared depletion baseline, not the result.** Pooled across 17
analyses, IsoGraph co-switch module genes are cis-QTL *depleted* for both QTL types
(`all_modules` eQTL OR 0.808, p=7.3e-138; sQTL OR 0.867, p=9.1e-31 —
`qtl_anchoring_meta.parquet`), the expected signature of coordinated/constrained network
genes carrying fewer common-variant cis-QTL. This shared OR < 1 baseline is the same for
expression and splicing and is **not** the estimand; it is reported only to motivate the
contrast.

**Splicing-specificity contrast — the result.** The paired sQTL-OR / eQTL-OR ratio exceeds
1 and is carried by the **phenotype-associated** modules
(`qtl_anchoring_meta_contrast.parquet`): `all_modules` **1.068** (95% CI 1.037–1.100,
p=1.3e-5, I²=0.29, Q=22.6, k=17), `pheno_sig_modules` **1.111** (1.048–1.177, p=3.6e-4,
I²=0.23, Q=14.3, k=12), `go_invisible_modules` 1.068 (0.993–1.148, **p=0.077, n.s.**,
I²=0.00, Q=9.4, k=11), and `go_visible_modules` 1.084 (1.000–1.174, p=0.050, I²=0.28,
Q=14.0, k=11). Splicing-QTL are spared ~7–11% relative to expression-QTL in co-switch
genes.

**The effect does NOT localise to the GO-invisible modules.** GO-invisible (1.068,
p=0.077) and GO-visible (1.084, p=0.050) are statistically indistinguishable — overlapping
CIs, and if anything GO-visible is the nominally larger of the two. Neither is a null and
neither is the site of the effect: it lives in the phenotype-associated set. The
GO-invisible arm's I²=0.00 is now a *homogeneous null*, not evidence of a consistent
effect, so the old "GO-invisible is homogeneous while GO-visible rides on heterogeneity"
contrast (I²=0.00 vs 0.68) is gone — on the refreshed inputs GO-visible's I² is 0.28.
**IsoGraph's DTU-without-DGE claim rests on the GO-invisible content gate
(`GO_INVISIBLE_GATE_SUMMARY.md`), not on these genetics.**

> **Provenance — never quote a pre-2026-08-29 number from this analysis.** An earlier
> version of this file reported `pheno_sig` 1.163 (p=3.6e-7) and `go_invisible` 1.172
> (p=2.3e-5, I²=0.00) and framed the latter as the headline. Those came from a
> `module_enrichment` table written *before* the fit it was joined to. Leiden re-assigns
> module ids on every fit (median same-id Jaccard between two cortex fits = 0.000), so the
> join did not merely use stale numbers — it handed each module another module's
> `pheno_fdr` and `n_go_terms`, meaning the GO-invisible/GO-visible partition the old
> contrast rested on was a relabelling. `all_modules` is bit-identical across the refresh,
> which is the control proving the pipeline itself did not change. Guarded since
> 2026-08-30 by `isograph_benchmark/real_data/partition_provenance.py`.

**Matched-baseline method effect — only IsoGraph shows it.** On the 13 GTEx tissues all
methods share (8-tissue common-subset contrast, `qtl_anchoring_meta_contrast_common.parquet`),
IsoGraph reproduces the specificity (`all_modules` **1.105**, p=5.8e-6, I²=0.00;
`pheno_sig` **1.108**, p=0.001, I²=0.00; `go_visible` 1.076, p=0.078; `go_invisible` 1.066,
p=0.11, n.s.) while the matched WGCNA baselines on identical switch features are **null
everywhere**: `wgcna_switch_only` pheno_sig 0.968 (p=0.46), go_invisible 1.021 (p=0.71),
go_visible 0.959 (p=0.30); `wgcna_multiplex` pheno_sig 0.980 (p=0.41), go_invisible 0.993
(p=0.86), go_visible 0.993 (p=0.77). **This is the primary internal control** — it holds
the switch features fixed and varies only the inference, which localises the effect far
more sharply than a module-content contrast can. Same features + classical inference loses
the splicing-genetic signal that IsoGraph's VAE + Leiden concentrates in the
phenotype-associated modules — the genetic-anchoring analog of the three-baseline result.
(One tissue, frontal_cortex, looked specific for `wgcna_switch_only` at 1.29 but does not
survive pooling; trust the meta, not one tissue.)

**Headline:** *Splicing-QTL are spared ~11% over expression-QTL in IsoGraph's
phenotype-associated co-switch modules (1.111, 95% CI 1.048–1.177, p=3.6e-4), and the
effect is absent in matched WGCNA baselines built on identical switch features — so
genetic anchoring of the isoform-regulation layer is an IsoGraph method effect, not a
property of the switch features. The effect is not specific to the GO-invisible modules.*

## Figure and table notes

- **Main figure — splicing-specificity contrast
  (`manuscript/_m/figures/figQtlSpecificity.{pdf,png}`, built by
  `manuscript/_h/qtl_specificity_figure.R`).** Three panels, no in-panel titles:
  - **(A)** raw pooled cis-QTL ORs for IsoGraph by module set, sQTL vs eQTL, reference
    line at OR = 1 — both depleted, sQTL above eQTL (motivates the contrast).
  - **(B)** splicing-specificity contrast (sQTL-OR / eQTL-OR) per module set for IsoGraph,
    forest-style point + 95% CI, reference line at ratio = 1 — carried by pheno-sig; the
    GO-invisible and GO-visible arms overlap each other and the caption must not present
    either as the site of the effect or as a null.
  - **(C)** matched-baseline method effect on the shared tissues: the contrast per module
    set coloured by method (IsoGraph / wgcna_switch_only / wgcna_multiplex) — only IsoGraph
    is positive; baselines straddle 1.
  - Key message: the splicing-specific genetic anchoring is carried by the
    phenotype-associated modules and is a method effect, not a feature effect.
- **Supplementary table:** `qtl_anchoring_meta.parquet`,
  `qtl_anchoring_meta_contrast{,_common}.parquet`, and `per_analysis.parquet` — full
  pooled ORs, contrasts, and per-tissue ledger (graph_method, module_set, xqtl_kind, k,
  ORs/ratios with CI, p_fe/p_re, Q, I²).

## Reproducibility information

- Analysis directory: `05_genetic_anchoring/_m/qtl_anchoring_meta/` (meta) and per-region
  `02_module_discovery/{brainseq,gtex}/<region>/_m/qtl_anchoring*.parquet`.
- Primary scripts: `isograph_benchmark/real_data/qtl_anchoring.py` (per-analysis, with
  `--method {isograph,wgcna_switch_only,wgcna_multiplex}`) and
  `isograph_benchmark/real_data/qtl_anchoring_meta.py` (IVW FE + DL RE pooling and
  contrasts).
- SLURM drivers: `05_genetic_anchoring/_h/01.qtl_anchoring.sh` (17-analysis array) and
  `05_genetic_anchoring/_h/02.qtl_anchoring_matched.sh` (2 methods × 13 GTEx tissues).
- Key parameters: power-matched logistic enrichment with cis-variant count / gene length /
  isoform count / intron-group size covariates; FDR 0.05; FE + DL-RE pooling; deterministic
  seeds.
- Outputs: `qtl_anchoring_meta.parquet`, `qtl_anchoring_meta_contrast.parquet`,
  `qtl_anchoring_meta_contrast_common.parquet`, `per_analysis.parquet`,
  `QTL_ANCHORING_META.md`; per-region `qtl_anchoring.{parquet,json}` + `QTL_ANCHORING.md`.
- Git commit: aligned with the matched-baseline run (2026-06-26).
- Compute environment: PSC Bridges-2 RM-shared; IsoGraph v0.1.5, project Python 3.12
  environment (`/ocean/projects/bio260021p/shared/opt/envs/isograph`).
- Missing reproducibility information: per-package versions are taken from the live
  environment, not a per-run lockfile.

## Limitations and integration notes

- The estimand is the **contrast**, not the raw OR. Co-switch genes are cis-QTL depleted
  for both QTL types (constraint baseline); reporting raw sQTL enrichment alone would be
  misleading, and the favorable single-tissue SCZD sQTL tail (OR 1.37) does not survive
  pooling — the honest signal is the splicing-vs-expression specificity, not absolute
  enrichment.
- The contrast standard error is conservative (it treats the sQTL and eQTL estimates as
  independent though they share the foreground genes), so the p-values are if anything
  understated.
- Scope: cis-sQTL anchors *member-gene splicing* to genetics, not the co-switching
  coordination itself; the claim is genetic anchoring of the modules' splicing biology, not
  a genetic basis for the network structure.
- High I² in some pooled rows flags between-tissue heterogeneity; prefer the RE estimate
  there. The headline `pheno_sig` contrast is low-heterogeneity (I²=0.23, and I²=0.00 on
  the 8-tissue common subset).
- **The GO-invisible modules are no longer a genetic subgroup claim.** They remain the
  content result (GO-invisible gate) and the set the coloc deep-dive draws on, but this
  analysis does not show splicing genetics concentrating in them. Any integration text must
  keep the content claim and the genetic claim on separate evidence.
- This is a **primary** orthogonal-validation analysis. Integrate with the module-trust
  funnel (the modules are reproducible), the GO-invisible gate (their biology is real
  isoform switching), and the three-baseline comparison (IsoGraph is not globally superior
  on module rates) — together these supply the complementary-DTU-without-DGE Results arc,
  with this analysis as the genetic-mechanism anchor.
