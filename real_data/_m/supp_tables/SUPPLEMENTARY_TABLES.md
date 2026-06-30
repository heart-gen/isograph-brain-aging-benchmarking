# Real-data supplementary tables

Manubot-ready legends for the real-data supplementary tables. Each table is a clean CSV in
this directory, regenerated from the committed analysis parquet ledgers by
`real_data/_h/assemble_supp_tables.py` (login node, no SLURM). Every number is copied
verbatim from the source ledger — do not hand-edit a CSV; rerun the assembler. Each table is
tied to one honest claim from the real-data analysis spine (AGENTS.md §§1–4).

| Table | File | Source ledger | Honest claim it supports |
| --- | --- | --- | --- |
| S1 | `tableS1_baseline_pooled.csv` | `baseline_comparison/baseline_comparison_pooled.parquet` | IsoGraph is **not** globally superior on module metrics |
| S2 | `tableS2_baseline_per_region.csv` | `baseline_comparison/baseline_comparison.parquet` | per-region ledger backing S1 |
| S3 | `tableS3_qtl_specificity_contrast.csv` | `qtl_anchoring_meta/qtl_anchoring_meta_contrast.parquet` | splicing genetics spared in disease/GO-invisible IsoGraph modules |
| S4 | `tableS4_qtl_specificity_matched_baseline.csv` | `qtl_anchoring_meta/qtl_anchoring_meta_contrast_common.parquet` | the specificity effect is IsoGraph-only (clean method effect) |
| S5 | `tableS5_qtl_raw_enrichment_or.csv` | `qtl_anchoring_meta/qtl_anchoring_meta.parquet` | shared cis-QTL depletion baseline the ratio removes |
| S6 | `tableS6_go_invisible_gate.csv` | `brainseq/caudate_sczd/_m/go_invisible_gate.parquet` | disease switch modules are real DTU-without-DGE biology |
| S7 | `tableS7_module_trust_funnel.csv` | `stability/_m/module_trust/*` | per-module trust scaffold (stability→drivers→replication→complementarity) |

---

## Table S1 — Three-baseline module comparison (pooled per method)

Pooled per-module rates across 17 analyses (1 SCZD + 3 BrainSEQ aging + 13 GTEx aging) for
the four module sets fit on the same samples: `isograph` (VAE+Leiden on switch+abundance),
`wgcna_switch_only` and `wgcna_multiplex` (classical WGCNA on the **identical** switch / 
switch+abundance features), and `wgcna_gene` (classical WGCNA on abundance only).
Columns: median module count and size, per-module phenotype-significant rate (`pheno_fdr ≤
0.1`), "both" rate (phenotype-sig AND GO-enriched), GO-enriched rate (`n_go_terms > 0`), and
raw totals. **Read rates, not totals** — totals scale with module count (IsoGraph runs finer:
median 35 vs 8–18.5 modules). The phenotype-sig rate is highest for the two switch-fed
methods (`wgcna_switch_only` 0.336, `isograph` 0.268) over the abundance-fed ones
(`wgcna_multiplex` 0.189, `wgcna_gene` 0.180): phenotype sensitivity comes from the switch
features, not the inference method. GO enrichment is abundance-dominated (`wgcna_gene` 0.885 ≫
`isograph` 0.217). The one clean method effect: on identical multiplex features `isograph`
(0.268) > `wgcna_multiplex` (0.189). Supports: IsoGraph's value is DTU-without-DGE content,
not better module-level enrichment.

## Table S2 — Three-baseline comparison, per cohort × region

The 66-row per-analysis ledger behind S1 (cohort, region, method, features, module count and
median size, GO-enriched / phenotype-sig / both counts and fractions). Use to confirm the
pooled rates are not driven by a single region and to read the per-region module-count
disparity that makes raw totals non-comparable.

## Table S3 — QTL splicing-specificity contrast (all-tissue meta)

Random-/fixed-effects meta-analysis of the **paired sQTL-OR / eQTL-OR ratio** within each
analysis, by module set and graph method. The ratio cancels the shared cis-QTL depletion of
constrained network genes (Table S5), isolating whether *splicing* genetics is spared. For
`isograph`, the ratio is > 1 and significant in all-modules (1.07, p=8e-6), phenotype-sig
(1.13, p=1.5e-4) and GO-invisible (1.13, p=2.6e-3) sets, and null for GO-visible (1.04, ns) —
an internal control. The matched WGCNA baselines (`wgcna_switch_only`, `wgcna_multiplex`) are
null in every set (ratios 0.97–1.02). Columns: method, module set, n analyses, fixed-effects
ratio + 95% CI, p, and heterogeneity I². Supports: splicing genetics is spared exactly in the
disease/GO-invisible IsoGraph modules, and only for IsoGraph.

## Table S4 — Matched-baseline QTL specificity contrast (8-tissue common set)

The same contrast restricted to the 8 tissues where all three graph methods have a result, so
`isograph` vs `wgcna_switch_only` vs `wgcna_multiplex` is a like-for-like method comparison on
identical tissues. `isograph` stays positive (pheno-sig 1.11, p=2.1e-3; GO-invisible 1.11,
p=0.011); both WGCNA baselines stay null (0.97–1.02). This is the clean method effect — the
single-tissue `wgcna_switch_only` frontal-cortex blip (1.29) does not survive pooling.
Supports: the splicing-specificity signal is a property of IsoGraph's inference, not the
switch features alone.

## Table S5 — Raw cis-QTL enrichment ORs (eQTL & sQTL meta)

Per-method, per-QTL-kind meta ORs (foreground module genes vs background) by module set. Both
eQTL and sQTL ORs are < 1 everywhere (e.g. `isograph` all-modules eQTL 0.82, p=5e-122; sQTL
0.88, p=1.6e-25): network/module genes are cis-QTL-depleted, as expected for constrained
genes. This shared depletion is the baseline the S3/S4 **ratio** removes — the point of the
specificity contrast is that sQTL depletion is *weaker* than eQTL depletion in the disease
modules, not that either is enriched. Supports: the raw-OR context that prevents
mis-reading the specificity ratio as enrichment.

## Table S6 — GO-invisible SCZD disease switch modules (BrainSEQ caudate)

Per-module ledger for the four SCZD-associated IsoGraph switch modules (`pheno_fdr ≤ 0.1`),
plus a pooled `_background` row over all switch transcripts. All four return zero GO terms
(`go_invisible = True`), yet nearly every member carries a real anticorrelated transcript
pair (e.g. M010 73/83, M026 21/23), with max switch strength 1.10–1.39 and 93–459 significant
switch transcripts. Driver functional-consequence fractions (CDS / coding-status / biotype /
UTR change) sit at or above the pooled background (0.84 / 0.67 / 0.74 / 0.61), so
GO-invisibility reflects GO's gene-level/abundance bias, not low module quality. Supports: the
disease switch signal is genuine, functionally consequential isoform regulation invisible to
pathway enrichment (biology gate PASS, complementary form).

## Table S7 — Per-region module trust funnel

One row per region × method (IsoGraph) summarising the four-question funnel: Q1 stability
(`n_modules`, `n_trusted`, `frac_trusted` above a size-matched permutation null at FDR<0.05),
Q2 drivers (`median_driver_rho`, `frac_positive_rho` of shared-gene switch-loading Spearman ρ
across split halves), Q3 cross-cohort aging replication (`n_replicates`/`n_replication_pairs`,
BrainSEQ↔GTEx, sign-concordant and jointly significant; BrainSEQ rows only), and Q4
complementarity (`median_frac_dtu_without_dge`, `median_frac_in_wgcna_age`). Across six
regions: 236/266 modules trusted; driver ρ medians 0.77–0.82 with 96–100% positive; 25/130
BrainSEQ modules replicate cross-cohort (caudate 3/45, DLPFC 9/35, hippocampus 13/50).
Supports: the per-module trust scaffold the biological claims rest on.

---

## Reproducibility

- Assembler: `real_data/_h/assemble_supp_tables.py`
  (`python real_data/_h/assemble_supp_tables.py`; login node, no SLURM).
- Inputs: the committed analysis parquet ledgers named in the table above.
- Outputs: `tableS1`–`tableS7` CSVs in this directory.
- Environment: project Python 3.12 (`/ocean/projects/bio260021p/shared/opt/envs/isograph`),
  pandas/pyarrow.
- The assembler only re-shapes and rounds; it computes no new statistics. Regenerate after any
  re-run of the upstream analyses (baseline_comparison, qtl_anchoring, go_invisible_gate,
  module_trust).
