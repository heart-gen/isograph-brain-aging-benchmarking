# Real-data supplementary tables

Manubot-ready legends for the real-data supplementary tables. Each table is a clean CSV in
this directory, regenerated from the committed analysis parquet ledgers by
`manuscript/_h/assemble_supp_tables.py` (login node, no SLURM). Every number is copied
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
| S8 | `deep_dive/deep_dive_panel.tsv` | `deep_dive/deep_dive_panel.parquet` | per-gene verdict for every colocalized gene (splicing-led vs expression-led) |
| S9 | `deep_dive/deep_dive_events.tsv` | `coloc/_m/coloc_isoform_events_combined.parquet` (+direction) | per-event anchor→switch→consequence for every colocalized gene |
| S10 | `deep_dive/deep_dive_rbp.tsv` | `_m/rbp/{rbp_switch_calls,rbp_regulon}.parquet` | per-gene switched + module-enriched RBP regulators |
| S11 | `deep_dive/deep_dive_exon_clinical.tsv` | per-region `clinical_consequence/exon_clinvar.parquet` | per-gene/exon switched-vs-constitutive ClinVar & CDS annotation |
| S12 | `deep_dive/deep_dive_literature.tsv` | `deep_dive/deep_dive_literature.parquet` | curated known-isoform-biology literature per resolved splicing-led gene (with Manubot citekeys) |

Tables S8–S12 (the per-gene deep-dive) live under `05_genetic_anchoring/_m/deep_dive/` and are regenerated
by `05_genetic_anchoring/_h/15.build_deep_dive.sh` (not the `assemble_supp_tables.py` assembler). Together they
let a reader reconstruct the SNCA-style mechanistic vignette for any colocalized gene without a
hand-written narrative: S9 gives the variant→junction→switch-pair→structural-consequence chain,
S10 the candidate RBP regulators, S11 the clinical (coding-vs-non-coding, ClinVar) read, S8
the one-line verdict, and S12 the known-isoform-biology literature (four genes with documented
disease isoform biology, eight flagged `novel_candidate`).

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
across split halves), Q3 cross-cohort aging concordance (`n_concordant`/`n_replication_pairs`,
BrainSEQ↔GTEx, sign-concordant and jointly significant in both cohorts; BrainSEQ rows
only; counted on the linear covariate-free arm — see REPLICATION_PERMUTATION.md for the
permutation null and the linear-vs-spline asymmetry), and Q4
complementarity (`median_frac_dtu_without_dge`, `median_frac_in_wgcna_age`). Across six
regions: 236/266 modules trusted; driver ρ medians 0.77–0.82 with 96–100% positive; 25/130
BrainSEQ modules replicate cross-cohort (caudate 3/45, DLPFC 9/35, hippocampus 13/50).
Supports: the per-module trust scaffold the biological claims rest on.

## Table S8 — Per-gene deep-dive panel (all colocalized genes)

One row per colocalized disease gene (n = 68), ranked splicing-led first. Columns: anchoring
trait(s), QTL kinds (sQTL/eQTL), max eCAVIAR CLPP and its tissue, gnomAD LOEUF, number of
resolved switch-pair events and a multi-locus flag, concordant trait(s), BrainSeq replication and
GO-invisible flags, and the verdict (splicing-led = an sQTL resolving onto a concordant IsoGraph
switch pair; expression-led = eQTL gene-level only; splicing-unresolved = sQTL not mapping onto
the switch pair). Of 68 genes, 12 are splicing-led — **all 12 GO-invisible** — 23
splicing-unresolved, and 33 expression-led; four splicing-led genes are multi-locus (PPP6R2,
SNCA, CDIP1, DLG1). Supports: the genetically-anchored isoform-switch set and its
verdict-by-gene provenance.

## Table S9 — Per-event genetic anchor → switch → consequence

One row per colocalized isoform event (n = 141 over 68 genes): gene, trait, case, QTL kind,
tissue, lead variant (rsID, ref/alt) and risk allele, signed risk-QTL effect and direction, the
LeafCutter junction, CLPP, GO-invisibility, the IsoGraph switch pair and whether the junction
maps into it (`junction_in_switch_pair`), the structural consequence, GTEx concordance, BrainSeq
region and replication, and a resolved-event sentence. This is the reader's reconstruction table:
filtering to a gene reproduces the variant→junction→switch chain used for the SNCA vignette.

## Table S10 — Per-gene RBP regulators

One row per gene × region × module × RBP (n = 973 over 54 genes) where an RBP motif is both
called switched in the gene and enriched in the gene's IsoGraph module at BH q < 0.05; columns
include the module, GO-invisibility, motif enrichment and q. Identifies candidate trans splice
regulators of each gene's switch (e.g. ADAR/CPEB2/PTBP2 for SNCA). Exploratory (motif presence,
not measured binding). Supports: the regulatory-logic layer of each vignette.

## Table S11 — Per-gene/exon clinical annotation

One row per gene × region × exon (n = 6,041 over 64 genes): coordinates and length, whether the
exon is differentially used (`switched`) or constitutive, whether it overlaps CDS
(`cds_overlap`), and its ClinVar pathogenic/likely-pathogenic (`n_plp`) and total (`n_clinvar`)
variant counts. Lets a reader do the SNCA-style clinical read for any gene — e.g. confirm that a
switched exon is non-coding and pathogenic-variant-free while the gene's P/LP burden sits in
shared constitutive coding exons. Supports: the clinical-consequence layer of each vignette.

## Table S12 — Per-gene literature (known isoform biology)

One row per resolved splicing-led gene (n = 12): a curated synthesis of known isoform biology in
the relevant disease (`literature`), Manubot citekeys (`references`), and a `curation` flag —
`documented` for the four genes with established disease isoform biology matching their resolved
switch (SNCA [@doi:10.3389/fgene.2019.00584; @doi:10.3390/genes9020063], DLG1/SAP97
[@doi:10.1038/tp.2015.154], CTSH [@doi:10.1038/s41386-023-01542-2], ARVCF
[@doi:10.1038/sj.mp.4001586]) and `novel_candidate` for the eight without established
disease-specific isoform literature (PPP6R2, GGNBP2, PGS1, CDIP1, PRRC2B, RTEL1, TBC1D15, TPCN1).
Curation is data (a dict in `gene_deep_dive.py`) so the table and vignette Section 6 regenerate
deterministically; no DOI/PMID is fabricated. Supports: the literature layer (layer 6) of each
resolved vignette and the deep-dive Results paragraph.

---

## Reproducibility

- Assembler: `manuscript/_h/assemble_supp_tables.py`
  (`python manuscript/_h/assemble_supp_tables.py`; login node, no SLURM).
- Inputs: the committed analysis parquet ledgers named in the table above.
- Outputs: `tableS1`–`tableS7` CSVs in this directory.
- Environment: project Python 3.12 (`/ocean/projects/bio260021p/shared/opt/envs/isograph`),
  pandas/pyarrow.
- The assembler only re-shapes and rounds; it computes no new statistics. Regenerate after any
  re-run of the upstream analyses (baseline_comparison, qtl_anchoring, go_invisible_gate,
  module_trust).
- Tables S8–S12 have a separate generator: `05_genetic_anchoring/_h/15.build_deep_dive.sh`, which runs
  `python -m isograph_benchmark.real_data.gene_deep_dive` (deterministic joins over the coloc,
  RBP, and clinical-consequence ledgers, plus the curated `_LITERATURE` dict) and writes
  `deep_dive_*.{tsv,parquet}` under
  `05_genetic_anchoring/_m/deep_dive/`. The same script extracts the SNCA transcript exons from GENCODE v47
  and renders the four genetic-anchoring figures. Regenerate after any re-run of the coloc,
  RBP-regulon, or clinical-consequence analyses.
