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
| S7 | `tableS7_module_trust_funnel.csv` | `stability/module_trust_tables/region_funnel.csv` | per-region funnel, **both methods** (stability→drivers→split-half aging→complementarity) |
| S7a | `tableS7a_split_half_module_ledger.csv` | `module_stability__*.parquet` | every module behind the 93/118 and 62/72 trusted counts, with its permutation null |
| S7b | `tableS7b_projection_module_ledger.csv` | `eigengene_projection/eigengene_projection_all.parquet` | every module behind the frozen-eigengene transfer counts, standardised **and** raw |
| S7c | `tableS7c_crosscohort_permutation.csv` | `replication_permutation__*__stats.json` | the matched-pair count under all three covariate modes × three statistics × two nulls |
| S7d | `tableS7d_functional_preservation.csv` | `functional_preservation__*__stats.json` | what matched cross-cohort pairs still share when the gene partition does not replicate |
| S7e | `tableS7e_resolution_sensitivity.csv` | `stability_summary.parquet` | split-half ARI/NMI per region across Leiden 0.5–20, production 2.0 flagged, WGCNA beside it |
| S7f | `tableS7f_projection_sign_scale.csv` | `eigengene_projection_all.parquet`, `switch_axis_alignment__*.parquet` | raw vs null-standardised projected-age sign agreement, with the switch-axis sign convention |
| S8 | `deep_dive/deep_dive_panel.tsv` | `deep_dive/deep_dive_panel.parquet` | per-gene verdict for every colocalized gene (splicing-led vs expression-led) |
| S9 | `deep_dive/deep_dive_events.tsv` | `coloc/_m/coloc_isoform_events_combined.parquet` (+direction) | per-event anchor→switch→consequence for every colocalized gene |
| S10 | `deep_dive/deep_dive_rbp.tsv` | `_m/rbp/{rbp_switch_calls,rbp_regulon}.parquet` | per-gene switched + module-enriched RBP regulators |
| S11 | `deep_dive/deep_dive_exon_clinical.tsv` | per-region `clinical_consequence/exon_clinvar.parquet` | per-gene/exon switched-vs-constitutive ClinVar & CDS annotation |
| S12 | `deep_dive/deep_dive_literature.tsv` | `deep_dive/deep_dive_literature.parquet` | curated known-isoform-biology literature per resolved splicing-led gene (with Manubot citekeys) |
| S13a | `tableS13a_axis_orthogonality.csv` | `axis_orthogonality_summary.parquet` | per-analysis switch vs abundance separation for all 17 analyses (the one per-analysis table that needs no deconvolution, so it spans the five GTEx regions S13 cannot) |
| S13 | `tableS13_composition_adjustment.csv` | `composition_adjustment.parquet` + `composition_adjustment_gtex.parquet` | how much of the DTU-without-DGE layer survives cell-type adjustment, with the gene-level overlap (`n_overlap`) and newly detected genes (`n_new`), and the marker counts with their denominator (`n_marker_tests`) (backs Fig 2c,d and `figCompositionAgeCoupling` D) |
| S14 | `tableS14_longread_orthogonal_confirmation.csv` | `switch_orthogonal_confirm/anchored_gene_confirmation.parquet` | per-gene long-read confirmation of the anchored switch pairs (backs S-real-8) |
| S15 | `tableS15_isa_concordance.csv` | `isa_concordance/*/summary.json` | independent-caller (satuRn) DTU concordance across all 17 analyses (backs S-real-9) |
| S16 | `tableS16_module_genetic_anchoring.csv` | `module_genetic_anchoring_meta/` | per-module splicing anchoring vs a size-matched permutation null — **a table on purpose, not a figure** |
| S17 | `tableS17_rbp_eclip_binding_support.csv` | `rbp/rbp_binding_support.parquet` | per-RBP eCLIP binding capacity, switched vs constitutive exons (backs S-real-6B). **ENCODE eCLIP is HepG2/K562, not brain:** this is binding capacity at alternative-exon sequence, not neuronal occupancy of these regulons |
| S18 | `tableS18_scz_convergence.csv` | `scz_age_projection/convergence.parquet` | SCZ-risk convergence per module with candidate RBP regulators (backs Fig 4E) |
| S19 | `tableS19_qtl_anchoring_sensitivity.csv` | `qtl_anchoring_meta/` + `qtl_anchoring_meta/sensitivity/<arm>/` | splicing-specificity contrast under the primary arm and the three pre-specified sensitivities — constraint-adjusted (gnomAD LOEUF + missense z + log expression, both covariate sets on the identical constraint-complete subset), threshold-free continuous (rank-INT of -log10 pval_beta), and SuSiE credible-set dose (backs Fig 3) |
| S20a | `tableS20a_coloc_convergence_global.csv` | `module_coloc_convergence/global.parquet` | per (trait, source) coloc concentration vs a size-matched null, and anchored-module enrichment under BOTH denominators (all module genes vs the CLPP-tested pool) — backs S-real-10A/C |
| S20b | `tableS20b_coloc_convergence_per_module.csv` | `module_coloc_convergence/convergence.parquet` | per-module colocalizing gene and locus counts for all five traits, with the `testable` flag and the leave-one-locus-out worst case — backs S-real-10B |
| S21 | `tableS21_module_context_heldout_dtu.csv` | `dtu_added_value/region_summary.parquet` + giant-module arm | module context against held-out satuRn DTU evidence, per analysis, before and after the co-expression comparator (backs figDtuAddedValue; manuscript Table S15) |
| S22 | `tableS22_allelic_imbalance_regions.csv` | `ase_junction_switch/<region>/allelic_summary.json` + `module_cis_control_summary.json` | within-donor allelic test per region, with its homozygous-at-lead calibration, estimator and between-donor agreement, and the module cis-control correlation (backs Fig 5a / Fig S19; manuscript Table S16) |
| S23 | `tableS23_module_cis_control.csv` | `ase_junction_switch/module_cis_control.parquet` | per-module cis-controlled-gene rate against age-association rank (backs Fig S19e; manuscript Table S17) |
| S24 | `tableS24_clpp_isoform_events.csv` | `coloc/coloc_isoform_events_combined.parquet` | every CLPP isoform event and whether its junction maps into the tissue-matched switch pair — the 76 events / 30 genes (backs Fig 5b,c; manuscript Table S18) |
| S25 | `tableS25_signal_coloc_nominations.csv` | `coloc_signal_susie/all_introns/{genes,cells_hierarchy}.parquet` | signal-level sQTL nominations (PP4_sQTL >= 0.8), all-introns arm, with the headline estimator — 127 rows, 119 gene x trait cells, 116 genes (manuscript Table S19) |
| S26 | `tableS26_coloc_modality_contrast.csv` | `coloc_signal_susie/all_introns/contrast.parquet` | paired sQTL vs eQTL colocalization contrast, the unbiased check on the sQTL-preferential selection (manuscript Table S20) |
| S27 | `tableS27_ldsc_partitioned.csv` | `ldsc/ldsc_partitioned.parquet` | partitioned heritability of the switch-derived QTL annotations, single and joint models (backs Fig 5e, which reads this CSV; manuscript Table S21) |
| S34 | `tableS34_cohort_description.csv` | `inputs/bundles/*/*/{samples.parquet,manifest.json}` | per-analysis discovery cohort: donors, diagnosis, sex, age (GTEx top-coded at 70), ancestry, library kit, RIN, PMI or ischemic time, death classification, feature counts around the expression filter and MuSiC reference (manuscript Table S29; Methods) |
| S35 | `tableS35_synthetic_scenarios.csv` | `configs/synthetic_grid.yaml` | synthetic scenario grid: dimensions, swept and fixed parameters, seeds and datasets per scenario (manuscript Table S30; Methods) |

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
raw totals. **Read rates, not totals** — totals scale with module count, and IsoGraph runs
several-fold finer (per-region counts in S2). *(Requoted 2026-09-12 to the 2026-09-09
matched-baseline re-fit; the earlier 0.336 / 0.268 / 0.217 are stale.)* The phenotype-sig rate
is highest for the switch-only baseline (`wgcna_switch_only` **0.389**) and then `isograph`
(**0.274**), over the multiplex and abundance-fed ones (`wgcna_multiplex` 0.187, `wgcna_gene`
0.180): phenotype sensitivity comes from the switch features, not the inference method. GO
enrichment is abundance-dominated (`wgcna_gene` 0.885 ≫ `isograph` 0.235). On identical
multiplex features `isograph` (0.274) > `wgcna_multiplex` (0.187), but that gap is
granularity-confounded. **The comparison is feature-matched, not pipeline-matched:** IsoGraph
residualizes its discovery covariates inside the fit, the WGCNA baselines regress RIN only.
Supports: IsoGraph is not a better module-level enricher; its value is the switch layer.

## Table S2 — Three-baseline comparison, per cohort × region

The 66-row per-analysis ledger behind S1 (cohort, region, method, features, module count and
median size, GO-enriched / phenotype-sig / both counts and fractions). Use to confirm the
pooled rates are not driven by a single region and to read the per-region module-count
disparity that makes raw totals non-comparable.

## Table S3 — QTL splicing-specificity contrast (all-tissue meta)

Random-/fixed-effects meta-analysis of the **paired sQTL-OR / eQTL-OR ratio** within each
analysis, by module set and graph method. The ratio cancels the shared cis-QTL depletion of
constrained network genes (Table S5), isolating whether *splicing* genetics is spared.
*(Requoted 2026-09-12 from `qtl_anchoring_meta_contrast.parquet`; the earlier GO-invisible
1.13 and GO-visible "internal control" readings predate the 2026-08-29 refresh and are wrong.)*
For `isograph`, the fixed-effects ratio is > 1 in all-modules (**1.068**, p = 1.3e-5, k = 17)
and phenotype-significant modules (**1.111**, 95% CI 1.048–1.177, p = 3.6e-4, I² = 0.23,
k = 12). It does **not** localize to GO-invisible modules (1.068, p = 0.077), and GO-visible
modules are no longer a null (1.084, p = 0.050). Across the matched WGCNA baselines no set clears
0.05: `wgcna_switch_only` phenotype-significant 1.025 (p = 0.58), `wgcna_multiplex`
phenotype-significant **1.046 (p = 0.053, I² = 0.75)**, the closest. Columns: method, module
set, k analyses, fixed- and random-effects ratio + 95% CI, p, Q and I². Supports, **at set level
only**: splicing QTL are relatively spared in phenotype-associated IsoGraph modules. This is
supporting evidence, not the paper's headline; the per-gene tests (GTEx signal-level,
BrainSEQ in-sample coloc) lean toward expression.

## Table S4 — Matched-baseline QTL specificity contrast (8-tissue common set)

The same contrast restricted to the tissues where all three graph methods have a result, so
`isograph` vs `wgcna_switch_only` vs `wgcna_multiplex` is a like-for-like comparison on
identical tissues. `isograph` stays positive in phenotype-significant modules (**1.112**,
p = 6.5e-4, k = 7, I² = 0) and all modules (1.110, p = 2e-6, k = 8); GO-invisible does not
clear (1.066, p = 0.11). The baselines do not clear 0.05 in any set: `wgcna_switch_only`
phenotype-significant 1.015 (p = 0.74); `wgcna_multiplex` phenotype-significant **1.044
(p = 0.070, I² = 0.80)**. Since the 2026-09-09 baseline re-fit the multiplex baseline sits near
the threshold rather than at a flat null, so read the matched-baseline control as
"IsoGraph-only at 0.05, narrowly", not as a clean method effect. Supports: the set-level contrast
is not reproduced by WGCNA on the identical switch features.

## Table S5 — Raw cis-QTL enrichment ORs (eQTL & sQTL meta)

Per-method, per-QTL-kind meta ORs (foreground module genes vs background) by module set. Both
eQTL and sQTL ORs are < 1 everywhere (e.g. `isograph` all-modules eQTL 0.82, p=5e-122; sQTL
0.88, p=1.6e-25): network/module genes are cis-QTL-depleted, as expected for constrained
genes. This shared depletion is the baseline the S3/S4 **ratio** removes — the point of the
specificity contrast is that sQTL depletion is *weaker* than eQTL depletion in the disease
modules, not that either is enriched. Supports: the raw-OR context that prevents
mis-reading the specificity ratio as enrichment.

## Table S6 — GO-invisible SCZD disease switch modules (BrainSEQ caudate)

Per-module ledger for the **eight** SCZD-associated IsoGraph switch modules (`pheno_fdr ≤ 0.1`),
plus a pooled `_background` row over all switch transcripts. Six return zero GO terms
(`go_invisible = True`) and two are GO-visible (M012, M011). Four of the six carry a real
anticorrelated transcript pair in nearly every member (M026 23/23, M022 29/29, M023 27/27, M010
65/83); M025 is partial (11/24) and M020 weak (4/30), with max switch strength 1.02–1.72 and
7–509 significant switch transcripts across the six. *(Corrected 2026-09-12: an earlier
version of this table listed four modules from a stale gate run.)* Driver functional-consequence fractions (CDS / coding-status / biotype /
UTR change) sit at or above the pooled background (0.84 / 0.67 / 0.74 / 0.61), so
GO-invisibility reflects GO's gene-level/abundance bias, not low module quality. Supports: the
disease switch signal is genuine, functionally consequential isoform regulation invisible to
pathway enrichment (biology gate PASS, complementary form).

## Table S7 — Per-region module trust funnel

One row per cohort × region × method. Q1 stability (`n_modules`, `n_trusted`,
`frac_trusted` above a size-matched permutation null at FDR<0.05, with the median
co-assignment density and its null); Q2 drivers (`median_driver_rho`, `frac_positive_rho`
of the shared-gene switch-loading Spearman ρ across split halves, IsoGraph only — WGCNA has
no transcript drivers); within-cohort split-half aging concordance
(`n_sign_concordant`/`n_both_age_sig` of `n_split_half_pairs`); and Q4 complementarity
(`median_frac_dtu_without_dge`, `median_frac_in_wgcna_age`, IsoGraph only). Across the six
regions: **93/118** IsoGraph modules trusted against **62/72** for the matched WGCNA
baseline — IsoGraph is *not* the higher trusted fraction, it is trustworthy at a much finer
granularity; driver ρ medians **0.72–0.82** with 94–100% of pairs positive; split-half age
direction agrees in **22/23** IsoGraph and **14/14** WGCNA both-significant pairs.
Cross-cohort counts are in S7b (transfer) and S7c (membership matching), not here, because
they are not per-region quantities. Supports: the per-module trust scaffold the biological
claims rest on.

*Regenerated 2026-09-28 from the Leiden-2.0 production fits. The earlier 250/266, 0.66–0.88
and 23/130 in this legend were computed at Leiden 5.0 and are stale; so is
`MODULE_TRUST_SUMMARY.md`, which has not been rebuilt.*

## Table S7a — Split-half module ledger

One row per module per method (190 rows): size, co-assignment density, the size-matched
permutation null mean, permutation *p*, BH *q* and the trusted flag, plus
`best_match_jaccard`. The Jaccard column is descriptive only — it is granularity-confounded,
which is why the trust criterion is the chance-calibrated co-assignment test and not this
column. Supports: the trusted counts in Fig 3a,b are auditable module by module.

## Table S7b — Cross-cohort eigengene projection ledger

One row per trusted module per direction (155 rows): the frozen-weight signed kME against
its type-matched null with BH *q*, the module's age correlation in each cohort, that
correlation expressed as a *z* against same-weight random projections, and the sign match.
The `raw_sign_match` and `raw_both_sig` columns are the point of the table: on the raw
projected correlations the three BrainSEQ→GTEx pairs give 44/44, 1/50 and 32/36 sign
agreement — a property of each target cohort's own age-correlated structure, not of the
modules. Supports: the transfer claim is reported on the standardised statistic, and the
raw one is shown rather than omitted.

## Table S7c — Cross-cohort matched-pair count against its permutation nulls

The full grid: 2 methods × 3 statistics (covariate-free Pearson, covariate-adjusted linear,
covariate-adjusted spline F) × 3 covariate modes × 2 nulls (`age`, Freedman–Lane on the age
association; `matching`, permuting which target module each source module is matched to —
the stricter). `is_published_statistic` flags the covariate-free Pearson arm the Results
text quotes. **Read the two arms together:** against the matching null the Pearson arm gives
1/38 (*p* = 0.91) for IsoGraph and 2/52 (*p* = 0.15) for WGCNA, while the stage's
pre-registered covariate-adjusted spline arm gives 0/38 (*p* = 1) and 5/52 (*p* = 0.16).
Neither clears. Per the decision rule in `REPLICATION_PERMUTATION.md` the word *replication*
must not be used for this arm; the honest wording is "matched modules with concordant age
effects". Supports: the membership-matching null result, reported in full rather than at its
most favourable setting.

## Table S7d — Functional preservation of matched cross-cohort pairs

Matched BrainSEQ↔GTEx module pairs against a null that re-pairs each module with a random
target module **from the same gene-count decile** (all of these similarities grow with
module size). Median gene Jaccard is 0.038 — the overlap these measures are asked to look
past. IsoGraph clears on both computable measures (GO-term Jaccard 0.097 vs 0.030,
*p* = 0.001; cell-type profile *r* = 0.208 vs 0.040, *p* = 0.006), while the WGCNA baseline
sits on its own null (GO Jaccard 0.136 vs 0.139, *p* = 1). Rows with `n_finite = 0` are
**untested, not negative**: transcript structure had no computable pair, and cell-type
profiles are produced only in the IsoGraph artifact tree. Supports: the cohorts recover
related biological programs at a higher level of organisation than gene identity — the
positive half of a paragraph whose first half is a null result.

## Table S7e — Split-half agreement across the Leiden resolution sweep

Mean adjusted Rand index and normalised mutual information between the two halves of each
donor split, over five seeds, at eight Leiden resolutions (0.5, 1, 2, 3, 5, 8, 12, 20).
Production is 2.0 (PI decision, 2026-09-16) and is flagged in `is_production`; the WGCNA
rows carry no resolution, because there is none to set, and are kept so the sweep is read
against the fixed baseline rather than against itself. Mean ARI across regions is 0.37–0.47
over the whole sweep (production 0.42), and the per-region best resolution is 0.5, 0.5, 1, 3, 20
and 20 — that is, no resolution is uniformly best and nothing singles out 2.0 as
tuned-to-fit. Per-region production values span ARI 0.19–0.57, against 0.29–0.71 for WGCNA
on the same splits, which is the granularity contrast the funnel reports throughout: WGCNA
agrees with itself more readily on a partition roughly 3.6× coarser. Supports: the
stability claim is a property of the representation, not of the resolution chosen for it.

## Table S7f — Raw against null-standardised projected-age sign agreement

One row per method × direction × region, with the switch-axis sign convention for the
IsoGraph rows. The raw and standardised arms are counted over different denominators and
both are written: `raw_sign_match_both_sig` over modules whose projected age correlation is
significant in both cohorts, `std_sign_match_all` over every age-testable module. The raw
arm is why the reported statistic is standardised — among raw both-significant modules the
three BrainSEQ→GTEx regions give 8/8, 8/8 and **0/13**: two regions agree unanimously and
one disagrees unanimously, which is a property of each target cohort's own age-correlated
structure (RNA quality, ischemic time, composition), not of the modules. About 45% of
orientable genes have their switch axis fitted with the opposite sign in the two cohorts
(2,505–2,616 of 5,450–5,739 per region; median |cosine| 0.50 over all shared genes, 0.80
among orientable ones), so a raw projected coefficient carries a sign the projection did
not earn. Supports: the standardisation in the projection result is a stated method
choice with its own evidence, not a post-hoc rescue.

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

## Table S13a — Switch against abundance axis separation, per analysis

One row per analysis (17): genes tested, the median absolute correlation between the
gene-abundance and isoform-switch coordinates with its quartiles, and the fraction of genes
below 0.1, below 0.3 and above 0.5. Median |*r*| runs 0.111 (BrainSEQ aging hippocampus) to
0.243 (GTEx nucleus accumbens), with 23–46% of genes under 0.1 and 1.1–16.8% over 0.5, over
12,042–13,222 genes per analysis. The quartiles and tail fractions are kept because the
claim is about the bulk of the distribution rather than a central value. Unlike S13, this
test needs no deconvolution, so it covers all 17 analyses including the five GTEx regions
with no matched snRNA reference. Supports: the switch coordinate is not a re-description of
abundance in any analysis, while being far from orthogonal in some.

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
