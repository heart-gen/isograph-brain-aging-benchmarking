# Cell Genomics impact additions — progress tracker

Working doc (uncommitted). Tracks the three high-impact biological additions from the
approved plan (`~/.claude/plans/review-this-repo-for-serialized-nygaard.md`). **Paused after
Item 1 per user.** Items 2 & 3 have code built and are ready to run.

Provenance: all three extend committed CLIs under `isograph_benchmark/real_data/` + reuse
on-disk artifacts. No `git add` done yet. Interpreters/SLURM conventions per `AGENTS.md`.

---

## Item 1 — Cell-type composition adjustment — ✅ DONE (result is consequential)

**What was built**
- `isograph_benchmark/real_data/celltype_composition.py` — subcommands `fractions` (BrNum→
  sample_id join of the committed MuSiC proportions from
  `../sex_context_brain/cell_proportion_estimate/_m/`; marker-depletion cut) and `meta`
  (with-vs-without rollup → `04_module_characterization/_m/COMPOSITION_ADJUSTMENT_SUMMARY.md`).
- `incremental_association.py` — added `--composition` flag: adds cell-type fractions as
  inference covariates (reference type dropped for the simplex), writes to
  `incremental_association_composition/` so the canonical baseline is preserved.
- SLURM wrapper `04_module_characterization/_h/10.celltype_composition_brainseq.sh` (array 1–4). **Ran clean
  as job 42833041.**

**Result (gene-level composition-unique = DTU-without-DGE signal, base → composition-adjusted)**

| region | base | adjusted | retained |
| --- | --- | --- | --- |
| SCZD caudate | 34 | **2** | 0.06 |
| aging caudate | 43 | **17** | 0.40 |
| aging DLPFC | 8 | **15** | 1.88 (adjustment sharpens) |
| aging hippocampus | 0 | 0 | — (already null) |

**Interpretation / manuscript impact**
- **SCZD switch signal is heavily composition-confounded (34→2).** This directly affects the
  north-star "34 SCZD genes" claim (`AGENTS.md`, abstract). Must be disclosed and the
  composition-adjusted number reported. **Confounder-vs-mediator caveat is central:** if
  disease *causes* the composition shift, adjustment over-corrects (removes real signal); if
  composition varies technically, adjustment is the correct control. Report both-ways and lean
  on layers that don't depend on it (module-level anchoring, GO-invisible gate).
- **A genuine composition-robust aging DTU layer survives** (caudate 17, DLPFC 15) — the
  aging story is more robust to composition than the SCZD disease story.
- Marker-depletion cut ~null (modules are not bags of cell-type markers).
- **Suggested reframing:** lead the DTU-without-DGE claim with the composition-robust *aging*
  layer; present SCZD with the adjusted count + the mediator caveat.

**Optional follow-on (Item 4):** GTEx MuSiC re-run for the aging replication arm — ✅ DONE
(partial replication: limbic/striatal robust, cortical collapses). See "Item 4" below.

---

## Item 2 — RBP binding evidence (ENCODE eCLIP) — ✅ DONE (within-gene contrast; a binding-*capacity* result, not evidence-backed neuronal occupancy)

**Ran:** `rbp_binding.py run` (job 42834007). First pass gave a saturated "100% supported"
(peak-dense RBPs blanket the transcriptome), so hardened to a **within-gene contrast**: binding
at each gene's SWITCHED exon vs its CONSTITUTIVE (shared) exon, exact McNemar over discordant
genes, BH-FDR. Neutralizes peak abundance (U2AF2/HNRNPU correctly read null).
Outputs `07_rbp_regulation/_m/rbp/{rbp_binding_calls,rbp_binding_support,rbp_binding_regulon}.parquet`
+ `RBP_BINDING_SUMMARY.md`.

**B1 fix (statistical, 2026-08-05).** The earlier headline pooled BH-FDR across all 192
`region×module×RBP` regulon rows — but 98.8% of gene×RBP binding facts are region-invariant, so
that family counts the *same* signal once per region (MATR3 ×5, QKI ×4 = one signal replicated),
inflating the "supported" count. Corrected to an **independence-respecting unit**: one two-sided
exact McNemar **per RBP** over its *unique nominated regulon genes* (deduped across regions,
bound if bound in ANY region's switch definition), BH across the testable RBPs (38 in the
canonical GC-matched arm, 39 in the deprecated flat arm). The per-regulon rows are retained
**descriptively** (`rbp_binding_regulon.parquet`, 193 rows GC / 192 flat), not as an FDR family.

**Result (per-RBP headline, GC-matched background, 38 testable RBPs):**
- **31/38 directionally preferential**; **25/38 binding-supported** (preferential AND BH q≤0.05),
  spanning **55 of 66** independent motif families.
- **Median switched−constitutive gap = 0.024** — small and near-universal. This is binding
  **capacity** at alternative vs constitutive exons, **not** factor-specific occupancy of the
  predicted regulons.
- Largest gaps: **PCBP2** (0.116), **MATR3** (0.105, q≈0), HNRNPC (0.078), HNRNPM (0.074),
  HNRNPU (0.062), TIAL1 (0.060), ELAVL1 (0.054), PUM2 (0.045), GRSF1 (0.043), QKI (0.040).

> **Which background.** The numbers above are the **GC-matched** arm
> (`RBP_BINDING_SUMMARY.md`, `rbp_binding_support.parquet`) — the canonical one. The
> **flat-background** arm (`*_flatbg`) gives 31/39 preferential, 17/39 supported, median gap
> 0.013, and is deprecated: it inflates AU-rich binders (ELAVL, CPEB, hnRNPD). This document
> quoted the flat arm as the headline until 2026-09-09; the same deprecated arm produced the
> "829 hits / 129 RBPs" front-page error corrected in `41e293f`. Quote the GC arm.

**Interpretation (corrected — do NOT call this "evidence-backed").** A subset of nominated factors
(led by MATR3) bind their regulon genes' switched exons at a modestly higher rate than the
constitutive exons — this confirms binding *capacity* at alternative-exon sequence, a supporting
detail for the Fig-4 mechanism sentence. It does **not** demonstrate that these factors selectively
occupy the switched sequence beyond a generic alternative-vs-constitutive skew, and it does **not**
support a blanket "GO-invisible modules are RBP regulons" claim. Caveats compound: ENCODE eCLIP is
HepG2/K562, not brain, so this is capacity, not neuronal occupancy. The neuronal-CLIP validation
program (Tier 1/2/3: PTBP2/TARDBP/NOVA2 cTag-CLIP + perturbation) that would have shown *brain*
occupancy is an **honest null** — cross-species liftover attrition leaves every context
underpowered (`descriptive_underpowered`/`not_estimable_sparse_discordance`), reported as a
characterized limitation, not a positive claim.

### (build notes)

**Context:** POSTAR3's host is dead (serves a GitHub 404); pivoted to **ENCODE eCLIP**
(cell-type caveat: HepG2/K562, not brain). Of 130 significant regulon RBPs, **39 have eCLIP**
— including the headliners KHDRBS1(SAM68), U2AF2, TARDBP(TDP-43), FUS, RBFOX2, HNRNPC/K,
ELAVL1, QKI, MBNL1, PUM1/2, IGF2BP1/2/3. Honest partial coverage.

**What was built**
- `isograph_benchmark/real_data/rbp_binding.py` — `fetch` (download GRCh38 eCLIP peak beds
  per RBP via the ENCODE experiment API → `inputs/raw/rbp_binding/<RBP>.bed.gz`) and `run`
  (overlap each gene's **switched exon interval** = symmetric-difference exonic space of the
  switch-pair transcripts, from the GENCODE v47 gtf cache, with peaks; join the binding flag
  into `rbp_regulon`; write `rbp_binding_calls.parquet`, `rbp_binding_regulon.parquet`,
  `RBP_BINDING_SUMMARY.md`).

**State**
- `fetch` **complete**: all **39/39** testable RBPs downloaded (130 MB) →
  `inputs/raw/rbp_binding/*.bed.gz` + `eclip_manifest.tsv`. (74k–1.96M peaks/RBP.)

**To resume**
1. `python -m isograph_benchmark.real_data.rbp_binding run` (moderate compute — SLURM: add a
   wrapper `07_rbp_regulation/_h/05.rbp_binding.sh`; interval overlap over 17 regions).
3. Read `07_rbp_regulation/_m/rbp/RBP_BINDING_SUMMARY.md`: the per-RBP binding-supported count
   (25/38, GC-matched background) and the small median switched−constitutive gap (0.024).
   Do **not** read `RBP_BINDING_SUMMARY_flatbg.md` for the headline. Supports the Fig 4 RBP
   sentence with binding *capacity* for a subset of factors — **not** an "evidence-backed"
   neuronal-occupancy upgrade (see B1 fix and the neuronal-CLIP honest null above).

---

## Item 3 — Module-level genetic anchoring — ✅ DONE (nuanced result)

**Ran:** `module_genetic_anchoring.py` array (job 42833717, 17 analyses, n_perm=1000) + `--meta`.
Outputs `05_genetic_anchoring/_m/module_genetic_anchoring_meta/` (`MODULE_GENETIC_ANCHORING_META.md`,
`_all.parquet`, `_meta.parquet`) + per-analysis `MODULE_GENETIC_ANCHORING.md`.

**Result (177 pheno-sig modules; per-module log sQTL-OR − log eQTL-OR vs size-matched perm null):**
- Median contrast ≈ 0; **14/177 (~8%) anchored at perm_p≤0.05** (vs ~9 expected) — modest excess.
- Real tail of anchored programs: sQTL OR 1.5–2.9 vs eQTL OR 0.3–0.9, several in neurodegen
  regions (SN M031/M036, hippocampus M026, amygdala M009; top ACC M011 p=0.001, HIP M026 p=0.001).
- **NOT concentrated in GO-invisible modules** per-module (GO-inv 7/124≈5.6%≈chance; GO-vis 7/53=13%).

**Interpretation (reviewer-relevant):** the perm null draws from the co-switch universe, so the
test asks "is *this module* more splicing-anchored than the average co-switch set." Answer: the
**GO-invisible-specific splicing anchoring is a POOLED-gene property (Finding 3), not a per-module
one.** Use this test to name **exemplar anchored programs** (the 14, neurodegen-weighted), NOT to
claim GO-invisible modules are anchored as units. Deeper follow-ons (eigenswitch×TOPMed eigen-QTL;
module-restricted S-LDSC) remain the way to test the *coordination* directly — not built.

### (superseded design note)

**What was built**
- `isograph_benchmark/real_data/module_genetic_anchoring.py` — per phenotype-associated
  module, the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the
  module's own member genes vs a **size-matched permutation null**, reusing
  `qtl_anchoring.matched_enrichment` + the GTEx sGenes/eGenes catalogs. No controlled
  genotypes. Writes `module_genetic_anchoring.parquet` + `MODULE_GENETIC_ANCHORING.md`.

**Rationale:** defends the "genetically anchored *programs*" title claim — the pooled QTL
result doesn't show any individual module is anchored; this tests each module as a unit.

**To resume**
1. `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis brainseq-sczd`
   (then `brainseq-aging`, `gtex-aging`). Permutation loop (default n_perm=1000 × 2 logits/
   module) is the heavy part → SLURM wrapper `05_genetic_anchoring/_h/04.module_anchoring.sh`.
2. Add a cross-analysis meta rollup (mirror `qtl_anchoring_meta`) if ≥1 module is anchored.
3. **Deeper follow-ons (need new controlled-genotype extraction, not built):** eigenswitch ×
   TOPMed genotype eigen-QTL (reuse `scz_age_projection.load_genotypes` + plink2 in eqtl env);
   module-restricted S-LDSC annotation (reuse `ldsc_annot_prep.py`).

---

## Item 4 — GTEx composition (aging replication arm) — ✅ DONE (partial replication)

**Why:** Item 1 showed the composition-robust aging DTU layer survives in BrainSEQ (caudate 17,
DLPFC 15). This replicates that de-confounding in the independent GTEx cohort, so the lead
"composition-robust aging" claim is not BrainSEQ-only. GTEx had no prior deconvolution, so this
runs MuSiC in-repo against the **same Tran/LIBD snRNA references** (seed 13) used for BrainSEQ.

**What was built (all reproducible, no `git add` yet)**
- `celltype_composition.py` extended: `export-gtex` subcommand (genes×samples raw-count bulk
  parquet per region), `GTEX_REF` region→reference map, sample_id-keyed `build_fractions`
  branch, `gtex-aging` in `fractions`, and a GTEx replication section in `meta`
  (`02_module_discovery/gtex/_m/composition/GTEX_COMPOSITION_SUMMARY.md` + a section appended to the
  shared `COMPOSITION_ADJUSTMENT_SUMMARY.md`).
- `04_module_characterization/_h/gtex_music_deconv.R` — self-contained MuSiC deconvolution (inlines the
  board-level cell-type mapping incl. striatal MSN D1/D2; reads the exported bulk parquet via
  arrow; writes `music-proportions-gtex-<region>.tsv` + `marker_stats_genes.gtex-<region>.csv`).
- `04_module_characterization/_h/11.gtex_composition.sh` — array 1–8: export → R MuSiC → fractions →
  `incremental_association --composition`, per region.

**Region coverage (honest):** only the **8 GTEx regions with a defensibly matched Tran
reference** are deconvolved — striatum→NAc (caudate/putamen/NAc basal ganglia), cortex→DLPFC
(cortex, frontal_cortex_ba9), and AMY/sACC/HPC direct (amygdala, ACC, hippocampus). The 5
without a matched panel (cerebellum, cerebellar_hemisphere, hypothalamus, spinal_cord,
substantia_nigra) are **not** deconvolved rather than forced against a mismatched reference.

**Result — region-dependent, honestly PARTIAL replication (8 regions, base → composition-adjusted comp_unique)**

| GTEx region | base | adj | retained | note |
| --- | --- | --- | --- | --- |
| amygdala | 7 | 6 | 0.86 | robust |
| hippocampus | 61 | 21 | 0.34 | robust |
| caudate_basal_ganglia | 32 | 10 | 0.31 | robust |
| putamen_basal_ganglia | 4 | 3 | 0.75 | robust |
| nucleus_accumbens_basal_ganglia | 1 | 4 | — | robust |
| anterior_cingulate_cortex_ba24 | 545 | 88 | 0.16 | retains (88) |
| **frontal_cortex_ba9** | 531 | **0** | 0.0 | **collapses** |
| **cortex** | 438 | **0** | 0.0 | **collapses** |

**6/8 retain ≥1 composition-robust gene; the 2 cortical regions collapse to 0.**

**Interpretation (do NOT oversell):** the **limbic/striatal** aging switch layer reproduces as
composition-robust across cohorts — a real independent-cohort win. But the two **cortical**
regions (frontal_cortex_ba9, cortex) collapse entirely, *even though BrainSEQ DLPFC survived*.
Cortical composition is the most age-coupled and these are the weakest reference matches (GTEx
cortex/BA9 → DLPFC snRNA), so the collapse is consistent with genuine composition-confounding
of cortical aging switches **and/or over-adjustment** (7–8 age-correlated fraction covariates
absorbing the age-spline variance) — the data can't cleanly separate them. **Manuscript use:**
cite GTEx as a *partial* composition replication (limbic/striatal robust; cortical
composition-entangled), a stated limitation, not a claim of universal robustness. This is
reviewer-defensive: the test is clearly not rigged to always "survive."

**Validated locally first (amygdala):** MuSiC 181/181 samples, 8 cell types summing to 1;
7→6 (Oligo dropped as simplex ref). Array 2–8 (job 42847015) completed clean; `meta` regenerated
`GTEX_COMPOSITION_SUMMARY.md` + the GTEx section of the shared `COMPOSITION_ADJUSTMENT_SUMMARY.md`.

---

## Resume checklist
- [x] Item 2: eCLIP fetch + overlap done; B1 fixed (per-RBP independence-respecting unit,
  25/38 binding-supported at a GC-matched background, median gap 0.024 = binding capacity);
  neuronal-CLIP tiers are an honest null.
- [x] Item 3: run module anchoring (SLURM) → meta → decide on eigen-QTL / S-LDSC depth.
      DONE — see §Item 3 above for the nuanced result; the eigen-switch QTL was rejected on
      power and module-level coloc convergence came back null in 10/10 (trait, cohort) cells.
- [ ] Fold Item 1 result into the manuscript (adjusted SCZD count + mediator caveat;
      composition-robust aging layer as the lead DTU claim).
- [ ] Decide commit/PR scope (nothing staged yet; never `git add -A`).
