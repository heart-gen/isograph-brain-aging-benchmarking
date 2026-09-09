# Manuscript Evidence, Requirements, and Strategy — IsoGraph Brain-Aging

Provenance-labeled manuscript-planning document produced by the repository-to-manuscript
skill. Three layers are kept visibly separate: **[EVIDENCE]** = what the repository
demonstrates, **[REQUIREMENT]** = what an external authority requires, **[RECOMMENDATION]**
= strategist judgment, **[UNRESOLVED]** = conflicts/gaps that could change the plan.

> Scope note: this repository is already at an advanced manuscript stage — figures are
> ordered (`manuscript/FIGURE_ORDERING.md`), a Results section is drafted
> (`manuscript/GENETIC_ANCHORING_RESULTS.md`), and a separate Manubot repo holds
> `content/*.md`. This plan therefore audits and *re-anchors* an existing plan rather than
> inventing one, and flags the two things still genuinely open: the Cell Genomics re-target
> mechanics. (The SCZ age-projection convergence layer, in flight when this plan was first drafted,
> is complete but the convergence result is **RETRACTED 2026-08-30** — the hypergeometric used the
> wrong background; against the coloc-tested pool it is null (P=0.40). See Results 6 and
> `module_coloc_convergence.py`. The module-disruption layers are unaffected.)

---

## Part I. Repository Evidence

### 1. Repository and Analysis Inventory

| Analysis | Question | Data | Main finding | Evidence strength | Validation | Reproducibility | Sources |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Synthetic benchmark | Does IsoGraph recover switch modules where truth is known? | Simulated paired grid | Module recovery favours IsoGraph in 7/15 scenarios (13,410 runs, 16 scenarios) — wins where switching dominates, loses when abundance-dominated or degraded; complete switch-gene detection | Strong (ground truth) | Internal (synthetic truth) | CLI + SLURM, seeds | `01_synthetic_benchmark/03_metrics/`, `fig1`, `tableS_benchmark_summary.csv` (supp) |
| Confound ablation | Does residualization buy robustness? | Synthetic confounds | `isograph_vae_residual` holds under composition/batch/depth; degradation is the honest limit; abundance fallback recovers it | Strong (synthetic) | Internal | `figS8/figS9` | AGENTS.md §3, `01_synthetic_benchmark/03_metrics/figures/` |
| Module trust funnel | Are real-data modules reproducible? | BrainSEQ + GTEx fits | 236/266 modules chance-trusted (6 regions); driver ρ 0.77–0.82; 25 aging replications vs ~6 abundance baseline | Strong | Split-half + perm null + cross-cohort | `03_module_trust/`, `figTrustFunnel`, `tableS7` | 
| GO-invisible gate | Are disease switch modules real DTU, not abundance/low-quality? | BrainSEQ caudate SCZD | 4/4 pheno-sig SCZD modules GO-invisible (M026/M020/M010/M023), all carry real anticorrelated switch pairs, consequence ≥ background | Strong (within cohort) | Within-gene perm null; internal control vs background | CLI + SLURM, seeds, post-refit 2026-06-29 | `GO_INVISIBLE_GATE_SUMMARY.md`, `go_invisible_gate.parquet` |
| QTL splicing-specificity | Are co-switch modules anchored to *splicing* genetics, method-specifically? | GTEx v11 sQTL/eQTL | sQTL/eQTL ratio **1.111 pheno-sig (p=3.6e-4)**, all 1.068 (p=1.3e-5); GO-invisible 1.068 (p=0.077) and GO-visible 1.084 (p=0.050) indistinguishable — **no GO-invisible localisation**; IsoGraph-only vs matched WGCNA baselines | Strong (matched-baseline null) | Matched-baseline null (primary control); IVW+DL meta; 2026-08-29 refresh — never quote earlier numbers | CLI + SLURM, 17 analyses | `QTL_ANCHORING_SUMMARY.md`, `qtl_anchoring_meta_contrast*.parquet` |
| Colocalization / deep-dive | Do disease variants resolve to specific isoform switches? | GTEx v11 SuSiE × 5 GWAS (SCZ/AD/PD/LBD/ALS), eCAVIAR | 141 events / 68 genes; 12 splicing-led, **all GO-invisible**; 9 GTEx-concordant; SNCA cross-disease (LBD+PD) same alt-first-exon | Moderate (modest posteriors; set-level coherence) | Set-level pattern; cross-disease concordance; literature layer | Deterministic per-gene join CLI | `GENETIC_ANCHORING_RESULTS.md`, `deep_dive/*.md` |
| Partitioned heritability (S-LDSC) | Does switch layer carry disease heritability? | baselineLD v2.2 × 5 traits | Aging sQTL annotation enriched all neurodeg traits (LBD 7.17×, PD 3.80×, AD 3.37×, ALS 2.96×, SCZ 1.74×); disease-SCZ is eQTL-led | Moderate–strong | Single- + joint-annotation models | `05_genetic_anchoring/_m/ldsc/` | `LDSC_SUMMARY.md` |
| GWAS/MAGMA resolution | Is module-GWAS enrichment a size artifact? | MAGMA per module | At res 5.0 giant-module artifact disappears (0/8 sig giant vs 18/43 res2.0, 79/99 gene-level WGCNA); SCZ signal survives 6 regions | Strong (artifact control) | Resolution sweep + WGCNA comparison | CLI + SLURM | `GWAS_RESOLUTION_SUMMARY.md`, `figGwasResolution` |
| Switch consequence | What do the switches do structurally? | Within-gene perm null | Only productive UTR (1.27×, 10/10) + CDS (1.04×, 10/10) enriched; degradative classes never; identical in GO-invisible | Strong | Within-gene perm null; Fisher-combined | CLI + SLURM | `SWITCH_CONSEQUENCE_SUMMARY.md` |
| RBP regulons | Do modules share candidate trans-regulators? | ATtRACT PWMs + MOODS | 714 module×RBP sig (hypergeometric q<0.05), 138 RBPs over 14 regions; 310 under the opportunity-adjusted GLM with only 43 in common; recurrent neuronal 3′UTR/splicing factors (KHDRBS1/ELAVL4/PPRC1 8/14…) | Moderate (motif prediction, not binding) | Per-module over-rep test; recurrence across regions | CLI + SLURM (mature/intronic/combined) | `RBP_REGULON_SUMMARY.md` |
| Clinical consequence | Is the switch layer constrained but non-coding? | gnomAD LOEUF + ClinVar | Switch genes more LoF-constrained (LOEUF 0.72 vs 0.94); switched exons LOWER ClinVar P/LP (ratio 0.18) | Strong | Fisher-combined MWU; CDS-only robustness | CLI + SLURM | `CLINICAL_CONSEQUENCE_META.md` |
| Three-baseline comparison | Is IsoGraph globally superior to WGCNA? | 17 analyses × 4 methods | **No.** Phenotype rate lives in switch features (both switch-fed win); GO abundance-dominated; one clean effect: isograph > wgcna_multiplex on identical features | Strong (scope-bounding) | Per-module rates; matched features | Login aggregation CLI | `BASELINE_COMPARISON_SUMMARY.md` |
| Abundance/switch separability | Is the switch channel non-redundant with abundance? | Incremental association | Axes ~orthogonal (median \|r\|≈0.13); 34 SCZD / 43 caudate composition-unique genes; NREP flat abundance (p=0.93) but switch p=7e-5 | Strong | De-confounded incremental test | CLI + SLURM | `figSeparation`, `abundance_structure_separation.py` |
| **SCZ age-projection (NEW)** | Do SCZ-risk loci converge on age-sensitive switch programs disrupted in disease? | BrainSEQ caudate_sczd + TOPMed genotypes | **[RETRACTED 2026-08-30]** The convergence hypergeometric used ALL module genes as the background. A gene can only colocalize if it was coloc-TESTED (under a GWAS peak with a QTL credible set), and anchored modules are *defined* by MAGMA SCZ enrichment, so their genes enter the tested pool preferentially. Against the tested pool the background is already 45% and the same 15/31 = 48% is **null (P=0.40)**. `module_coloc_convergence.py` repeats this for all five traits with a size-matched permutation null: **no concentration anywhere (P=0.19-1.00), no module survives FDR (min q=0.23)**. The convergence claim is withdrawn; the module-disruption layers (B 4/10, D 3/10, C 10/10) and the named candidate regulators are unaffected and stand on their own. Disruption layers: disruption B 4/10, D 3/10, C 10/10; named regulators (M002→SNRNP70/ZCRB1, M008→ZC3H10/RBM14/CELF5, M006→DDX58/ADAR/YTHDC1); single-locus genotype concordance **null** 28/62 P=0.81 → supplementary | Established (module-level) | Convergence enrich + preservation + age-deviation + direction concordance; genotype layer null/underpowered (n~62) | CLI + SLURM (18), seeds, 2000 perm | `scz_age_projection.py`, `SCZ_AGE_PROJECTION.md`, `convergence.parquet`, `mechanism_rbp.parquet` |

### 2. Established Findings

#### Finding 1: IsoGraph is a complementary DTU-without-DGE layer, not a globally superior method
**Evidence status:** Established.
**Supporting analyses:** three-baseline comparison; de-confounded gene-level test; abundance/switch separability.
**Key quantitative evidence:** pooled per-module phenotype-sig rate wgcna_switch_only 0.336 > isograph 0.268 > wgcna_multiplex 0.189 ≈ wgcna_gene 0.180; GO rate abundance-dominated (wgcna_gene 0.885 ≫ isograph 0.217); one clean method effect isograph 0.268 > wgcna_multiplex 0.189 on identical features; 34 SCZD / 43 caudate composition-unique genes.
**Validation:** matched-feature baselines; per-module rates (not totals).
**Limitations:** IsoGraph does not win enrichment/phenotype rates outright — this finding *bounds* the paper's claim.
**Repository sources:** `BASELINE_COMPARISON_SUMMARY.md`, `figSeparation`, AGENTS.md north-star.

#### Finding 2: The phenotype-associated disease switch modules are real, GO-invisible isoform switching
**Evidence status:** Established (single cohort).
**Key quantitative evidence:** 4/4 SCZD pheno-sig modules GO-invisible; 21–73 members/module carry real anticorrelated pairs; consequence fractions ≥ pooled background (CDS 0.84 bg; M026 1.00).
**Validation:** within-gene permutation null; internal control vs background switch population.
**Limitations:** one cohort/region; drivers heterogeneous within a module (shared switch axis, not shared GO process).
**Repository sources:** `GO_INVISIBLE_GATE_SUMMARY.md`.

#### Finding 3: The phenotype-associated switch layer is genetically anchored to splicing, method-specifically
**Evidence status:** Established (strongest orthogonal result).
**Key quantitative evidence:** sQTL/eQTL specificity contrast **1.111 pheno-sig (95% CI 1.048–1.177, p=3.6e-4, I²=0.23)**, all_modules 1.068 (p=1.3e-5). IsoGraph-only on the 8-tissue common subset: pheno_sig 1.108 (p=0.001, I²=0.00) and all_modules 1.105 (p=5.8e-6) vs matched WGCNA baselines null (0.968 / 0.980, p≥0.41).
**SNCA is orthogonally confirmed; CTSH is not (2026-09-03).** BrainSEQ short-read junction
PSI on the exact Fig 4A contrast gives minor-form usage 0.189 (DLPFC, n=222) and 0.234
(caudate, n=238) against a pre-registered 0.05 threshold, versus 0.29% in ONT long-read —
so the long-read failure was an assay limitation and SNCA may carry a main figure. CTSH
reaches only 0.016–0.020 in hippocampus and stays off any main figure, per the same
pre-registered rule. CLI `junction_coloc_confirm.py`, wrapper `06_switch_mechanism/_h/12`.

**Both negative controls are now written into `GENETIC_ANCHORING_RESULTS.md`** (section
"Two negative controls bound the genetic claim", added 2026-09-03): the per-gene modality
contrast is null in all four gene pools, and module-level coloc convergence is null in 10/10
(trait, cohort) cells. Fig 3's claim is set-level and survives both; neither may be quoted as
corroborating it.

**Scope bound — the effect does not localise to the GO-invisible modules.** GO-invisible (1.068, p=0.077) and GO-visible (1.084, p=0.050) are indistinguishable on the 2026-08-29 refresh, so this finding is about the *phenotype-associated* set, not the DTU-without-DGE subset. The GO-invisible claim rests on the content gate (Finding 2 / S-real-3) instead. **Never quote a pre-2026-08-29 anchoring number** — the retired 1.172 / I²=0.00 came from a stale `module_enrichment` join that relabelled the GO partition.
**Validation:** the matched-baseline null is the primary internal control (identical switch features, inference varied); IVW-FE + DL-RE meta over 17 analyses.
**Limitations:** estimand is the *contrast* not raw OR (co-switch genes cis-QTL depleted for both); cis-sQTL anchors member-gene splicing, not the coordination itself; bulk GTEx under-samples cell-type-specific splicing.
**Repository sources:** `QTL_ANCHORING_SUMMARY.md`.

#### Finding 4: Disease variants colocalize onto GO-invisible isoform switches; SNCA is the coherent exemplar
**Evidence status:** Supported at set level; suggestive per-gene.
**Key quantitative evidence:** 12/12 splicing-led genes GO-invisible; SNCA risk alleles for LBD (rs7680557) and PD (rs1471483) both raise usage of the same alt-first-exon junction on one switch pair.
**Validation:** GTEx-tissue-matched concordance (9 events, 2-bp junction floor); cross-disease concordance; curated literature layer (4 genes with matching known isoform biology).
**Limitations:** eCAVIAR posteriors modest (max 0.39); per-gene claims suggestive — defensible claim is the *set-level pattern*.
**Repository sources:** `GENETIC_ANCHORING_RESULTS.md`, `deep_dive/`.

#### Finding 5: IsoGraph's real-data modules are reproducible and replicate aging associations
**Evidence status:** Established.
**Key quantitative evidence:** 236/266 chance-trusted; driver ρ 0.77–0.82; 25 cross-cohort aging replications vs ~6 for matched abundance baseline.
**Validation:** split-half stability, size-matched permutation null, cross-cohort replication.
**Repository sources:** `figTrustFunnel`, `tableS7`, AGENTS.md §4.

#### Finding 6: On synthetic ground truth IsoGraph recovers switch modules and residualization buys confound-robustness
**Evidence status:** Established.
**Key quantitative evidence:** module recovery favours IsoGraph in 7 of 15 scenarios — a conditional win (switching-dominated) and not a general one; complete switch-gene detection; figS8 residualization holds under composition/batch/depth.
**Limitations:** RNA-degradation is the honest exception (figS9 abundance fallback).
**Repository sources:** `fig1`, `tableS_benchmark_summary` (supp), `figS8/figS9`.

### 3. Negative, Null, and Sensitivity-Dependent Findings

| Finding | Analysis | Interpretation | Effect on manuscript claims | Sources |
| --- | --- | --- | --- | --- |
| Raw cis-QTL OR < 1 for co-switch genes (both sQTL+eQTL) | qtl_anchoring | Constraint baseline of network genes, NOT the result | Forces reporting the *contrast*, not raw enrichment | `QTL_ANCHORING_SUMMARY.md` |
| sQTL intron-direction concordance at/below perm null (\|rho\| 0.29–0.34, p≈0.76–1.0) | sqtl_concordance | Diagnosed null — single lead sQTL tags near-constant per-tx sign | Directionality carried by coloc instead; reported as underpowered negative | AGENTS.md §6 |
| Disease-SCZ heritability is eQTL-led (sQTL 0.71× ns; eQTL 2.04×) | S-LDSC | Honest asymmetry: SCZ resolves via abundance more than splicing | Splicing-anchoring headline is aging/neurodeg; SCZ disease is the qualifier | `LDSC_SUMMARY.md` |
| SCZD single-tissue sQTL tail (OR 1.37) does not survive pooling | qtl_anchoring_meta | Single-tissue fluke | Trust meta, not one tissue | `QTL_ANCHORING_SUMMARY.md` |
| IsoGraph lowest "both" (pheno+GO) rate (0.074) | baseline_comparison | GO-enrichment low by construction | Do NOT claim "pathways WGCNA misses" | `BASELINE_COMPARISON_SUMMARY.md` |
| Cross-cohort GO replication favours abundance WGCNA (0.077 vs 0.0/0.048) | baseline_comparison | Abundance/GO advantage | Bounds the claim | `BASELINE_COMPARISON_SUMMARY.md` |
| **SCZ single-locus genotype concordance null/underpowered** (28/62, P=0.81) | scz_age_projection layer A | Demoted from headline to supplementary | Headline is convergence + module disruption, not per-locus genotype | `scz_age_projection.py`, `A_gwas_directed_concordance.parquet` |

### 4. Repository Inconsistencies

| Issue | Conflicting sources | Consequence | Required resolution |
| --- | --- | --- | --- |
| GO-invisible gate module count: 8 modules (2 GO-visible) vs 4 modules (all GO-invisible) | Pre-refit prose vs post-refit `go_invisible_gate.parquet` (2026-06-29) | Wrong number in older narrative | **Resolved in-repo:** cite the regenerated parquet (4, all GO-invisible). Ensure no drafted prose still says 8/2. |
| Coloc gene count phrasing: "68 genes" vs "141 events in 68 genes" vs "12+23+33" split | `GENETIC_ANCHORING_RESULTS.md` internal | 12+23+33 = 68 ✓ but "141 colocalized isoform events" vs Fig4D "68 colocalized genes" needs one canonical count in caption | Reconcile event-count vs gene-count wording in Fig 4 caption before submission |
| Target journal: AGENTS.md header + §5 say Cell Genomics; `FIGURE_ORDERING.md` header still says "Nature Methods" | AGENTS.md vs FIGURE_ORDERING.md | Stale label | Update FIGURE_ORDERING.md header to Cell Genomics |
| Fig 3 headline (full 17-analysis) vs common-subset | `qtl_anchoring_meta_contrast{,_common}.parquet` (pheno-sig 1.111 full vs 1.108 on the 8 tissues common to all methods) | Both correct for different subsets | State which subset each number is from in the caption |

### 5. Reproducibility Status

| Analysis | Code | Results | Logs | Environment | Version info | Status |
| --- | --- | --- | --- | --- | --- | --- |
| All real-data layers | Committed CLIs under `isograph_benchmark/real_data/` + SLURM `_h/` | Parquet + `.md` summaries | SLURM logs | isograph env (Py 3.12) | IsoGraph v0.1.5; **no per-run lockfile** | Reproducible; version pinning is per-live-env not lockfile |
| Synthetic benchmark | Committed | fig+table | yes | isograph env | v0.1.5 | Reproducible |
| SCZ age-projection | Committed CLI + wrapper 18 | Complete; committed 8f1d315 | SLURM log | isograph+eqtl envs | v0.1.5 | Genotype dosages controlled-access (gitignored) — layer A not fully re-runnable without TOPMed access |

**[UNRESOLVED — reproducibility]** No per-run dependency lockfile exists; versions read from the live environment. Low-risk but a Cell Genomics Key Resources Table will want explicit versions.

---

## Part II. External Requirements

### 6. Sources Consulted

| Source | Type | Date accessed | Scope |
| --- | --- | --- | --- |
| AGENTS.md §5 re-target notes | Repository planning note | 2026-07-20 | Internal record of Cell Genomics re-target tasks |
| Cell Genomics *Article types* (cell.com) | Target-journal instructions | 2026-07-20 (via search index; pages 403 to direct fetch) | Word/display-item limits per article type |
| Cell Genomics *Information for authors / Journal policies* (cell.com) | Target-journal instructions | 2026-07-20 (via search index) | Data & code availability policy |
| Cell Press *STAR Methods* guide + *Resource availability* asset (cell.com) | Reporting-format guideline | 2026-07-20 (via search index) | STAR Methods section schema + Key Resources Table |
| Manubot documentation | Tooling | Not re-verified this session | Manuscript already Manubot-based |

> Retrieval note: cell.com returns HTTP 403 to direct WebFetch; the requirements below were
> recovered from the Cell Press search index (title + snippet extraction), not from the rendered
> pages. They are accurate as indexed on 2026-07-20 but should be eyeballed against the live pages
> at final formatting. URLs: article types `/cell-genomics/information-for-authors/article-types`;
> policies `/cell-genomics/information-for-authors/journal-policies`; STAR Methods
> `/information-for-authors/star-authors-guide`; resource availability
> `/pb-assets/journals/assets/info-for-authors/resource-availability.html`.

### 7. Mandatory Requirements (VERIFIED against Cell Genomics guidelines, 2026-07-20)

**Article length & display items (research Article — the target type):**

| Requirement | Detail | Repository status | Action needed |
| --- | --- | --- | --- |
| Main-text word limit | **< 8,000 words** excluding references, STAR Methods, supplemental (longer OK after editor discussion) | Drafted; total not yet counted | Count words; trim to <8k |
| Display items (figures **+ tables**) | **≤ 7** combined | **5 planned: Fig 1–4 + Table 1 (=QTL contrast)**; +Table 2 (genes) → 6; +SCZ panel → 7. Synthetic benchmark table DEMOTED to supp (`tableS_benchmark_summary.csv`, done). | Within limit; see §15 |
| Summary/abstract | Structured summary required | Written for NM framing | Re-lead for CG biology payoff |
| Alt article types (for reference) | Short Article ≤4,000 w (incl. legends, excl. STAR Methods); Technology/Resource <7,000 w | — | Article is the right type |

**STAR Methods (no word limit; typeset with main text; NOT copyedited):** required section order —

| Required section | Required subsections | Repository status | Action needed |
| --- | --- | --- | --- |
| **Key Resources Table** (required) | tools, software, cell lines, **original source data for computational studies**; each item must also appear in Method Details | Not built | Assemble KRT: IsoGraph v0.1.5, MOODS, ATtRACT, plink2, GTEx v11 sQTL/eQTL, BrainSEQ, TOPMed, gnomAD, ClinVar, baselineLD v2.2, GWAS sumstats |
| **Resource Availability** (required) | **Lead contact**, **Materials availability**, **Data and code availability** (must include accession numbers + DOIs) | Placeholders only | Name lead contact; write DACA statement with accessions + minted DOIs |
| Method Details (required) | narrative methods; every KRT item referenced here | Drafted per-analysis, not STAR-formatted | Reformat to STAR Method Details |
| Experimental Model and Subject Details (when appropriate) | cohorts, sample sizes, demographics, ethics | Numbers exist; not STAR-blocked | Add cohort/ethics block |
| Quantification and Statistical Analysis (when appropriate) | tests, seeds, thresholds | Fully specified in summaries | Consolidate into this section |

**Data & code availability policy (mandatory before acceptance):**

| Requirement | Detail | Repository status | Action needed |
| --- | --- | --- | --- |
| Code in DOI-minting repo | All original code deposited (e.g. Zenodo) or in supplement **before acceptance**; DOI reported in paper; public by publication date; any embargo lifted | Zenodo staging exists (`zenodo/MANIFEST.tsv`); DOIs are placeholders | Mint Zenodo DOI(s); report in DACA |
| Data deposited/accessible | All data in a repository public by publication or shared by lead contact on request; accessions in DACA | BrainSEQ/GTEx are external controlled/public; TOPMed controlled-access | List accessions; state controlled-access route for TOPMed genotypes |
| DACA statement completeness | Must include **all** accession numbers and DOIs | Not written | Draft comprehensive DACA |

### 8. Recommended Reporting Practices

| Practice | Source | Mandatory or recommended | Repository status |
| --- | --- | --- | --- |
| Deposit code + reproducible pipeline | General open-science | Recommended (often mandatory) | Met — committed CLIs + SLURM wrappers |
| Deposit heavy artifacts off-git | Repository policy | Recommended | Met — Zenodo-bound (`zenodo/MANIFEST.tsv`), git-LFS lean results |
| Controlled-access genotype handling | Data-use policy | Mandatory | Met — TOPMed dosages gitignored, never committed |

### 9. Requirements Not Yet Evaluated

- **Reference/citation style + supplemental-information caps** — the Cell Press reference format and any supplemental figure/table count limit were not recovered from the search index; verify on the live *Revise your manuscript* / *Final file requirements* pages before formatting.
- **Highlights / eTOC blurb / graphical abstract** — whether Cell Genomics research Articles require Highlights (typically 3–4 bullets ≤85 chars) and/or a graphical abstract was not confirmed; check the article-types page directly.
- **Human-subjects / data-use statements** — exact IRB/DUC language for BrainSEQ + GTEx + TOPMed; the *policy exists* (Experimental Model and Subject Details) but the required wording was not audited.
- **DOI minting** for Zenodo/benchmark repo/protocols.io — placeholders only (policy is verified as mandatory; the DOIs themselves are pending).
- **Inclusion & diversity / ethics statements** (Cell Press standard) — not audited.

---

## Part III. Strategic Recommendations

### 10. Recommended Manuscript Framing

**Recommendation:** Keep the existing **complementary-layer** framing, and lead the biological payoff with the **genetic-anchoring convergence**, not method benchmarking.

**Central claim:** IsoGraph surfaces a reproducible, genetically-anchored **isoform-switch regulatory layer** — invisible to abundance/GO pipelines — whose disease-relevant modules are spared at splicing-QTL specifically and onto which disease variants (e.g. SNCA across LBD+PD) colocalize.

**Evidence supporting the recommendation:**
- The strongest, most-controlled result is Finding 3 (splicing-specificity contrast, homogeneous across tissues at I²=0.00, with the matched-baseline null as its internal control + matched-baseline method effect).
- Findings 2+4 give the mechanism (real GO-invisible DTU that disease variants resolve onto).
- Finding 1 keeps it honest (complementary, not superior).

**External requirements affecting the recommendation:** No *requirement* dictates framing. The verified ≤7-display-item limit (figures + tables) is satisfied by the current 5-item plan (Fig 1–4 + the QTL-contrast table), with room for a 6th (SCZ convergence) panel; the biology-led framing is a [RECOMMENDATION] matched to a genomics-journal audience, not a journal rule.

**Reasoning:** This converts a "yet another network method" into a genetics-anchored discovery of a regulatory layer, differentiating it from a WGCNA benchmark and matching a genomics-journal audience.

**Tradeoffs:** Rests on coloc posteriors that are modest; the defensible unit is the *set-level* pattern, and the paper must repeatedly say so.

**Conditions that would change the recommendation:** The SCZ age-projection convergence is **RETRACTED** (P=0.40 against the coloc-tested background; the P=0.0058 used an all-genes background inflated by shared ascertainment), so the "risk loci converge on age-sensitive programs" sub-claim moves OUT of main text. What remains main-text-worthy is that the age-sensitive modules are disrupted in disease (B/C/D), which is a different and independent claim. The remaining swing factor is coloc robustness — the per-gene genetics stays set-level.

**Confidence:** Moderate-to-high (evidence is strong and self-controlled; per-gene genetics is suggestive).

### 11. Alternative Framings

#### Alternative Framing: Methods-first (benchmark-led)
**Central claim:** IsoGraph is a new isoform-switch network method that beats WGCNA.
**Supporting evidence:** synthetic Fig 1 (7/15 scenarios on module recovery), trust funnel.
**Advantages:** clean, quantitative, Nature-Methods-shaped.
**Weaknesses:** directly contradicted by Finding 1 (not globally superior); would overclaim.
**Why not preferred:** the repository's own scope-bounding analysis forbids a superiority headline.
**Evidence that could make it preferable:** none available — this is a settled scope decision.

#### Alternative Framing: SNCA/synucleinopathy-led (single-gene vignette)
**Central claim:** A shared SNCA alt-first-exon switch underlies LBD+PD risk.
**Supporting evidence:** SNCA cross-disease concordance.
**Advantages:** concrete, memorable.
**Weaknesses:** single locus, modest CLPP; overweights one gene.
**Why not preferred:** the strength is set-level coherence, not one locus.
**Evidence that could make it preferable:** high-posterior functional validation of the SNCA switch (not in repo).

### 12. Proposed Claims Hierarchy

| Proposed claim | Repository evidence | Evidence status | Narrative role | Principal caveat |
| --- | --- | --- | --- | --- |
| Splicing-QTL specifically anchor the GO-invisible switch layer (IsoGraph-only) | Finding 3 | Established | **Headline (Fig 3)** | Estimand is the contrast |
| Disease variants colocalize onto GO-invisible switches (SNCA exemplar) | Finding 4 | Supported/set-level | **Payoff (Fig 4)** | Modest posteriors |
| Disease switch modules are real GO-invisible DTU | Finding 2 | Established | Mechanism bridge | One cohort |
| Modules are reproducible + replicate aging | Finding 5 | Established | Trust (Fig 2) | — |
| Method recovers switch modules on truth | Finding 6 | Established | Foundation (Fig 1) | Synthetic |
| IsoGraph is complementary, not superior | Finding 1 | Established | Honest bound (Supp) | Deliberately limiting |
| ~~SCZ-risk loci converge on age-sensitive programs~~ (RETRACTED) | SCZ projection | **Not supported** — P=0.40 vs the coloc-tested background (was P=0.0058 vs an all-genes background) | Fig 4E now carries the size-matched null instead (done 2026-09-03); full version is S-real-10 | Ascertainment: MAGMA anchoring and coloc testing select on the same GWAS |
| ~~Switch genes prefer splicing over expression colocalization, per gene~~ (NULL) | Per-gene coloc modality contrast | **Not supported** — splicing share of discordant genes ~31% in all four gene pools; locus-matched switch 20/46 vs non-switch 88/188, Fisher exact P=0.88; conditional posterior null in the switch arm (P=0.39) | Report as a negative control in Results; **no figure** | Within-gene design means module selection cancels; the expression-favouring skew is GTEx eQTL discovery power. Does **not** contradict Fig 3 — different estimand — but does not corroborate it either |
| Age-sensitive switch modules are disrupted in SCZ | SCZ projection | **Established (module-level)** B 4/10, D 3/10, C 10/10 | Main-text Results 6 | Independent of the retracted convergence test |

### 13. Recommended Results Outline

Order is dependency-driven: method works → modules trustworthy → modules are real GO-invisible DTU → genetically anchored → variants resolve onto them → (new) risk loci converge in disease. Honest bound sits in Supp.

```markdown
### Results 1: IsoGraph recovers isoform-switch modules on synthetic ground truth (Fig 1)
- Evidence: module recovery favours IsoGraph in 7/15 scenarios; complete switch-gene detection; figS8/S9 confound robustness.
- Role: foundation. Placement: first. Transition: "but does it recover *trustworthy* structure on real brain?"
- Confidence: High.

### Results 2: The real-data modules are per-module trustworthy and replicate aging (Fig 2)
- Evidence: 236/266 chance-trusted; ρ 0.77–0.82; 25 vs 6 aging replications.
- **Declare the model with the number.** The replication count is `25/130 matched module
  pairs, complement covariate mode, linear age model` (matching permutation p = 0.014,
  null 16.47 ± 3.50). `complement` is the pre-specified mode — the IsoGraph fit already
  residualized discovery covariates, so the `full` arm is double adjustment and its
  collapse to 2/130 is expected arithmetic, not a refutation. State the mode in the
  sentence, not only in Methods.
- **The linear age model is a justified choice, not a convenience.** The df = 3 spline arm
  gives 15/130 (p = 0.15) at the same covariate mode, and a reader is entitled to know
  whether that is the spline correcting a misspecified line or the line being cheaper. It
  is the latter, tested directly: across all 266 modules in all six cohort × region fits,
  **not one** shows detectable curvature in its age trajectory (Wald second-difference
  contrast inside the spline's own fit; 0 at nominal p < 0.05 where chance alone predicts
  13.3; smallest p = 0.71; median χ² = 0.007 on 1 df). The contrast is well powered
  precisely because the shared uncertainty of the projected effects cancels. Nor does the
  spline find different modules: 4/266 are spline-only significant against 78/266
  linear-only — the 4 is the chance rate. Aging is linear in these data at module
  resolution, so the extra 2 df buy nothing and cost power.
- **Placement:** spline arm → Supplement as a declared sensitivity, cited from the Results
  sentence in one clause ("the count is lower under a df = 3 spline, which these data do
  not support — Supp. Table X"). Do **not** present it as a co-equal arm; that would imply
  a model choice the data reject.
- Evidence CLI: `isograph_benchmark/real_data/age_model_curvature.py` (+
  `03_module_trust/_h/14.age_model_curvature.sh`) →
  `03_module_trust/_m/stability/module_trust/age_model_curvature__isograph.{parquet,json}`.
  Refits nothing; reads the committed `age_spline`/`age_linear` tables.
- Role: answers "fine partition = noise?". Transition: "what *are* the disease-associated ones?"
- Confidence: High.

### Results 3: The disease switch modules are real, GO-invisible DTU (S-real-3 support)
- Evidence: 4/4 SCZD pheno-sig GO-invisible; consequence ≥ background; productive UTR/CDS remodeling (S-real-5).
- Role: mechanism bridge. Transition: "are they genetically real?"
- Confidence: High (one cohort — say so).

### Results 4: Splicing-QTL specifically anchor the phenotype-associated layer — an IsoGraph method effect (Fig 3) [HEADLINE]
- Evidence: contrast 1.111 pheno-sig (p=3.6e-4) / 1.068 all_modules (p=1.3e-5); IsoGraph-only vs matched WGCNA (null, p≥0.41) — the primary internal control; S-LDSC heritability.
- **Do not claim GO-invisible localisation here:** GO-invisible (1.068, p=0.077) and GO-visible (1.084, p=0.050) are indistinguishable. The DTU-without-DGE claim is Results 3's content gate, not this genetic one.
- Role: headline orthogonal validation. Transition: "which variants, on which switches?"
- Confidence: High.

### Results 5: Disease variants resolve onto specific isoform switches (Fig 4)
- Evidence: 12 splicing-led all GO-invisible; SNCA LBD+PD alt-first-exon; RBP regulons (S-real-6); clinical consequence (S-real-7).
- Role: biological payoff. Transition: "do risk loci converge as a program in disease?"
- Confidence: Moderate (set-level).

### Results 6 (NEW): Age-sensitive switch programs are disrupted in schizophrenia

_(Retitled 2026-08-30. The former title asserted convergence of SCZ-risk loci on these
programs; that test is retracted — see the RETRACTED note below. What survives is module
disruption in disease, which does not depend on it.)_
- Evidence: ~~convergence hypergeom P=0.0058 (15/32 loci in anchored modules)~~ **RETRACTED 2026-08-30, P=0.40 on the coloc-tested background**; disruption B 4/10, D 3/10, C 10/10; named regulators (M002 SNRNP70/ZCRB1; M008 ZC3H10/RBM14/CELF5; M006 DDX58/ADAR/YTHDC1). Single-locus genotype concordance null (28/62, P=0.81) → Supp.
- Role: **superseded.** The convergence result was retracted 2026-08-30 (ascertainment); it is not significant on the coloc-tested denominator. Fig 4E now carries the *null* instead, as of 2026-09-03. Do not promote this to a Results unit.
- Confidence: Moderate — module-level convergence is significant; per-locus genotype resolution is underpowered (n~62), disclose as such.

### Results (bound): The aging switch layer is composition-robust; the schizophrenia layer is not (Fig 5)
- Evidence: MuSiC cell-type adjustment of the de-confounded gene-level switch test.
  Aging caudate 43 → 17 composition-unique genes (40% retained), aging DLPFC 8 → 15
  (adjustment is not uniformly attenuating — this is the internal control), aging
  hippocampus 0 → 0 (contributes nothing either way; say so rather than counting it as a
  pass). SCZD caudate **34 → 2 (6%)**. Module member genes are not marker bags
  (0/3,716 depleted, 33/3,716 enriched). GTEx replication is region-dependent: 6/8
  deconvolved regions retain a composition-robust signal, the two cortical regions
  collapse to zero.
- **Decision (P1-6): concede the disease arm as composition-entangled, and say why that is
  a biological hypothesis rather than a failure.** The aging claim is the one that
  survives adjustment and is the one to lead with. For schizophrenia, write the two
  readings explicitly and do not adjudicate between them from bulk data:
  - *Confounder reading:* cell-type proportions differ between cases and controls for
    reasons incidental to the switch layer (sampling, dissection, agonal state), and
    adjustment is the correct control. Then 2/34 is the honest count.
  - *Mediator reading:* the composition shift **is** part of the disease — schizophrenia
    post-mortem cortex and striatum carry real differences in interneuron and glial
    proportions — in which case adjusting for it removes signal that lies on the causal
    path, and 34 → 2 is over-adjustment, not deflation of a false positive.
  Bulk RNA-seq cannot separate these: proportion and per-cell composition are confounded
  in the same measurement by construction, and no covariate arrangement inside a bulk
  design breaks that. This is a limit of the assay, not of the model.
- **Future work, stated in the Discussion:** resolving it requires the switch layer to be
  measured *within* cell types — single-nucleus (or single-cell) RNA-seq of the same
  regions, with enough depth or long reads to quantify isoform usage per nucleus. If the
  switches persist within matched cell types, the mediator reading holds and the disease
  layer is real; if they disappear once cell type is fixed, the confounder reading holds
  and the disease claim should not have been made. Name this as the decisive experiment
  rather than gesturing at "further work" — it is a specific, falsifiable follow-up, and
  saying so converts the collapse from a weakness into a scoped next step.
- **Do not** make any SCZD-specific composition-independent claim on these data. The
  disease-side claims that stand are the ones that do not route through composition:
  the GO-invisible content gate (Results 3), the module-level genetic anchoring
  (Results 4), and module disruption in disease (Results 6).
- Role: the confound reviewers will press hardest, answered in the open. Placement: main
  (Fig 5), with the mediation discussion in the Discussion.
- Confidence: Moderate — direction clear, counts small enough that a handful of genes
  moves the percentage substantially; disclose the n.

### Results (bound): IsoGraph is complementary, not globally superior (S-real-1)
- Evidence: three-baseline rates; abundance/switch separability.
- Role: honest scope. Placement: Supplement, referenced from Discussion.
- Confidence: High.
```

### 14. Recommended Methods Outline (STAR-target)

```markdown
### Methods: IsoGraph model + module inference
- Evidence: run_models.py, vae.py, residualize.py; res 5.0, giant-cap off, seed 13.
- Requirement: STAR "Method Details"; report VAE arch, single-LR + estimability, Leiden res.
- Missing: per-run lockfile.

### Methods: Confound residualization policy (discovery-only)
- Evidence: covariate-decouple (raw feature_scores) commit c958174; residualization_qc.parquet.
- Requirement: state the discovery-vs-inference split explicitly.

### Methods: Trust funnel / stability
- Evidence: 03_module_trust/; split-half + perm null.

### Methods: QTL anchoring + contrast meta
- Evidence: qtl_anchoring.py (+--method), qtl_anchoring_meta.py; power-matched logistic; IVW+DL.
- Requirement: define the contrast estimand; GTEx v11 access date (KRT).

### Methods: Colocalization + deep-dive
- Evidence: coloc (SuSiE + eCAVIAR), gene_deep_dive.py; 5 GWAS sources need citekeys.
- Missing: pinned citekeys for eCAVIAR/GTEx v11/SuSiE/gnomAD ([citation needed] markers).

### Methods: S-LDSC; MAGMA resolution; switch/clinical consequence; RBP regulons
- Evidence: respective CLIs + wrappers; ATtRACT+MOODS for RBP.

### Methods: SCZ age-projection (NEW)
- Evidence: scz_age_projection.py + wrapper 18; TOPMed genotypes controlled-access.
- Requirement: data-use statement; note dosages not shared.
```

### 15. Figure and Table Recommendations

Figure order already follows the argument (`FIGURE_ORDERING.md`). **Main display items (≤7 cap): Fig 1–4 + main Table 1 (QTL contrast) = 5.** With the synthetic table demoted, the biology tables become the main-text Table 1/2. Deltas only:

| Proposed item | Evidence | Claim supported | Recommendation | Reasoning | Confidence | Placement |
| --- | --- | --- | --- | --- | --- | --- |
| Fig 1–4 as ordered | benchmark + trust + QTL + anchoring | Findings 6/5/3/4 | **Keep as-is** | Dependency-correct, biology-led | High | Main |
| ~~Synthetic benchmark summary~~ (`tableS_benchmark_summary.csv`) | 01_synthetic_benchmark/03_metrics | Finding 6 | **MOVED TO SUPPLEMENT — DONE** | Demoted from main; regenerated clean (6 core scenarios × 6 methods, no NA rows) and renamed via `synthetic_benchmark.R` (`ISOGRAPH_TABLES_ONLY=1`). References updated in `FIGURE_ORDERING.md` + `01_synthetic_benchmark/03_metrics/README.md`. | High | **Supp ✓** |
| **Main Table 1 — sQTL/eQTL splicing-specificity contrast (NEW, built)** | qtl_anchoring_meta contrast | Finding 3 | **PRIMARY main biology table** | The p-value-bearing anchor: contrast **1.111 pheno-sig (p=3.6e-4)** / 1.068 all_modules (p=1.3e-5), IsoGraph-only vs matched WGCNA (null everywhere, p≥0.41). Does **not** localise to GO-invisible (1.068, p=0.077 vs GO-visible 1.084, p=0.050) — the legend in `assemble_main_tables.py` says so explicitly and must not be softened. **Reproducible:** `manuscript/_h/assemble_main_tables.py` → `table2_qtl_specificity_contrast.{csv,md}` (file keeps `table2_` stem; main-text number = Table 1). | High | **Main** |
| **Main Table 2 — Splicing-led colocalized genes (NEW, built; reframed)** | deep-dive panel + literature | Finding 4 | **Main companion OR keep in Supp (S8/S9)** | Per-gene resolution; **CLPP posteriors are individually modest** — coloc threshold is eCAVIAR CLPP≥0.01 ("strong" ≥0.05), only **4/12 clear 0.05** and only CTSH (0.39) is substantial. Caption states the claim is *set-level coherence*, NOT per-locus significance (which lives in the contrast table + S-LDSC). Max CLPP carries confidence stars (`*` >0.01, `**` >0.05, `***` >0.10 → 1×`***`, 3×`**`, 8×`*`). Reproducible: same builder → `table3_splicing_led_genes.{csv,md}`. | Moderate | **Main or Supp** |
| S-real-1 (baseline rates) | baseline_comparison | Finding 1 | **Keep in Supp** | Bounds, not advances | High | Supp |
| **SCZ convergence panel** | SCZ projection (**RETRACTED**, P=0.40 on the correct background) | Results 6 | ~~Remove from Fig 4E~~ — **DONE 2026-09-03**: Fig 4E now shows the size-matched-null test for all five traits (null in 10/10 cells), lifted from S-real-10; the full three-panel version stays as S-real-10 (`figColocConvergence`) | A count-only panel cannot support convergence: it shows neither the size-matched null nor the tested-pool denominator | — | **Fig 4E (null panel) + S-real-10** |
| RBP intronic/combined scope | new intronic scan | S-real-6 mechanism | **Fold into S-real-6 or its table (S10)** | Strengthens regulon call with intronic niche | Moderate | Supp table |

> **Note — synthetic Table 1 → Supplement: DONE.** Regenerated clean (six core accuracy scenarios ×
> six main methods, no NA rows — the multiplex/scale scenarios that produced NA are correctly
> excluded) and renamed to `01_synthetic_benchmark/03_metrics/_m/tableS_benchmark_summary.csv` via a new
> `ISOGRAPH_TABLES_ONLY=1` fast path in `synthetic_benchmark.R` (no figure churn). Old
> `table1_benchmark_summary.csv` `git rm`'d; references updated in `FIGURE_ORDERING.md` +
> `01_synthetic_benchmark/03_metrics/README.md`.
>
> **Display-item math (≤7 cap, figures + tables).** Recommended main set = Fig 1–4 + **main Table 1
> (QTL contrast)** = **5**; add **main Table 2 (splicing-led genes)** → 6; conditional SCZ panel → 7.
> If both biology tables + the SCZ panel are wanted, that is exactly 7 — at the cap, so either keep
> the genes table in the supplement (it already exists as S8/S9) or fold the SCZ result into Fig 4
> rather than a standalone item. Recommendation: **QTL-contrast table main (statistical anchor),
> genes table main only if a per-gene table is wanted for the biology payoff; otherwise Supp.**
>
> **CLPP honesty (why Table 2 leads, not Table 3).** The colocalization posteriors behind the
> splicing-led gene set are individually weak (eCAVIAR CLPP≥0.01 inclusion; only 4/12 ≥0.05; only
> CTSH substantial at 0.39). The genetic-anchoring *significance* therefore rests on the set-level
> QTL specificity contrast (Table 2, p=1.5e-4 / 2.6e-3) and S-LDSC partitioned heritability — not
> on per-gene CLPP. Any main-text sentence must not present the per-gene colocalizations as
> "significant"; frame them as resolved candidates whose strength is cross-disease/GO-invisible
> coherence. This corrects the earlier draft caption.

### 16. Introduction Recommendations

```markdown
### Introduction Unit: The gap — abundance pipelines miss isoform regulation
- Scope: brain aging + psychiatric/neurodegenerative disease; bulk RNA-seq; DTU vs DGE.
- Literature: WGCNA/co-expression; sQTL vs eQTL; isoform switching in brain.
- Argument: a complementary switch layer could carry disease genetics invisible to DGE/GO.
- Avoid: claiming superiority over WGCNA.

### Introduction Unit: Why genetic anchoring is the right test
- Argument: motif/module coherence alone is weak; QTL/coloc anchoring makes the layer falsifiable.
- Connect to central framing: sets up the two-internal-null design.
```

### 17. Discussion Recommendations

```markdown
### Discussion Unit: What the complementary layer adds
- Interpret Findings 1–4 jointly; emphasize DTU-without-DGE content + genetic specificity.
- Alternative explanation to rebut: "it's just abundance re-labeled" (answered by separability + GO-invisible gate).

### Discussion Unit: Limitations (disclosed, not hidden)
- Degradation is the residualization exception (figS8/S9).
- Coloc posteriors modest → set-level claims only.
- RBP motifs are predictions, not measured binding.
- SCZ resolves eQTL-led (splicing headline is aging/neurodeg).
- sQTL intron-direction concordance is an underpowered null (coloc carries direction).
- Single-cohort GO-invisible gate; bulk GTEx under-samples cell-type splicing.

### Discussion Unit: Speculation (bounded)
- SNCA alt-first-exon as a shared synucleinopathy 5′-regulatory mechanism — flag as hypothesis for functional follow-up.
```

### 18. Supplementary and Omitted Analyses

| Analysis | Evidence status | Recommendation | Reason |
| --- | --- | --- | --- |
| Three-baseline rates | Established | Supp (S-real-1) | Bounds claim |
| Switch/RBP/clinical consequence | Established/moderate | Supp (S-real-5/6/7) | Mechanism support |
| sQTL intron-direction concordance | Null (diagnosed) | Supp text / omit from tables | Underpowered by construction |
| RBP intronic scope | Moderate | Supp table addendum to S-real-6 | Extends niche, not a headline |
| SCZ single-locus genotype concordance | Null/underpowered | Supp | Demoted per prior decision |
| Module-level degradation QC flag (#10) | Not done | Omit; disclose as limitation | Explicit WON'T-DO decision |

### 19. Additional Work

| Proposed work | Category | Evidence gap addressed | Effect | Priority |
| --- | --- | --- | --- | --- |
| ~~Decide Results-6 / SCZ-panel placement~~ — **CLOSED 2026-09-03** | Impact-enhancing | Disease-convergence sub-claim **withdrawn** (ascertainment) | Fig 4E now shows the size-matched-null test; S-real-10 keeps the full treatment | — |
| Retrieve + reconcile Cell Genomics author instructions | Submission-critical | Unverified requirements | Unblocks formatting | **High** |
| STAR Methods conversion + Key Resources Table | Submission-critical | Format compliance | Required for submission | **High** |
| Pin `[citation needed]` citekeys (eCAVIAR, GTEx v11, SuSiE, gnomAD, baselineLD) | Submission-critical | Missing citations | Removes placeholders | High |
| Mint Zenodo/protocols.io DOIs | Submission-critical | Data/code availability | Required | High |
| Reconcile stale counts (8→4 modules; 141 events vs 68 genes; NM→CG header) | Reviewer-defense | Internal inconsistency | Prevents reviewer confusion | Medium |
| Per-run dependency lockfile | Reviewer-defense | Version reproducibility | Strengthens KRT | Low |
| RBP experimental validation | Out of scope | Motif ≠ binding | Would upgrade "candidate" | Out of scope |

### 20. Final Strategic Recommendation

**[EVIDENCE] Current repository state:** A near-complete, self-controlled body of work establishes IsoGraph as a reproducible, complementary DTU-without-DGE layer whose phenotype-associated modules are specifically anchored to splicing genetics (contrast 1.111, p=3.6e-4; matched-baseline null, IsoGraph-only) and onto which disease variants colocalize (SNCA cross-disease). The DTU-without-DGE character of those modules is a *content* result (GO-invisible gate), not a genetic one — the two claims are now carried by separate evidence. Figures are ordered, a Results section is drafted, and the supplement is wired. The one honest bound (not globally superior to WGCNA) is settled and belongs in the supplement.

**[REQUIREMENT] Submission constraints (verified 2026-07-20):** Cell Genomics research Article — **<8,000 words** (excl. refs/STAR Methods/supp), **≤7 display items counting figures AND tables** (current plan = Fig 1–4 + QTL-contrast table = 5, with room for a 6th SCZ panel; synthetic benchmark table demoted to supp), structured summary. **STAR Methods** required with a **Key Resources Table** + a **Resource Availability** block (Lead contact / Materials availability / Data and code availability). **Mandatory before acceptance:** all original code in a DOI-minting repo (Zenodo) with the DOI reported, and a data-and-code-availability statement listing every accession + DOI. Still to verify on the live pages: reference style, supplemental caps, Highlights/graphical-abstract, and inclusion/ethics wording.

**[RECOMMENDATION] Preferred manuscript strategy:** Biology-led, genetics-anchored complementary-layer framing; Results ordered method→trust→GO-invisible DTU→splicing-QTL anchoring (headline)→variant resolution→module disruption in disease (B/C/D; the *convergence* test is retracted, P=0.40); keep the superiority bound in the supplement; state the set-level nature of the genetic claims repeatedly.

**[RETRACTED 2026-08-30; was RESOLVED 2026-07-20] Former highest-priority issue:** The SCZ age-projection convergence is **not supported** (hypergeom P=0.40 against the coloc-tested pool; the reported P=0.0058 used an all-genes background confounded by shared ascertainment). "SCZ-risk loci converge on age-sensitive switch programs" is **withdrawn** as a main-text unit. The surviving main-text claim is the weaker and separate "age-sensitive switch modules are disrupted in SCZ" (B/C/D), with named candidate RBP regulators and the null single-locus genotype layer disclosed as supplementary. Committed to `main` (8f1d315).

**[UNRESOLVED] Highest-priority open issue now:** Display-item budget + the Cell Genomics re-target mechanics (STAR Methods, Key Resources Table, DACA with minted DOIs). With Results 6 promoted, main items could reach Fig 1–4 + Table 2 + (Table 3?) + SCZ panel — decide which of {Table 3, SCZ-as-standalone} stays main vs supplement to hold ≤7.
