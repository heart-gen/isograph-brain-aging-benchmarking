# Validation plan — what would move the isoform-switch layer up a journal tier

**Status:** draft for co-author discussion, 2026-09-09.
**Scope:** the two experiments that would convert a well-controlled set-level correlation into a
mechanism. Everything quoted here is recomputed from committed repository outputs; provenance
for each number is in the appendix.

---

## 0. Where the evidence actually stands

### 0.1 What carries the paper

| Result | Effect | Control that makes it credible |
| --- | --- | --- |
| Splicing-specificity contrast, phenotype-associated modules | **1.111** (95% CI 1.048–1.177), p = 3.6e-4, I² = 0.23 | Matched WGCNA baselines on *identical* features: null (switch-only 1.015 p = 0.74; multiplex 1.044 p = 0.070) |
| Same contrast, IsoGraph-only method effect (13 shared tissues) | **1.112**, p = 6.5e-4 | Isolates VAE + Leiden inference from the feature matrix |
| Long-read confirmation, genetically anchored pairs (eCAVIAR layer) | switch-like **0.453 vs matched null 0.252**, p = 5e-4 (n = 53 detected pairs) | Abundance-decile-matched IsoGraph switch pairs from *other* genes, *same samples*, *same code* |
| …restricted to usably expressed anchored isoforms | **0.600 vs 0.304**, p = 5e-4 (n = 30) | as above |
| Long-read confirmation, anchored pairs from the **signal-level nominations** (2026-09-11) | switch-like **0.565 vs 0.238**, p = 5e-4 (n = 124 detected pairs, 24 genes); usable only **0.672 vs 0.298** (n = 61) | as above; a separate arm, not pooled with the eCAVIAR pairs. UNC13A and PICALM cannot enter: neither has a switch pair in the tissue where it colocalizes |
| Module trust funnel | **250/266** chance-trusted (94%), 6 regions | Split-half permutation null |
| Cross-cohort aging replication | **23/130**, matching-permutation p = 0.034 | Independent cohorts, different quantifiers |
| Switch consequence | UTR-remodeled **1.279** (9/9 regions), CDS-changed **1.044** (9/9); NMD, biotype and coding-status all *depleted* | Within-gene permutation null |

The consequence row is the one that dictates assay design, and it is usually skimmed: **these are
productive UTR/CDS remodeling events, not decay**. An assay that only measures protein-coding
change is aimed at the smaller of the two enriched classes.

### 0.2 The three counter-lines

These are not reviewer speculation — they are the project's own results, and they converge on the
same target from three directions.

| # | Counter-line | Number | What it threatens |
| --- | --- | --- | --- |
| **1** | Per-gene sQTL-vs-eQTL coloc contrast is null | Locus-matched switch 20/46 vs non-switch 88/188, **Fisher p = 0.88**; splicing share ~31% in all four gene pools | The set-level effect has no per-gene counterpart, yet the proposed mechanism is per-gene |
| **2** | The GO-invisible arm is null | **1.068**, p = 0.077 (was 1.172, p = 2.3e-5 before the 2026-08-29 input refresh) | The distinctive framing collapses to the blander "phenotype-associated modules" |
| **3** | S-LDSC points at expression, not splicing | sQTL coefficient clears nominal p < 0.05 in **1 of 6** trait-contexts (AD, p = 0.0355) with **no multiple-testing correction applied**; eQTL clears in 4 of 6. Joint model, SCZ: sQTL p = 0.931 vs eQTL p = **0.0064** | For a splicing paper, the partitioned-heritability evidence favours the abundance annotation |

Counter-line 3 is the least discussed internally and the most dangerous externally, because
S-LDSC is the size-robust test the field trusts and it is currently summarised in-repo as
supporting the switch layer without stating that the splicing arm does not survive correction.

**Counter-line 3 is a GTEx power problem, not an unfixable one.** The first draft called it
unrescuable; that was too strong. It cannot be fixed *within GTEx*, but BrainSEQ carries 452–500
genotyped donors per region against GTEx brain's ~181–300, in the same cohort the switches were
discovered in — see §0.4.2. Disclosure plus the §3.4 framing fix is the fallback if that fails,
not the first move.

---

## 0.4 Experiment 0 — exhaust the computational genetics first *(added 2026-09-10)*

**This should run before either wet-lab arm, and most of it is already computed.** The plan as
first drafted went straight to the bench; that was premature.

### 0.4.1 The anchored gene list was built on the wrong statistic — now re-anchored

The 12-gene list in §1.1 comes from **CLPP** (eCAVIAR-style), where the best non-CTSH value is
0.093. *(Rewritten 2026-09-11.)* The first version of this section re-anchored on the
`coloc.abf` posteriors in `coloc_modality_contrast/` and read the result as splicing-specific
recovery of UNC13A and PICALM. The co-author review asked for signal-level colocalization, an
event-level audit, and sQTL-preferential rather than splice-specific language. All three are now
implemented, and the UNC13A / PICALM reading did not survive them.

**Estimator hierarchy: `coloc.susie` > `coloc.abf` > CLPP.** `coloc.susie` on re-fine-mapped GTEx
v11 QTL (a re-fit credible set is kept only where it agrees with GTEx's own), `coloc.abf` where
either side does not fine-map (flagged with a `fallback_reason`, never dropped), and CLPP kept as
orthogonal sensitivity evidence. Only 230 of 579 GWAS loci fine-map and 7.0% of (cell, modality)
rows are scored at signal level, so the fallback is the majority estimator.

| | CLPP (§1.1 list) | signal-level hierarchy, all-introns arm (PP4 sQTL / eQTL, estimator) |
| --- | --- | --- |
| SNCA (LBD) | 0.038 | **0.974 / 0.113**, `coloc.susie`, 11/13 tissues |
| PGS1 (ALS) | 0.093 | **0.976 / 0.165**, `coloc.abf` (no GWAS credible set), 13/13 |
| CDIP1 (SCZ) | 0.028 | 0.949 / 0.854, `coloc.susie` |
| TPCN1 (AD) | 0.011 | 0.861 / 0.819, `coloc.susie` |
| CTSH (AD) | 0.386 | 0.989 / 0.993, `coloc.abf` — colocalizes for *both* modalities |

**42 gene × trait cells (40 genes) reach PP4_sQTL ≥ 0.8, and 14 are sQTL-preferential**
(PP4_eQTL < 0.5). That is sQTL-preferential colocalization under prespecified thresholds — strong
evidence for sQTL colocalization without comparable evidence for eQTL colocalization — and **not**
splice-specific mediation: a low eQTL PP4 can reflect no eQTL, no power, a different causal
architecture, or multiple signals.

| Gene | Trait | PP4 sQTL | PP4 eQTL | estimator | call holds from p12 |
| --- | --- | --- | --- | --- | --- |
| GABBR2 | SCZ | 0.976 | 0.085 | susie | 1e-6 |
| PGS1 | ALS | 0.976 | 0.165 | abf | 1e-6 |
| PLCB2 | SCZ | 0.975 | 0.157 | abf | 5e-6 |
| SNCA | LBD | 0.974 | 0.113 | susie | 5e-6 |
| SNCA | PD | 0.974 | 0.382 | susie | 5e-6 |
| TXNDC15 | ALS | 0.969 | 0.379 | abf | 5e-6 |
| **UNC13A** | ALS | 0.961 | 0.377 | susie | 5e-6 |
| WIPI2 | ALS | 0.893 | 0.478 | abf | 5e-6 |
| KLC1 | SCZ | 0.882 | 0.033 | susie | 1e-5 |
| NDUFS3 | AD | 0.864 | 0.099 | abf | 1e-5 |
| ASB3 | SCZ | 0.859 | 0.478 | abf | 1e-5 |
| SPAG9 | AD | 0.839 | 0.298 | abf | 1e-5 |
| **PICALM** | AD | 0.818 | 0.491 | abf | 1e-5 |
| AZI2 | PD | 0.803 | 0.145 | abf | 1e-5 |

**UNC13A and PICALM are not recovered mechanisms.** The event audit
(`05_genetic_anchoring/_m/locus_event_audit/susie_all_introns/`) walks GWAS signal → sQTL signal →
intron phenotype → driver transcript against curated, coordinate-verified literature events:

| Tier | Loci |
| --- | --- |
| `known_mechanism_recovered` | none |
| `context_distinct_splice_colocalization` | **UNC13A** — colocalizes at chr19:17,630,750-17,632,782 in cerebellum and cerebellar hemisphere; the TDP-43 cryptic-exon intron 20 (chr19:17,641,557-17,642,844) is testable in exactly those two tissues and does not colocalize |
| `disease_locus_splice_linked` | **PICALM** — `coloc.abf` only on the primary grid (locus over the 12,000-SNP guard), at an intron ~45 kb from the TWAS-anchored AD event, which could not promote in any case. In the 30,000-SNP sensitivity arm it reaches `context_distinct` at PP4 0.813: LD-robust, but prior-sensitive (0.303 at p12 = 1e-6) |
| `novel_splice_led_candidate` | the other 40, including SNCA |

**Caveat to state plainly:** conditioning on PP4_sQTL ≥ 0.8 and then observing PP4_eQTL < 0.5 is a
selection, so 14/42 is not an unbiased splicing-specificity estimate. The unbiased test is the
paired within-gene comparison, and it shows no splicing preference at signal level either
(representative arm, pooled: 21 splicing-only vs 43 expression-only discordant genes, McNemar
P = 0.008 in the *expression* direction; conditional-posterior Wilcoxon P = 0.33). These 14 are
*locus nominations*, which is a different and legitimate job.

### 0.4.2 BrainSEQ junction QTLs — ~~the highest-value item on this page~~ run; switch-QTL arm complete, junction layer dropped

> **Outcome (2026-09-11/12).** The switch-QTL arm was mapped (`S_g` vs `A_g`; 420/360/362
> all-samples and 229/169/175 EA-only donors), both pre-specified checks passed, and in-sample-LD
> coloc followed. It runs *against* splicing: switch-QTL signal is largely shared with abundance
> (π₁(A | S) 0.77–0.89 vs π₁(S | A) 0.29–0.41) and per-gene coloc favours abundance, 5
> switch-only vs 33 abundance-only (P = 4.3e-6). The discovery cohort does not rescue
> counter-lines 1 and 3. The junction-level layer was run and shelved, then **dropped on
> 2026-09-12** (PI decision): it is an event-naming interpretability layer on swQTL-positive genes
> only, its within-gene STAR arm found zero QTLs, and restoring it with a local denominator would
> not change the genetic claim now that splicing specificity is supporting evidence only. The text
> below is the original rationale, kept as the record. See `reports/pi/08_signal_level_genetics.md`.

Counter-line 3 was described in the first draft as unrescuable because it is a GTEx bulk sQTL
power problem. That is true of GTEx and **false of the project as a whole**: BrainSEQ genotypes
are already in the repository (`inputs/raw/brainseq/genotypes/TOPMed_LIBD.{pgen,pvar,psam}`),
alongside the PSI/junction tables the SNCA test already used.

| | GTEx brain | BrainSEQ |
| --- | --- | --- |
| n per region | ~181–300 | junction phenotypes 452–500 (caudate 487, DLPFC 500, hippocampus 452); **switch-QTL arm 420 / 360 / 362 genotyped donors** (caudate / DLPFC / hippocampus, `Age > 13`, not dropped) |
| Cohort | different from discovery | **same cohort the switches were discovered in** — same-tissue genetic anchoring, *not* independent replication |
| Ancestry | predominantly European | roughly half African-ancestry — anything colocalized against the EUR GWAS needs the EA-only arm |
| Regions | broad, shallow | caudate / DLPFC / hippocampus, matched to the aging arm |

Mapping junction QTLs in BrainSEQ gives ~1.7–2× the donors, in the same tissue, on the same
junctions, in the cohort that generated the hypothesis. It is the single most direct attack on
both counter-line 1 (per-gene power) and counter-line 3 (S-LDSC splicing power), it needs no
new samples, and it is compute-only.

### 0.4.3 SMR + HEIDI — corroboration beneath coloc, not a co-equal test

*(Revised 2026-09-11 per co-author review; built as `05_genetic_anchoring/_h/27.smr_heidi.sh`.)*
Scoped to the signal-level nominations, on GTEx v11 QTL with 1000G EUR LD. What it adds, and what
it may not be read as:

- **`b_SMR` is a signed ratio-type estimate** relating genetically predicted molecular phenotype
  to disease under a single-causal-variant, no-pleiotropy model. It does **not** establish causal
  direction in the biological sense and cannot distinguish causality from horizontal pleiotropy.
  For an intron-excision probe the sign is also compositional within its LeafCutter cluster.
- **HEIDI** tests compatibility with one shared causal variant versus distinct linked variants.
  **Failing to reject HEIDI is not proof of sharing**, especially at GTEx brain n.
- **A HEIDI rejection does not overrule a strong multi-signal colocalization**, and SMR
  significance does not promote a locus coloc did not support. Disagreements are reported as
  disagreements (`agreement` categories in `SMR_HEIDI.md`).
- **BrainSEQ QTL are used only from the EA-only mapping**: the cohort is roughly half
  African-ancestry, and its LD does not match the EUR GWAS or the 1000G EUR panel. That arm
  exists since 2026-09-11 (`05_genetic_anchoring/_m/brainseq_switch_qtl/ea_only/`), and SMR +
  HEIDI has been run on it (`_m/smr_heidi/brainseq/ea_only/`).

### 0.4.4 Revised order

1. ~~Re-anchor the locus list on `coloc.abf` PP4 rather than CLPP.~~ **Done 2026-09-11 at signal
   level** (§0.4.1), with the event audit and the long-read check re-run on the new nominations
   (§0.1).
2. ~~Map BrainSEQ switch and junction QTLs; repeat the anchoring and the per-gene modality contrast
   there.~~ **Done 2026-09-11/12** (§0.4.2): per-gene coloc favours abundance (5 vs 33,
   P = 4.3e-6); the junction layer was run, shelved, and dropped on 2026-09-12.
3. ~~SMR + HEIDI on both QTL sources.~~ **Done on GTEx and BrainSEQ EA-only, 2026-09-11/12**
   (§0.4.3). GTEx: **30/42** nominations have at least one tissue where the gene's pre-designated
   primary sQTL probe is SMR-significant with HEIDI not rejected (31/42 counting either modality;
   17/42 if the probe must be the exact intron that coloc headlined in that tissue); 9 have no
   instrument at 5e-8, 2 are HEIDI-rejected and 1 is instrumented but null. *(An earlier 29/42 in
   this plan matched none of these definitions and is superseded.)* SNCA/PD disagrees, and the
   disagreement is a HEIDI rejection rather than an absence of SMR signal (HEIDI rejects in all 8
   instrumented sQTL tissues, 6 of which clear the family threshold).
4. **Only then** decide how much wet-lab is still needed. Steps 1–3 may reduce Experiment A to
   a single confirmatory reporter assay on SNCA, or make it unnecessary for the resource claim.

---

## 1. Experiment A — functional validation of a switch

**Question.** Does an anchored isoform switch *do* something to the transcript, or is it a
correlated measurement?

### 1.1 The targeting trap

The 12 genetically anchored, splicing-led genes split into two disjoint groups, and the split runs
exactly opposite to what a naive target pick would assume.

| Gene | Traits | max CLPP | max anchored IF | Long-read confirmed | Group |
| --- | --- | --- | --- | --- | --- |
| **PGS1** | ALS | **0.093** | 0.818 | yes | **dual evidence** |
| **PPP6R2** | ALS, SCZ | **0.065** | 0.904 | yes | **dual evidence** |
| CDIP1 | SCZ | 0.028 | 0.306 | yes | dual evidence |
| DLG1 | SCZ | 0.026 | 0.400 | yes | dual evidence |
| TBC1D15 | PD | 0.017 | 0.894 | yes | dual evidence |
| PRRC2B | SCZ | 0.017 | 0.281 | yes | dual evidence |
| TPCN1 | AD | 0.011 | 0.894 | yes | dual evidence |
| GGNBP2 | ALS | 0.011 | 0.930 | yes | dual evidence |
| ARVCF | SCZ | 0.010 | 0.405 | yes | dual evidence |
| **CTSH** | AD | **0.386** | **0.004** | **no** | genetics only |
| **SNCA** | LBD, PD | 0.038 | **0.003** | yes (not at usable abundance) | genetics only |
| RTEL1 | SCZ | 0.010 | 0.545 | no | neither |

**The gene with by far the strongest colocalization — CTSH, CLPP 0.386, an order of magnitude
above every other gene — has an anchored isoform at 0.4% usage and fails orthogonal
confirmation.** SNCA is the same shape: the narratively attractive locus sits at 0.3% long-read
usage. Meanwhile the genes whose switches are unambiguously real (IF 0.28–0.93, all long-read
confirmed) carry CLPP of 0.010–0.093.

Picking targets on genetic evidence alone lands on the two genes that cannot be assayed. Picking
on expression alone lands on genes a reviewer will call genetically weak. **Validate the
intersection and say why.**

### 1.2 Recommended targets

- **Primary: PGS1** (CLPP 0.093, IF 0.818, ALS) and **PPP6R2** (CLPP 0.065, IF 0.904, ALS + SCZ).
  Best joint genetic and expression evidence; PPP6R2 is additionally cross-trait.
- **Secondary: DLG1** (CLPP 0.026, IF 0.400, SCZ) — a synaptic scaffold, so a functional readout
  has interpretable neuronal meaning.
- **Promoted to primary: SNCA.** *(Corrected 2026-09-10 — the earlier draft demoted SNCA on
  the ONT 0.29% figure, which is the number this project had already concluded was an assay
  limitation. Quoting it as a reason not to assay the gene was an error.)* The short-read
  junction test puts minor-form usage at **0.188 in DLPFC** and **0.234 in caudate**, with
  **99.5% of 222 donors** above the 5% threshold — a well-measured, well-used switch. Paired
  with signal-level `coloc.susie` PP4_sQTL **0.974** vs PP4_eQTL **0.113** for LBD in 11/13
  tissues (§0.4.1; the call sits just below 0.8 at p12 = 1e-6), SNCA has the best joint genetic
  and expression evidence of any locus here.

### 1.3 Assay, chosen from the consequence data

Because UTR remodeling is the dominant enriched consequence (1.279, 9/9 regions) and NMD is
*depleted*, the readout should be UTR-dependent function, not decay and not protein truncation:

1. **Isoform-resolved quantification** in an independent brain panel — targeted amplicon or
   Nanopore cDNA over the two competing isoforms. Establishes the switch exists outside BrainSEQ
   and GTEx. *This is the minimum viable experiment and could stand alone.*
2. **UTR reporter.** Clone the alternative 3′/5′ UTR of each target into a dual-luciferase or
   destabilised-GFP backbone; measure steady-state output and, with actinomycin-D chase,
   transcript half-life. Read out **relative** difference between the two UTRs, in a neuronal
   line (SH-SY5Y or iPSC-derived neurons).
3. **Allelic direction — use the existing phASER ASE data, not a new assay.**
   *(Added 2026-09-10 on co-author suggestion; supersedes the population-level allelic test the
   first draft proposed.)*

   `/ocean/projects/bio260021p/shared/resources/processed-data/ase-files/{caudate,dlpfc,hippocampus}/`
   holds phASER output for BrainSEQ in exactly the three aging regions — 578 GB, including
   `haplotypes`, `variant_connections`, `haplotypic_counts`, `site_ase`, `allelic_counts` and
   `gene_ae`. The phasing, which is the expensive part, is already done.

   **`gene_ae` is keyed on ENST, not ENSG** — one row per *transcript* per donor with
   `aCount`/`bCount`/`log2_aFC` and a `gw_phased` flag. That is allele-specific *isoform*
   quantification, which is the axis a switch lives on.

   The test: for a switch pair (T₁, T₂) in the same gene and the same genome-wide phase block,
   ask whether the haplotype ratio **differs between the two transcripts**. A cis variant that
   drives the switch shifts T₁ and T₂ in opposite directions on the same haplotype; a variant
   that only drives expression shifts both together. **Each donor is its own control**, so the
   within-donor design removes population-structure confounding and controls much donor-level
   trans and environmental variation — a cleaner design than any population-level allelic
   association. It does not remove everything, and three things survive it: bulk cell composition
   does **not** cancel if the allele-specific effects themselves differ by cell type; residual
   reference-mapping bias remains; and random allelic imbalance remains. WASP filtering was
   applied when the ASE data were generated (`qc/wasp_qc.tsv`, per-sample `wasp_vcfs/`), which
   mitigates the mapping bias rather than eliminating it.

   Three caveats, all real:
   - **Per-donor depth is thin.** Allelic counts need reads over heterozygous sites; spot-checking
     SNCA transcripts in one caudate donor gives `totalCount` of 0–2 per transcript. Power comes
     from aggregating heterozygotes across n ≈ 450 donors per region (beta-binomial or
     phASER-POP style), not from any single sample.
   - **Use `gw_phased = 1` rows only.** Cross-transcript comparison requires haplotype A to mean
     the same haplotype for both transcripts.
   - **For a specific junction, prefer junction-spanning allelic reads** (`allelic_counts` /
     `site_ase` intersected with the junction) over whole-transcript counts, which inherit the
     same read-assignment ambiguity that makes isoform quantification hard in the first place.

   Do **not** use `gene_ae` aggregated to the gene, or the gene-level `aFC`: allelic imbalance in
   total gene output is the *eQTL* axis and would be the wrong test for a switch.

### 1.4 What the junction test did and did not establish

The SNCA/CTSH junction analysis used `min(median PSI, 1 − median PSI)` — deliberately
convention-free, because PSI orientation is undocumented for the LIBD tables. It asks whether
**both forms are actually used**. That is an existence-and-usage result, and it is conclusive:
SNCA's anchored switch is real, well-measured and abundant; CTSH's is not.

It is **not** a genotype-association test and **not** a functional test. So for SNCA the open
question is no longer "is the switch real" — it is "does the risk allele drive it, and does the
resulting UTR change do anything". Step 3 below (allelic direction) and the reporter assay are
what remain, and they are a smaller ask than the earlier draft implied.

### 1.5 Pre-registered decision rule

Write this before running, in the manner already used for SNCA/CTSH:

> A target is **functionally validated** if (a) both isoforms are detected in the independent
> panel at ≥ 5% usage, **and** (b) the two UTRs differ in reporter output or half-life at
> p < 0.05 with the direction stated in advance. A target is **not** validated if either fails;
> report it and keep it off the main figure.

Two of three targets validating supports a main-figure mechanism panel. Zero of three means the
paper stays a set-level correlation and should be framed that way — which is publishable, just
not at the tier this experiment is meant to reach.

---

## 2. Future grant aim — snRNA-seq and the composition question

> **Scoped out of this manuscript, 2026-09-10.** At the n = 40 that would make it a real study,
> single-nucleus long read is too expensive to serve as functional validation for this paper,
> and — per the power analysis in §2.3 — an aging-correlation design at that n is *worse
> powered than the bulk experiment it would adjudicate*. It is a good standalone aim and a poor
> validation step.
>
> **For this manuscript:** disclose the composition question as a stated limitation, with the
> decisive experiment named in the Discussion. That is already the agreed handling for the SCZD
> arm (P1 round-1 item 6).
>
> **As a grant aim:** the design below is the one to propose — n = 40 donors sampled at the
> extremes of the cell-type proportion distribution (20 per group), short-read cell-type labels
> paired with long-read isoform quantification on the same nuclei, powered to ΔPSI ≈ 0.04–0.13
> depending on between-donor variance. The sections that follow are written to be liftable into
> an aims page.

**Question.** Do the switches persist *within* matched cell types, or are they a proportion shift?

### 2.1 What is at stake (preliminary data for the aim)

| Contrast | comp_unique base → adjusted | Retained |
| --- | --- | --- |
| BrainSEQ SCZD caudate | 34 → **2** | **6%** |
| BrainSEQ aging caudate | 43 → 17 | 40% |
| BrainSEQ aging DLPFC | 8 → 15 | 188% (adjustment *increased* it) |
| BrainSEQ aging hippocampus | 0 → 0 | undefined — contributes nothing, must not be counted as passing |
| GTEx frontal cortex BA9 | 531 → **0** | **0%** |
| GTEx cortex | 438 → **0** | **0%** |
| GTEx anterior cingulate BA24 | 545 → 88 | 16% |
| GTEx hippocampus | 61 → 21 | 34% |
| GTEx caudate | 32 → 10 | 31% |

The disease arm is already conceded as composition-entangled. The sharper problem is that **both
GTEx cortical regions collapse to exactly zero** while limbic and striatal regions retain 16–34%.
Bulk deconvolution cannot say whether that is real regional biology or over-adjustment, because
proportion and per-cell composition are confounded in the same measurement by construction.

### 2.2 The technical constraint that decides the design

**Standard 3′-biased snRNA-seq cannot measure isoform usage.** A 10x 3′ library reads the last
few hundred bases of a transcript; most of the switches here are UTR or internal-CDS remodeling
events that a 3′ tag cannot distinguish. Running conventional snRNA-seq and hoping to recover
isoform ratios would produce a null that means nothing.

Viable options, in order of preference:

1. **10x 3′ + long-read on the same libraries** (PacBio Kinnex / MAS-ISO-seq, or ONT). Cell-type
   labels from the short reads, isoform quantification from the long reads on the identical
   nuclei. This is the design that answers the question.
2. **Targeted capture** of the specific junctions for the anchored gene set. Far cheaper, but
   only tests the genes already nominated — it cannot re-derive the layer, so it settles the
   named loci and not the composition question in general.
3. **Single-nucleus long-read alone.** Cleanest, most expensive, lowest throughput per donor.

This project has an existing 3′ bias sensitivity finding; the same physics applies here and
should be cited in the methods rather than rediscovered.

### 2.3 Power — and why the obvious design fails

Detectable ΔPSI between two groups (two-sided α = 0.05, 80% power), by donors per group and
between-donor PSI standard deviation:

| n / group | SD = 0.05 | SD = 0.08 | SD = 0.10 | SD = 0.15 |
| --- | --- | --- | --- | --- |
| 10 | 0.063 | 0.100 | 0.125 | 0.188 |
| 20 | 0.044 | 0.071 | 0.089 | 0.133 |
| 30 | 0.036 | 0.058 | 0.072 | 0.109 |
| **40 (20/group — the proposed grant n)** | **0.031** | **0.050** | **0.063** | **0.094** |
| 50 | 0.028 | 0.045 | 0.056 | 0.084 |

Age as a **continuous** predictor is the trap:

| Donors | Detectable \|r\| | R² |
| --- | --- | --- |
| 20 | 0.591 | 0.349 |
| 40 | 0.431 | 0.185 |
| 60 | 0.355 | 0.126 |
| 100 | 0.277 | 0.077 |
| **238 (the bulk BrainSEQ caudate n)** | **0.181** | 0.033 |

**A 40-donor snRNA-seq study detects \|r\| ≥ 0.43; the bulk study that generated the hypothesis
detects \|r\| ≥ 0.18.** Any realistically sized single-nucleus aging study is *worse powered than
the experiment it is meant to adjudicate*, so a null result would be uninterpretable. Do not run
this as an age-correlation study.

**Design for the aim:** frame it as a *composition* test, not an aging test. Sample donors at the
extremes of the cell-type proportion distribution (or case/control for the SCZD arm), and ask
whether the switch persists within matched cell types across that contrast. Extreme sampling buys
back effect size that continuous age cannot. n = 20–30 per group detects ΔPSI ≈ 0.04–0.09 at
plausible SDs, which brackets the bulk effects.

### 2.4 Pre-registered decision rule

> The disease/aging switch layer is **cell-intrinsic** if the anchored switches show ΔPSI in the
> stated direction within at least one matched cell type at FDR < 0.05, with the cell-type
> proportion shift explicitly modelled. It is **compositional** if the switches are present
> across the proportion contrast in pseudobulk but vanish within every matched cell type at
> adequate power. Report the power achieved per cell type; a cell type with fewer than 20 donors
> contributing nuclei is reported as untested, not as negative.

---

## 3. Recommendation

### 3.1 If you run one thing

**Experiment 0 (§0.4) — the computational genetics, not the bench.** *(Updated 2026-09-11.)* The
re-anchoring is done, at signal level: it replaces a CLPP list topping out at 0.093 with 42
nominations at PP4_sQTL ≥ 0.8, 14 of them sQTL-preferential, and an event audit that recovers no
known mechanism — UNC13A is context-distinct from its cryptic exon, and PICALM is not signal-level
on the primary grid. What remains is mapping BrainSEQ switch QTL — same-tissue genetic anchoring in
the discovery cohort, not independent replication — which costs compute and attacks both of the
counter-lines a reviewer will press hardest. Neither needs a sample.

### 3.2 If you run two

Add **Experiment A step 3** — allelic direction from the existing BrainSEQ phASER ASE data
(§1.2). Transcript-level haplotypic counts already exist for caudate, DLPFC and hippocampus, so
this is a re-analysis rather than an experiment, and each donor acts as its own control. It is
the step that ties a switch to its *risk allele* rather than to aging. The reporter assay and
the snRNA-seq arm come after, sized by what remains open.

### 3.3 What is deliberately not in scope

**Single-nucleus long read is a grant aim, not a validation step** (§2). At n = 40 it costs more
than the rest of this page combined and, as an aging-correlation design, is worse powered than
the bulk study it would adjudicate. The composition question is disclosed as a limitation for
this manuscript and proposed as future work — which is already the agreed handling for the SCZD
arm. §2 is kept here in aims-page form so the preliminary data, the technical constraint and the
power table travel together when it is written up.

### 3.4 Correction to the first draft

This plan originally led with the bench and demoted SNCA on the ONT 0.29% usage figure. Both
were wrong. The 0.29% is the assay limitation this project had already diagnosed; short read
puts the same switch at 18.8% in 99.5% of donors. And the computational genetics — coloc.abf,
BrainSEQ junction QTLs, SMR — were not exhausted before proposing wet-lab work. They should be.

### 3.5 Framing fixes that cost nothing and should happen regardless

**All three applied 2026-09-12.**

1. ~~`MANUSCRIPT_PLAN.md` §12 still lists **"Splicing-QTL specifically anchor the GO-invisible
   switch layer"** as *Established* and as the Fig 3 headline.~~ **Fixed.** The row now reads
   "The phenotype-associated switch layer is genetically anchored to splicing,
   method-specifically", with the non-localisation (1.068, p = 0.077 vs GO-visible 1.084) as its
   caveat, so the claims table and `FIGURE_ORDERING.md` no longer disagree about the headline.
2. ~~State the S-LDSC splicing arm as it is: one of six nominal, uncorrected.~~ **Fixed.** The
   Table 3 legend used to cite partitioned heritability as set-level support; it now says the
   S-LDSC splicing arm does not support the claim and gives the 1-of-6 / AD p = 0.0355 /
   uncorrected numbers inline. The fix went into the generator,
   `manuscript/_h/assemble_main_tables.py`, not the emitted `.md` — regenerate the tables.
3. ~~`reports/pi/00_OVERVIEW.md` §8 is stale — 89% should be 94%, and SNCA is no longer
   "uncertain".~~ **Fixed** — though §8 itself had already been corrected; the stale 89% /
   25-of-130 numbers were surviving in the stage table and the §03 narrative instead, and both
   now carry the post-2026-09-09 split-half re-fit (250/266, 23/130, perm P = 0.034).

---

## Appendix — provenance

| Number | Source |
| --- | --- |
| 1.111 / 1.112 contrasts, baseline nulls | `05_genetic_anchoring/_m/qtl_anchoring_meta/QTL_ANCHORING_META.md` (re-run 2026-09-09 after the matched-baseline re-fit) |
| 0.453 / 0.252, 0.600 / 0.304 | `06_switch_mechanism/_m/switch_orthogonal_confirm/anchored_summary.json` |
| 12-gene anchored table, CLPP, IF | `06_switch_mechanism/_m/switch_orthogonal_confirm/anchored_gene_confirmation.parquet` |
| UTR 1.279 / CDS 1.044 / NMD depleted | `06_switch_mechanism/_m/switch_consequence_meta.parquet`, `analysis_class = aging, stratum = all` |
| 250/266, 23/130 p = 0.034 | `03_module_trust/_m/stability/module_trust/`, re-run 2026-09-09 on corrected split halves |
| SNCA 0.189 / 0.234 | `06_switch_mechanism/_m/junction_coloc_confirm/` |
| Composition table | `04_module_characterization/_m/COMPOSITION_ADJUSTMENT_SUMMARY.md` |
| Per-gene coloc Fisher p = 0.88 | `05_genetic_anchoring/_m/coloc_modality_contrast/` |
| S-LDSC tau p-values | `05_genetic_anchoring/_m/ldsc/LDSC_SUMMARY.md` |
| BrainSEQ phASER ASE (transcript-level `gene_ae`, haplotypic counts, phasing) | `/ocean/projects/bio260021p/shared/resources/processed-data/ase-files/{caudate,dlpfc,hippocampus}/` — 578 GB; see `ASE_GENERATION.md` |
| Power tables | Computed 2026-09-09; two-sample t and Fisher-z, α = 0.05 two-sided, power = 0.80 |
| *(superseded)* coloc.abf PP4 (40 strong hits, 13 sQTL-preferential) | `05_genetic_anchoring/_m/coloc_modality_contrast/genes.parquet`, switch arm — 1,647 cells / 1,156 genes over GTEx v11 all-pairs |
| Signal-level nominations (42, 14 sQTL-preferential), descriptors, paired contrast | `05_genetic_anchoring/_m/coloc_signal_susie/all_introns/{genes,cells_hierarchy}.parquet`; contrast in `05_genetic_anchoring/_m/coloc_signal_susie/contrast.parquet` (representative arm) — regenerated 2026-09-11 |
| Event-audit tiers, curated-event testability | `05_genetic_anchoring/_m/locus_event_audit/susie_all_introns/`; PICALM sensitivity tier in `susie_all_introns__max_snps_30000/` |
| SMR + HEIDI agreement (30/42 primary sQTL probe; 31/42 either modality; 17/42 exact coloc intron; SNCA/PD disagreement) | `05_genetic_anchoring/_m/smr_heidi/gtex/{smr_results.parquet,SMR_HEIDI.md}` |
| Signal-level long-read 0.565 / 0.238, 0.672 / 0.298 | `06_switch_mechanism/_m/switch_orthogonal_confirm/signal_coloc/anchored_summary.json` |
| SNCA minor-form usage 0.188, 99.5% of 222 donors | `06_switch_mechanism/_m/junction_coloc_confirm/junction_confirm.parquet` |
| BrainSEQ genotypes | `inputs/raw/brainseq/genotypes/TOPMed_LIBD.{pgen,pvar,psam}` (TOPMed-imputed, PLINK2) |
| GTEx v11 sQTL SuSiE fine-mapping | `/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11/GTEx_Analysis_v11_sQTL/*.SuSiE_summary.parquet` |
