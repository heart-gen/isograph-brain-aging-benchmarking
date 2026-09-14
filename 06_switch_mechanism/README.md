# 06 — Switch mechanism

**Question:** are the switches real measurements rather than model artifacts, what do
they do to the protein, and do the results survive sensitivity analysis?

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01–02 | `switch_consequence`, `switch_consequence_meta` | Structural consequence of the switch axis under a within-gene permutation null; cross-region meta |
| 03–04 | `validate_switch_splicing_{brainseq,gtex}` | Orthogonal PSI / junction validation of switch pairs |
| 05 | `switch_orthogonal_confirm` | Anchored + global-null confirmation against the compositional-closure baseline. `--events signal` scores the `coloc.susie` nominations (from `coloc_isoform_events --layer signal`) into `signal_coloc/`; the default CLPP arm and its S-real-8 / S14 outputs are unchanged |
| 06 | `isa_concordance` | Concordance with satuRn / ISA differential transcript usage |
| 07 | `longread_switch_confirm` | ONT DLPFC long-read confirmation (Aguzzoli-Heberle 2024) |
| 08–09 | `download_clinical`, `clinical_consequence` | gnomAD LOEUF constraint + ClinVar pathogenic density |
| 10–11 | `scz_confound_sensitivity`, `switch_feature_sensitivity` | Medication/toxicology/smoking availability audit; five feature-construction axes (one region per call; `--aggregate` writes the cross-region report and runs the quantification axis once, beside a same-quantifier reference) |
| 12 | `junction_coloc_confirm` | Short-read BrainSEQ junction confirmation of the anchored SNCA/CTSH switches |
| 13 | `switch_feature_sensitivity_gtex` | Step 11 over all 13 GTEx regions (array). GTEx's published setting is the **unfiltered** transcript matrix, so its expression axis runs from no filter upward |
| 14 | `switch_feature_refit` | Full production refit per preprocessing setting (BrainSEQ caudate, 12 fits incl. the published-setting noise floor) — lets the network change, which step 11 holds fixed |

## Key results

- Consequence is **productive remodeling, not decay**: of 9 structural classes only
  UTR-remodeled (1.27×, 10/10 regions) and CDS-remodeled (1.04×, 10/10) are enriched;
  NMD routing, biotype switch and coding-status loss are depleted. Identical in
  GO-invisible and GO-visible modules.
- Switch genes are more LoF-constrained than genome-wide (median LOEUF 0.72 vs 0.94),
  while switched exons carry *lower* ClinVar P/LP density (ratio 0.18) — expected
  alternative-exon biology.
- **The long-read switch-like rate must never be quoted without its null, and not on the
  all-pairs denominator.** Step 07 reports `pair_switch_like_rate` = **0.084** (5,173/61,302
  in BA9; 1,341/15,878 in cortex) — that denominator is *every* prespecified switch pair,
  including the ~73% never detected in n = 12 ONT samples, so it is a detection rate times a
  concordance rate and cannot be compared to anything. Step 05 supplies the comparison the
  number was missing, on the **detected-pair** denominator and against pairs drawn from the
  *same genes*, scored on the *same samples* by the *same code*, matched on abundance decile:

  | Set | switch-like | matched null | P |
  |---|---|---|---|
  | All IsoGraph switch pairs (18,268 detected) | 0.647 | 0.639 | 0.022 |
  | Genetically anchored pairs (53 detected) | 0.453 | 0.252 | 5e-4 |
  | …restricted to usable abundance (30) | 0.600 | 0.304 | 5e-4 |

  The global arm barely clears its null because **compositional closure** puts the null at
  0.639 before any biology: within-gene isoform fractions sum to one, so any two isoforms of
  a gene are negatively correlated by construction. The anchored subset is the result; the
  global switch-like rate is not. See `_m/switch_orthogonal_confirm/GLOBAL_NULL.md`.

- **SNCA's anchored switch is confirmed in short read; CTSH's is not** (step 12). On the
  exact contrast Fig 4A draws — the anchored proximal first exon against the canonical
  distal one — minor-form usage is **0.189 in DLPFC (n = 222)** and **0.234 in caudate
  (n = 238)**, far above the pre-registered 0.05 threshold. The ONT long-read check (step
  07) had put the anchored isoform at **0.29%** and failed it; short read measures the
  junction the sQTL actually tags rather than a whole-transcript proxy, at ~20x the n, so
  the long-read result is an assay limitation rather than a refutation. **CTSH reaches only
  0.016–0.020** in hippocampus and therefore stays off any main figure.

  Two things about step 12 are worth knowing before reusing this machinery. First, the
  statistic is `min(median PSI, 1 - median PSI)`: PSI orientation is undocumented for the
  LIBD tables, and this is invariant to it. Second, the analysis n is the **aging bundle**
  (238 adult controls in hippocampus), not the 452 sample columns in the PSI file — the
  same sample definition as every other aging analysis here.

## Two sensitivity results that must be reported as stated

- **rRNA rate is not a confounder.** It has no marginal association with diagnosis
  (d = −0.009, p = 0.93); its association appears only after conditioning on the
  published covariates (partial r = +0.108), which is the collider signature. The
  modules *most* correlated with it survive and the least correlated fall — the opposite
  of technical-artifact removal. Do not add it to the inference model on confounding
  grounds. Report as a sensitivity.
- **Per-gene switch-age effects do not transfer across cohorts at all, but do across tissues
  within a cohort** (steps 11/13 `--aggregate`, 2026-09-12). Matched regions across cohorts
  (Salmon vs RSEM) give Pearson 0.007 caudate, −0.002 hippocampus, 0.0004 DLPFC/BA9 (sign
  concordance 0.50–0.51). The new same-quantifier reference — different tissues within one
  cohort — gives 0.31–0.51 in GTEx and 0.18–0.23 in BrainSEQ (sign 0.56–0.66). Those pairs share
  most donors, which inflates them, and no same-quantifier cross-cohort pair exists, so quantifier
  and cohort remain confounded; but the cross-cohort figure is not merely low, it is zero, far below
  what tissue differences alone produce. Module-level replication (stage 03) survives this only as
  a count against a matching null, never as correlated per-gene effects.
- **The switch representation is robust to the pseudocount everywhere, and to transcript
  filtering only where the published matrix is already filtered.** The harness holds the module
  partition fixed and ran on BrainSEQ caudate plus all 13 GTEx regions. Pseudocount 0.1–2.0:
  minimum module age-effect correlation with the published effects 0.99 (BrainSEQ) and 0.81–0.98
  (GTEx), at most one sign flip among published-significant modules. Expression filter and
  minor-isoform threshold: 0.89–0.90 with no flips in BrainSEQ, whose fit already applies a count
  filter, but 0.07–0.69 in GTEx at the strictest setting, with up to 12 flips (BA9). GTEx is fit on
  the **unfiltered** transcript matrix, so its filter axis removes up to ~47% of switch genes and
  rewrites the coordinates of those left (median per-gene |r| 0.25 at count > 10 in ≥ 90%).
  The choice of expression filter is therefore a real analytic degree of freedom for GTEx and must
  be stated as fixed a priori.
- **Letting the network refit changes the picture: module identity is preprocessing-sensitive,
  the age signal less so** (step 14, BrainSEQ caudate, 12 production fits). Refitting the
  published setting reproduces the committed partition (ARI 0.97, 19/20 age associations), so
  refit noise is small. Every perturbed setting moves the partition far more: ARI 0.66–0.77 for
  pseudocount, 0.44–0.60 for the expression filter, 0.65 / 0.43 / 0.29 for minor-isoform usage
  ≥ 0.01 / 0.05 / 0.10. Of the 20 published FDR-significant module–age associations, 15–20 recur
  with the same sign in the best-matching refit module, but that match is loose — median best
  Jaccard 0.02–0.27 — so the aging signal reappears in substantially redrawn modules rather than in
  the same modules. The fixed-partition harness above therefore overstates stability for anything
  that names a specific module; the defensible statement is at the level of the switch layer's age
  association, not module membership. `_m/switch_feature_sensitivity/refit/brainseq/caudate/REFIT_SENSITIVITY.md`.

**CLIs:** `isograph_benchmark/real_data/{switch_consequence,switch_consequence_meta,validate_switch_splicing,switch_orthogonal_confirm,isa_concordance,longread_switch_confirm,clinical_consequence,clinical_consequence_meta,scz_confound_sensitivity,switch_feature_sensitivity,switch_feature_refit,junction_coloc_confirm}.py`.

## Display items

S-real-5 `figSwitchConsequence`, S-real-7 `figClinicalConsequence`.


## phASER allelic direction — screened, and where the follow-on must run

`ase_switch_direction.py --stage screen` asks whether the existing BrainSEQ phASER release
can support an allele-specific test of an isoform switch. It cannot, and the reason is not
the one the QC tables suggest.

**Depth is not the blocker.** The phASER QC tables report a median `gene_ae` depth of 1-2
reads, but that is the whole-transcript-interval aggregation. Restricted to
isoform-discriminating exonic sequence, the donors that carry signal carry plenty: SNCA has
243 donors at median 77 reads. Of 37 caudate switch pairs with any informative site, 28
have enough donors at adequate depth.

**The blocker is structural.** Only 5 of 37 pairs are informative on BOTH isoforms; 25 fail
on that alone with donors and depth to spare. A pair where one isoform is a truncation of
the other has no unique exonic sequence on that side, so there is nothing to contrast
however deeply it is sequenced. Three pairs (two genes) clear the screen — too few to carry
a claim.

**The rescue, and where it has to run.** Exonic-unique segments are a one-sided
discriminator; a discriminating *junction* is two-sided by construction. Counting reads
that cross a discriminating junction while carrying a phased heterozygous site would test
many of those 25 pairs. That needs read-level alignments and, specifically, the
**WASP-tagged ASE-grade BAMs** — reference mapping bias at heterozygous sites is the
dominant confound in a within-donor allelic contrast, and `ASE_GENERATION.md` is explicit
that a non-WASP alignment is not ASE-final.

Those BAMs are **not on Bridges-2**. The validated manifest's `bam_path` and `source_file`
both point at Northwestern's **Quest**, `/gpfs/projects/b1042/HEART-GeN-Lab/ase-processing/`,
and neither PSC project share holds a BAM or CRAM. CRAM is not the obstacle — samtools,
pysam and GATK all read it and the GRCh38 reference is available — the obstacle is that the
reads and their WASP tags exist only on Quest. **Run the junction-level arm on Quest
(b1042), where the alignments and the storage are.**
