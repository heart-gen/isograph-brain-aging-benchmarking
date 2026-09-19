# 06 — Switch mechanism

**Question:** are the switches real measurements rather than model artifacts, what do
they do to the protein, and do the results survive sensitivity analysis?

## Order

Run the whole stage with `bash 06_switch_mechanism/_h/run_stage.sh` (add `--dry-run` to print
the plan). The leading number of a wrapper is its tier; steps in one tier run in parallel.

| Step | Wrapper | Waits on | Produces |
|---|---|---|---|
| 01a | `switch_consequence` | stage 03 | Structural consequence of the switch axis under a within-gene permutation null |
| 01b–01c | `validate_switch_splicing_{brainseq,gtex}` | stage 03; `05_genetic_anchoring/_h/05c` events; GTEx junction usage | Orthogonal PSI / junction validation of switch pairs |
| 01d | `isa_concordance` | stage 02 | Concordance with satuRn / ISA differential transcript usage |
| 01e | `longread_switch_confirm` | stage 02 | ONT DLPFC long-read confirmation (Aguzzoli-Heberle 2024) |
| 01f | `download_clinical` | — (**login node**) | ClinVar + gnomAD v4.1 references into `inputs/raw/clinical/`; skips files already present |
| 01g | `scz_confound_sensitivity` | stage 02 | Medication/toxicology/smoking availability audit |
| 01h | `switch_feature_sensitivity` | stage 02 | Five feature-construction axes on BrainSEQ caudate, partition held fixed |
| 01i | `switch_feature_sensitivity_gtex` | stage 02 | 01h over all 13 GTEx regions (array). GTEx's published setting is the **unfiltered** transcript matrix, so its expression axis runs from no filter upward |
| 01j | `switch_feature_refit` | stage 02 | Full production refit per preprocessing setting (BrainSEQ caudate, 12 fits incl. the published-setting noise floor) — lets the network change, which 01h holds fixed |
| 01k | `junction_coloc_confirm` | `05_genetic_anchoring/_h/05c` events | Short-read BrainSEQ junction confirmation of the anchored SNCA/CTSH switches |
| 02a | `switch_consequence_meta` | 01a | Cross-region consequence meta |
| 02b | `switch_orthogonal_confirm` | 01e; stage-05 coloc events | `--mode anchored` then `--mode global-null` against the compositional-closure baseline, then `--events signal` for the `coloc.susie` nominations (from `05_genetic_anchoring/_h/08a`) into `signal_coloc/` — one after another, they share a directory. The CLPP arm backs S-real-8 / S14 |
| 02c | `clinical_consequence` | 01a, 01f | gnomAD LOEUF constraint + ClinVar pathogenic density, per region |
| 02d | `switch_feature_sensitivity_aggregate` | 01h, 01i | Cross-region report; runs the quantification axis once, beside a same-quantifier reference |
| 02e | `switch_feature_refit_aggregate` | 01j | `refit/brainseq/caudate/REFIT_SENSITIVITY.md` |
| 03a | `clinical_consequence_meta` | 02c | Cross-region clinical-consequence rollup |

## Key results

- Consequence is **productive remodeling, not decay**: of 9 structural classes only
  UTR-remodeled (1.27×, 10/10 regions) and CDS-remodeled (1.04×, 10/10) are enriched;
  NMD routing, biotype switch and coding-status loss are depleted. Identical in
  GO-invisible and GO-visible modules.
- Switch genes are more LoF-constrained than genome-wide (median LOEUF 0.72 vs 0.94),
  while switched exons carry *lower* ClinVar P/LP density (ratio 0.18) — expected
  alternative-exon biology.
- **The long-read switch-like rate must never be quoted without its null, and not on the
  all-pairs denominator.** Step 01e reports `pair_switch_like_rate` = **0.084** (5,173/61,302
  in BA9; 1,341/15,878 in cortex) — that denominator is *every* prespecified switch pair,
  including the ~73% never detected in n = 12 ONT samples, so it is a detection rate times a
  concordance rate and cannot be compared to anything. Step 02b supplies the comparison the
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

- **SNCA's anchored switch is confirmed in short read; CTSH's is not** (step 01k). On the
  exact contrast Fig 4A draws — the anchored proximal first exon against the canonical
  distal one — minor-form usage is **0.189 in DLPFC (n = 222)** and **0.234 in caudate
  (n = 238)**, far above the pre-registered 0.05 threshold. The ONT long-read check (step
  01e) had put the anchored isoform at **0.29%** and failed it; short read measures the
  junction the sQTL actually tags rather than a whole-transcript proxy, at ~20x the n, so
  the long-read result is an assay limitation rather than a refutation. **CTSH reaches only
  0.016–0.020** in hippocampus and therefore stays off any main figure.

  Two things about step 01k are worth knowing before reusing this machinery. First, the
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
  within a cohort** (steps 01h/01i, aggregated by 02d, 2026-09-12). Matched regions across cohorts
  (Salmon vs RSEM) give Pearson 0.007 caudate, −0.002 hippocampus, 0.0004 DLPFC/BA9 (sign
  concordance 0.50–0.51). The new same-quantifier reference — different tissues within one
  cohort — gives 0.31–0.51 in GTEx and 0.18–0.23 in BrainSEQ (sign 0.56–0.66). Those pairs share
  most donors, which inflates them, and no same-quantifier cross-cohort pair exists, so quantifier
  and cohort remain confounded; but the cross-cohort figure is not merely low, it is zero, far below
  what tissue differences alone produce. Module-level replication (stage 04) survives this only as
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
  the age signal less so** (step 01j, BrainSEQ caudate, 12 production fits). Refitting the
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

## Allele-aware junction recount — the two-sided rescue, on Quest (PI item 10a)

`ase_junction_switch.py` counts fragments that cross an isoform-**specific** splice
junction *and* carry a phased heterozygous site, from the WASP-tagged STAR BAMs. A junction
is two-sided by construction, so the pairs the exonic screen lost for want of a second side
become testable. The pre-registered gate is unchanged and imported from
`ase_switch_direction.py`: per pair, ≥ 30 donors with ≥ 5 informative fragments,
informative on both isoforms. The arm proceeds only if ≥ 30 pairs of the **gate family**
pass (pairs in genes with an all_samples DLPFC switch-QTL, q < 0.05). Otherwise it is
reported as the pre-specified negative result.

**Runs on Quest only; nothing bulk moves from Bridges-2.** The BAMs have moved from the old
`ase-processing` path to `/projects/b1042/HEART-GeN-Lab/brainseq_data/bsp2_dlpfc/star-wasp/`
(498 DLPFC BAMs, read in place). The phASER release, including the per-sample VCFs with the
genome-wide phase `PW`, is at `/projects/b1213/resources/processed-data/ase-files/`. The
population-phased genotypes are at `/projects/b1213/resources/processed-data/genotypes/ase/`.
All paths are in `configs/data_sources.yaml` under `quest:`. The switch pairs and QTL
leads come through git-LFS. Only `junction_allelic_counts.parquet` and the report travel back.

| Step | Wrapper | Writes (`_m/ase_junction_switch/<region>/`) |
|---|---|---|
| targets | `01l.ase_junction_targets.sh` | `pairs.parquet`, `junctions.parquet`, `lead_genotypes.parquet`, `regions.bed`, `samples.tsv` |
| count | `02f.ase_junction_count.sh` (array, one BAM per task) | `counts/<sample>.{pair_counts,pair_sites,site_counts}.tsv.gz`, `.qc.json` (gitignored) |
| concordance | `python -m …ase_junction_switch --stage concordance --sample R#####` | `counts/<sample>.concordance.json` |
| screen | `03b.ase_junction_screen.sh` | `junction_allelic_counts.parquet`, `pair_feasibility.parquet`, `screen_summary.json`, `ASE_JUNCTION_SCREEN.md` |
| test (step 4) | `04a.ase_junction_allelic.sh` (`ase_junction_allelic.py`) | `allelic_test.parquet`, `allelic_donor_counts.parquet`, `allelic_summary.json`, `ASE_JUNCTION_ALLELIC.md` |

Run all four for one region with
`bash 06_switch_mechanism/_h/ase_junction_quest.sh --region <dlpfc|caudate|hippocampus>`.
It sizes the count array from that region's manifest. `targets` drops any BAM that fails
`samtools quickcheck`, because the BAMs are staged onto Quest one region at a time
(DLPFC `bsp2_dlpfc/`, caudate `bsp3/`, hippocampus `bsp2_hippo/`). The Bridges-2 DAG lists
the four as manual steps and never submits them.

A fragment counts only if no mate fails WASP (`vW` ≠ 1), every junction it carries inside
the pair's span is an intron of the isoform it is assigned to, and all its phased sites
agree on the haplotype.

**Step 4, the allelic test (`ase_junction_allelic.py`), runs only where the gate passed**
(it exits cleanly as the pre-specified negative otherwise). It reads the lead-ALT haplotype
straight off each donor's phased genotype at the all_samples lead. That works because `PW`
is anchored to the same population phasing: at the DLPFC gate-family leads, `PW` equals the
lead `GT` in 243/243 heterozygous calls checked. Per pair, a beta-binomial GLMM on
donor × haplotype units (T1 of T1 + T2 fragments, a random donor intercept, LRT on the ALT
effect `beta`) gives the within-donor log odds of T1 on ALT vs REF. A random intercept, not a
fixed one: with two units per donor a fixed intercept biases `beta` away from zero. Pairs
need ≥ 10 lead-heterozygous donors informative on both haplotypes. BH runs over the fitted
gate-family pairs. Donors homozygous at the lead are the built-in null (same model, arbitrary
haplotype). A model-free stratified score test and the between-donor dosage-vs-T1-fraction
Spearman, from the same reads, are reported beside it. `beta` is ALT-oriented;
`--risk-alleles` (a few-KB rsID → risk-allele TSV from Bridges-2) re-orients it. Each pair also carries its IsoGraph module, the gene's module role, the module's age association and its module polarity r(T1) − r(T2). The module-direction check asks whether `beta` tracks that polarity across a gene's pairs, against a within-gene permutation null: does the lead move the isoforms along the module's own switch axis? Caveats: composition cancels only if allelic effects do not differ by cell type, and
some mapping bias survives WASP.

### After the git transfer: mapping `beta` to disease direction (risk-allele table)

`beta` is oriented to the switch-QTL lead's **ALT** allele, which says nothing about disease.
Turning it into "the risk allele shifts the isoforms toward T1/T2, and toward or away from the
module's direction" needs each lead's GWAS risk allele. The GWAS summary statistics live on
Bridges-2 (ALS is not on Quest), so this step runs there after the tables come over.

1. **Quest → Bridges-2.** Commit `_m/ase_junction_switch/<region>/allelic_test.parquet` (plus
   `allelic_donor_counts.parquet`, `allelic_summary.json`, `ASE_JUNCTION_ALLELIC.md`; parquet is
   git-LFS), push, then on Bridges-2 `git pull && git lfs pull`. `allelic_test.parquet`
   already carries what orientation needs: `variant_id_all` (lead rsID), `lead_ref`,
   `lead_alt`, `palindromic`, `coloc_traits`, `beta`, `module_polarity`, `module_age_trait`,
   `module_age_effect`.
2. **Build the risk-allele table on Bridges-2.** For the coloc-nominated rows, take each
   (trait in `coloc_traits`, `variant_id_all`) and call `coloc_direction._gwas_risk(trait,
   rsids, tmp_dir)`, which returns `rsid, risk_allele, gwas_other_allele, risk_beta, gwas_p`.
   Write one TSV per trait (for example `_m/ase_junction_switch/risk_alleles.<trait>.tsv`).
   A lead can colocalize with several traits whose risk alleles differ, so do not collapse
   traits into one file.
3. **Orient.** Either compute it on Bridges-2 directly from `allelic_test.parquet`, or commit
   the TSV and re-run `04a` on Quest with `--risk-alleles <tsv>`, one trait per run (the flag
   keeps one risk allele per rsID). The rule is the same either way:
   `beta_risk = +beta` if `risk_allele == lead_alt`, `-beta` if `risk_allele == lead_ref`,
   otherwise missing. Then `risk_along_module = beta_risk * sign(module_polarity)`: > 0 means
   the risk allele shifts the pair's isoforms the way the module score rises. Read it against
   the module's age association only where `module_age_trait == Age_linear`. An `Age_spline`
   effect is a slope at one age point (p25), not an overall direction.
4. **Check before believing a sign:**
   - **Allele match.** Drop any lead whose GWAS alleles are neither `lead_ref` nor
     `lead_alt`, and any `palindromic` (A/T, C/G) lead unless the GWAS frequency confirms the
     strand. Report how many were dropped.
   - **Same variant.** The risk allele must be looked up at the switch-QTL lead
     (`variant_id_all`, all_samples arm), not at coloc's best SNP. If the GWAS lacks the lead,
     use an EUR LD proxy and carry the allele across through the phased haplotype, not by
     frequency; say so for each gene.
   - **Arm.** Coloc ran on the ea_only arm. Where the ea_only lead differs (`lead_differs_ea`
     in `pair_feasibility.parquet`), the within-donor test is about a different variant than
     the one coloc paired with the GWAS, so orientation needs LD between the two leads.
   - **Significance.** Only pairs with `qval < 0.05` carry a direction. Coloc-nominated pairs
     outside the gate family have no q, and their p-values are nominal.
