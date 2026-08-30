# Neuronal CLIP validation of IsoGraph co-switch regulons

**Status:** analysis specification; implementation and data acquisition not started by
this document  
**Target:** Cell Genomics real-data validation  
**Date frozen:** 2026-07-30  
**Primary seed:** 13

## Analysis intent

We will test whether switch-relevant RNA sequence predicted to bind a prespecified
neuronal splicing factor is more frequently occupied by that exact factor in human
neuronal or human-brain eCLIP than matched non-switch sequence, supporting the claim
that a subset of IsoGraph co-switch modules has direct neuronal RBP-binding support.

This is a validation of a candidate **trans-regulatory mechanism**. It does not replace
the sQTL/eQTL specificity result, establish that an RBP causes module-wide
coordination, or imply that IsoGraph is globally superior to WGCNA.

## 1. Biological Claim

The intended claim is:

> Some independently nominated IsoGraph co-switch regulons are physically occupied
> by their exact neuronal RNA-binding protein at the sequence that distinguishes the
> switched transcript pair.

The claim must remain RBP- and context-specific. A peak anywhere in a member gene is
not validation. A cancer-cell-line peak is evidence of binding capacity, not
brain-specific occupancy. CLIP alone establishes physical binding support, not the
direction or causal effect of regulation.

The following stronger claims are out of scope unless separately demonstrated:

- one RBP coordinates an entire module;
- binding causes the observed age or disease association;
- a bulk-cortex peak is specific to a neuronal cell type;
- neuronal CLIP support is enriched specifically in GO-invisible modules;
- IsoGraph concentrates neuronal binding better than all WGCNA baselines.

## 2. Primary Hypothesis

For the prespecified exact-RBP predictions, the odds of a reproducible human neural
eCLIP peak overlapping the switch-relevant regulatory window are greater than the
odds for matched control windows from the same gene.

The primary estimand is the common matched odds ratio for:

`switch-relevant window versus matched within-gene control window`

across Tier-1 human neuronal or human-brain eCLIP datasets. Each peak map is used only
for its assayed RBP. The primary test is two-sided even though the expected direction
is odds ratio greater than 1.

### Frozen RBP family

The primary family is:

- TARDBP/TDP-43
- NOVA1
- NOVA2
- RBFOX2
- PTBP2

Aliases must be resolved in configuration rather than by fuzzy name matching. G3BP1,
FUS, KHDRBS1/SAM68, SNRNP70, and other RBPs are exploratory unless they are added in a
new dated specification before their neuronal CLIP results are inspected.

### Independent nomination rule

An RBP prediction is eligible only when it was nominated independently by the
existing motif analysis:

1. source artifact:
   `07_rbp_regulation/_m/rbp/rbp_regulon_intronic.parquet`;
2. exact RBP is one of the five frozen factors above;
3. module-by-RBP motif enrichment has BH `q <= 0.05`;
4. the gene has an RBP site-switch call in
   `07_rbp_regulation/_m/rbp/rbp_switch_calls_intronic.parquet`;
5. the gene's transcript pair is present in the corresponding canonical
   `structure_switch_pairs.parquet`.

The current frozen artifact contains 17 significant intronic module-by-RBP
nominations across these five factors. The implementation must write the exact
candidate rows, input file hashes, and creation date to a candidate manifest before
downloading or reading neuronal CLIP peaks.

The intronic scope is primary because all five factors are canonical splicing
regulators and because NOVA1/NOVA2 nominations are lost when mature and intronic
motif presence is unioned. Mature-transcript and combined-scope nominations are
secondary families and must not be merged into the primary multiple-testing family.

## 3. Statistical Unit and Design

### Unit of analysis

The genomic unit is a unique:

`gene × transcript pair × RBP × CLIP biological context`

Repeated appearances of the same gene or transcript pair across brain regions are
not independent. Candidate intervals must therefore be collapsed within each
RBP-context analysis, and uncertainty must be clustered or block-bootstrapped by
gene.

CLIP reads, peaks, and genomic windows are not biological replicates. Biological
replicates or donors are used to define a reproducible peak map; they are not counted
as independent observations in the interval-enrichment test. The primary inference
generalizes over nominated gene-switch targets in the assayed contexts, not over a
population of human donors.

### Primary human datasets

| RBP | Context | Accession | Access and primary use |
| --- | --- | --- | --- |
| PTBP2 | Human BA4 cortex; human iPSC-derived neurons | [GSE206650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206650), within GSE206661 | Public. The two context-labeled GEO consensus bigBeds are byte-identical and are therefore a shared sensitivity artifact only. Reconstruct BA4 and iPSC-neuron binding separately from each context's three IP and three size-matched-input replicates before either is primary-eligible. |
| TDP-43 | Human NGN2 neurons and directly converted neurons | [GSE276985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276985) | Public. Unstressed NGN2 is primary-eligible. Directly converted-neuron replicates contain one total enriched window and zero at q <= 0.05, so that context is `not estimable` rather than negative. Stressed NGN2 is a sensitivity, not pooled with baseline. |
| TDP-43 | Human iPSC-derived motor neurons, two control lines | [EGAS00001005880](https://ega-archive.org/studies/EGAS00001005880), EGAD00001008426 | Optional controlled-access extension; analyze lines separately and as a study consensus if obtained. It is excluded from the primary pooled estimate and missing access does not block completion. |
| NOVA1, NOVA2, RBFOX2 | Human iPSC-derived motor neurons, two control lines | [EGAS00001005880](https://ega-archive.org/studies/EGAS00001005880), EGAD00001008428 | Optional controlled-access extension. It can broaden exact-RBP support but is excluded from the primary pooled estimate and missing access does not block completion. |

Human contexts must never be unioned into one peak track. Study-provided
replicate-reproducible peaks are primary. If a study provides only replicate-level
peaks, a consensus requires support in at least two biological replicates or the
study's documented IDR/consensus rule.

### Orthogonal evidence tiers

The evidence tiers are analyzed and reported separately:

1. **Tier 1:** exact RBP, human neuronal or human-brain eCLIP;
2. **Tier 2:** exact RBP, cell-type-specific mouse cTag-CLIP;
3. **Tier 3:** exact RBP, bulk human or mouse brain CLIP/HITS-CLIP;
4. **Tier 4:** exact RBP, ENCODE K562/HepG2 eCLIP;
5. **Tier 5:** motif prediction only.

Tier 2 uses NOVA2 cTag-CLIP from
[GSE103315/GSE103316](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103316):
Emx1-positive cortical excitatory neurons, Gad2-positive cortical inhibitory neurons,
and Pcp2-positive cerebellar Purkinje cells. Emx1 and Gad2 may be contrasted within
cortex; Pcp2 must be interpreted separately because both tissue and cell type differ.

Tier 3 may include PTBP2 embryonic mouse-neocortex HITS-CLIP
([GSE47564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47564)) and NOVA1/2
embryonic mouse-cortex HITS-CLIP
([GSE69710/GSE69711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69711)).
FUS datasets GSE37190 and GSE40651 remain out of the frozen primary RBP family.

Tier 4 is the existing ENCODE analysis. Its K562/HepG2 results must not influence
candidate selection, substitute for missing Tier-1 data, or be pooled into the
primary estimate.

### Cross-species rule

Mouse evidence requires reciprocal mapping of the human switch window to the mouse
assembly and back to hg38, with preserved strand and interval class. The tested
universe is restricted to reciprocally mapped intervals; unmappable human sequence is
`not testable`, not unbound. Primary and control windows must be subjected to the same
orthology and mappability filters.

## 4. Primary Comparison

### Switch-relevant case windows

All coordinates use GENCODE v47 on GRCh38 and strand-aware BED conversion. For each
transcript pair, construct non-overlapping case windows by structural class:

- internal switched exon: the exon sequence present in only one transcript plus up
  to 100 nt of each adjacent intron;
- alternative donor or acceptor: the alternative splice-site segment plus 100 nt of
  intronic flank;
- transcript-specific first exon: the transcript-specific exon plus the adjacent
  donor-side intronic flank;
- transcript-specific last exon: the transcript-specific exon plus the adjacent
  acceptor-side intronic flank;
- other exonic sequence present in one transcript but absent from the other:
  exonic symmetric difference, classified as CDS or UTR where possible.

The primary intronic flank is **100 nt**, matching the completed intronic motif
analysis. A 250-nt flank is a prespecified sensitivity only. Window size must not be
changed after seeing overlap results.

Overlapping case components are merged within a transcript pair before calculating
callable nucleotide opportunity. A peak is switch-localized only when it intersects a
case window on the correct strand where strand is available.

### Matched controls

The primary control is non-switch sequence from the same gene. Each case is matched
to up to 10 deterministic control windows with:

- the same exonic/intronic and CDS/UTR class;
- equal total callable length within 10%;
- GC content within 5 percentage points;
- similar splice-site distance;
- similar uniquely mappable fraction;
- detectable input-library coverage in the same CLIP context.

Controls must exclude every case window from every eligible transcript pair for that
gene. First/last-exon cases must be matched to the corresponding terminal-exon class;
they must not be matched to an arbitrary internal constitutive exon.

If no valid within-gene control exists, the pair is excluded from the primary test.
It may enter a secondary background-gene analysis using expression-, isoform-count-,
transcript-number-, intron-length-, GC-, mappability-, and sequence-opportunity-matched
switch genes from the same IsoGraph analysis.

### Specificity contrast

As a necessary secondary control, fit the interaction:

`switch-relevant interval × independently nominated RBP-gene prediction`

using the same exact-RBP peak map. This asks whether switch localization is stronger
for nominated predictions than for otherwise matched, non-nominated switch events.
It protects against concluding validation merely because an RBP binds introns or
alternative exons broadly.

## 5. Covariates and Exclusions

Most confounding is handled by matching. Any residual model may include only:

- callable interval length;
- GC fraction;
- uniquely mappable fraction;
- CLIP input coverage or matched expression;
- sequence class;
- distance to the nearest splice site;
- dataset/RBP context.

Age, diagnosis, sex, and tissue-composition covariates from the IsoGraph discovery
cohorts do not enter this external interval-overlap model. The discovery modules and
switch pairs are fixed before CLIP validation.

Exclude or mark not testable:

- transcript IDs absent from GENCODE v47;
- transcript pairs mapping to different genes or chromosomes;
- ambiguous or non-strand-resolved intervals when strand is required;
- genes without detectable input coverage in the assayed CLIP context;
- windows without a valid primary match;
- blacklisted, low-mappability, or assembly-gap sequence;
- mouse intervals failing reciprocal liftOver;
- datasets lacking an exact RBP match or a defensible input/peak definition.

Peak absence in an unexpressed gene must never be coded as evidence against binding.

## 6. Expected Effect Direction

The primary matched odds ratio is expected to exceed 1. The absolute paired
difference in bound fraction must be reported with the odds ratio so that statistical
significance cannot obscure a negligible effect.

No signed splicing direction is inferred from a CLIP peak alone. A perturbation is
called directionally concordant only when all of the following were defined without
using the perturbation result:

1. a signed exon/transcript event maps unambiguously to the IsoGraph pair;
2. the RBP's expected positional effect on inclusion or exclusion is prespecified;
3. the module/trait direction and direction of RBP activity are independently
   identifiable.

Otherwise, the correct label is **direct and perturbation-responsive**, not
directionally concordant.

## 7. Minimal Primary Analysis

### Primary model

For each exact RBP and human neural context:

1. collapse repeated region-level appearances to unique gene-transcript pairs;
2. create case/control matched sets;
3. score overlap with the context-specific reproducible eCLIP peak map;
4. fit conditional logistic regression with case status as the exposure and peak
   overlap as the outcome, stratified by matched set;
5. obtain uncertainty by a gene-block bootstrap with 2,000 draws and seed 13;
6. report odds ratio, 95% confidence interval, two-sided p-value, absolute paired
   bound-rate difference, number of genes, matched sets, and informative discordant
   sets.

The single primary endpoint is a one-stage stratified estimate across the public
Tier-1 dataset-by-RBP contexts: separately reconstructed PTBP2 BA4 and iPSC-neuron
maps and baseline NGN2 TDP-43. Strata remain at the dataset/RBP/matched-set level and
uncertainty is gene-blocked. Dataset-specific and RBP-specific estimates must always
be shown beside the pooled value. EGA contexts, if obtained, form a separate optional
extension and are not added to the prespecified primary pool.

A context with fewer than 10 informative discordant matched sets receives an effect
estimate and exact interval but no confirmatory p-value. Optional controlled-access
datasets are reported separately when available and otherwise omitted without
changing the primary denominator.

### Multiple testing

- One pooled primary test: two-sided alpha 0.05.
- Public-primary RBP-specific tests: BH correction across TARDBP and PTBP2.
- Optional EGA RBP-specific tests, if run: a separate BH family across the exact
  RBP tests available in the two controlled datasets.
- Module-by-RBP follow-up tests: BH correction across the frozen eligible
  module-by-RBP family.
- Evidence tiers and mature/combined motif scopes are separate secondary families.

### Success and failure gates

The broad neuronal-binding claim is supported only when:

- the pooled Tier-1 odds ratio is greater than 1 with a 95% confidence interval
  excluding 1;
- both public-primary RBPs have concordant estimates greater than 1;
- leave-one-RBP-out pooled estimates remain greater than 1; and
- the effect survives input-coverage and matched-background sensitivities.

An individual module is **binding-supported** only when it has at least three
supported member genes, a module-by-RBP BH `q <= 0.05`, and enrichment greater than
1. A single illustrative locus cannot upgrade a module.

If the pooled confidence interval includes 1, enrichment appears only for
non-switch gene sequence, or the result disappears after expression/opportunity
matching, the primary hypothesis is not supported. A Tier-4-only signal is also a
failure of neuronal validation. Such a result leaves the existing genetic anchoring
intact and the RBP layer motif-supported only.

## 8. Sensitivity Analyses

The following are necessary and prespecified:

1. **Window width:** primary 100-nt intronic flank versus 250 nt and 50 nt.
2. **Localization:** any case-window overlap versus exact motif-instance overlap and
   motif instance padded by 50 nt.
3. **Peak representation:** context-valid study-provided consensus peaks versus a
   uniform input-aware reprocessing where raw data and compute permit. For GSE206650,
   context-specific IP/input reconstruction is required and the duplicated published
   consensus is sensitivity-only.
4. **Replicate rule:** consensus peaks versus peaks supported separately in each
   biological replicate.
5. **Control source:** strict within-gene controls versus matched background switch
   genes.
6. **Detectability:** alternative input-coverage thresholds; genes near the threshold
   excluded.
7. **Outcome:** binary peak overlap versus peak-covered nucleotide fraction and
   input-normalized signal.
8. **Duplicate genes:** one representative transcript pair per gene, selected by
   largest pre-CLIP switch strength, versus all unique pairs with gene blocking.
9. **Structural class:** internal exon/splice-site events only; terminal-exon events
   analyzed separately.
10. **Module class:** GO-invisible versus GO-visible interaction. This is secondary;
    a nonsignificant interaction must not be described as GO-invisible specificity.
11. **Phenotype restriction:** phenotype-associated modules and the SCZD
    GO-invisible modules analyzed as a separate targeted family only after the same
    motif-nomination rule is applied without inspecting neuronal peaks.
12. **Cross-species mapping:** reciprocal-liftOver threshold sweep and conserved
    splice-junction-only analysis for mouse evidence.
13. **Influence:** leave-one-dataset, leave-one-RBP, and leave-one-gene-out analyses.
14. **Peak-density negative control:** rotate or permute peaks within expressed genes
    while preserving chromosome, gene, strand, and peak length; 2,000 draws, seed 13.
15. **Matched-method control:** compare concentration in
    `wgcna_switch_only` and `wgcna_multiplex` modules on identical switch features as
    exploratory method specificity. This analysis cannot support global superiority.

No sensitivity may replace the primary result. A signal that appears only after
widening windows or relaxing matching is exploratory.

## 9. Orthogonal Validation

### Cell-type-specific mouse CLIP

NOVA2 cTag-CLIP is used to ask whether binding support differs between Emx1-positive
and Gad2-positive cortical neurons. The test is a cell-type-by-case-window
interaction using reciprocally mapped intervals. Purkinje-cell results are reported
as a separate cerebellar context, not as a third level of the cortical comparison.

Mouse cTag-CLIP can corroborate exact-RBP occupancy and suggest cell-type preference.
It cannot establish that the same occupancy occurs in adult human cortex.

### Perturbation

Use perturbation only when a switch event maps unambiguously:

- PTBP2 depletion/splicing data within GSE206661;
- NOVA2 conditional-knockout RNA-seq in
  [GSE103314/GSE103316](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103314);
- other exact-RBP neuronal perturbations only after accession, contrast, and
  thresholds are added to the frozen configuration.

A direct perturbation-supported event requires:

- exact-RBP CLIP overlap in the switch-relevant window;
- event-level FDR `<= 0.05`;
- absolute delta-PSI `>= 0.10`, or a prespecified transcript-usage effect of
  comparable magnitude;
- consistent event mapping between the perturbation annotation and GENCODE v47.

Test enrichment of perturbation-responsive events among CLIP-supported predictions
against the same matched switch background. Do not count individual exons as
independent when they arise from the same gene; block or collapse by gene.

### Neuronal aging and stress

GSE276985 age-stratified mouse-brain TDP-43 eCLIP and stressed human NGN2 neurons are
secondary context-modification analyses. They may test whether occupancy changes
with age/stress at already supported loci. They must not be pooled with baseline
occupancy to manufacture primary support.

## 10. Main Figure-Worthy Result

Only a successful analysis should enter a main or supplementary figure. The anchor
panel is:

- a forest plot of matched odds ratios for each Tier-1 RBP/context and the pooled
  estimate, colored by human tissue or neuronal model;
- an adjacent heat map of binding-supported module-by-RBP pairs, annotated for
  GO-invisibility and phenotype association;
- one locus panel showing the two transcript structures, switch-relevant window,
  exact-RBP neural eCLIP peak, and perturbation delta-PSI when available.

The forest plot is the inferential result; the locus track is an illustration. If
only one RBP succeeds, report an RBP-specific supplementary result rather than a
general neuronal-regulon panel. Promotion must also respect the Cell Genomics
display-item limit already documented in `MANUSCRIPT_PLAN.md`.

## 11. Reviewer Objections and Responses

| Likely objection | Prespecified response |
| --- | --- |
| Peaks reflect gene expression, not switch-specific binding. | Use input-detectable genes, strict within-gene matched windows, and an input-coverage sensitivity. |
| Long introns or exons have more opportunity for a peak. | Match callable length, sequence class, splice distance, GC, and mappability; report covered fraction as a sensitivity. |
| The candidate list was selected after seeing CLIP. | Freeze the five RBPs and exact motif-derived candidate manifest, with hashes, before neuronal peak ingestion. |
| Repeated genes across regions inflate sample size. | Collapse identical gene-transcript pairs and block all uncertainty by gene. |
| Reads or peaks are being treated as biological replicates. | Replicates define consensus peak maps; inference is over unique switch targets, with the limited donor generalizability stated explicitly. |
| Cancer-cell-line eCLIP is not neural evidence. | Keep ENCODE as Tier 4 and exclude it from the primary estimate. |
| Mouse binding is not human binding. | Keep mouse results orthogonal, require reciprocal mapping, and label unmappable sequence not testable. |
| Different studies use different peak callers. | Use study consensus peaks for the primary within-study analysis and uniform reprocessing only as a sensitivity. |
| A binding peak does not prove regulation. | Require matched perturbation response for functional support and reserve directional language for prespecified signed events. |
| Any peak anywhere in a target gene is being counted. | Require intersection with a structurally defined switch-relevant window and compare it with matched sequence from the same gene. |
| One or two loci drive the result. | Require gene-blocked inference, leave-one-gene-out analyses, at least three genes per supported module, and breadth across at least two RBPs for the global claim. |
| The result is generic alternative-exon binding. | Test the nominated-by-switch-window interaction against non-nominated matched switch events. |
| GO-invisible specificity or method superiority is post hoc. | Treat GO class and matched WGCNA comparisons as secondary interactions and retain the complementary, not globally superior, framing. |

## 12. Implementation Notes

### Reproducible interfaces

Implementation requires committed, parametrized modules under
`isograph_benchmark/real_data/`, not inline result-generating scripts:

- `neuronal_clip_validation.py`: manifest validation, switch-window construction,
  matching, peak overlap, primary and module-level tests;
- `neuronal_clip_meta.py`: dataset/RBP synthesis, sensitivity summaries, report
  generation;
- `neuronal_clip_perturbation.py`: event mapping and perturbation support.

All dataset IDs, URLs, assemblies, RBP aliases, peak rules, access states, matching
tolerances, window widths, thresholds, and seed belong in a new
`configs/neuronal_clip.yaml`, which is the single source of truth.

Public-data acquisition is submitted with
`inputs/_h/download_neuronal_clip.sh`; the wrapper calls the
configuration-driven downloader with safe resume, checksum verification, and archive
extraction.

The implemented pre-CLIP window stage is run with:

```bash
sbatch 07_rbp_regulation/_h/09.neuronal_clip_windows.sh
```

The wrapper calls `neuronal_clip_validation.py prepare-windows` and writes under
`07_rbp_regulation/_m/neuronal_clip/windows/`. It preserves every frozen nomination in
`candidate_window_membership.parquet`, collapses them to canonical RBP/gene/transcript
pairs, constructs pair-specific intronic flanks at 50, 100, and 250 nt, and matches
case windows without replacement to structurally and sequence-matched windows from
the motif-negative transcript. Mappability and context-specific input coverage remain
explicitly deferred to the callable-window stage; no IP peaks enter pre-CLIP matching.

The post-window motif-opportunity gate is run with:

```bash
sbatch 07_rbp_regulation/_h/10.neuronal_clip_motif_qc.sh
```

This stage writes `candidate_eligibility.parquet` plus strand-aware exploratory
NOVA-family YCAY instances and clusters under
`07_rbp_regulation/_m/neuronal_clip/motif_opportunity/`. It reproduces the finding that the
current NOVA1 presence/absence nomination is almost identical to whether only one
isoform has intronic sequence. Current NOVA1 rows are therefore retained for
provenance but labeled `not_testable_current_candidate_definition` and excluded from
confirmatory overlap. The YCAY coordinates do not rescue those nominations and must
not be reported as NOVA1-specific evidence.

An independent, opportunity-controlled NOVA-family re-nomination is run with:

```bash
sbatch 07_rbp_regulation/_h/12.nova_family_renomination.sh
```

This stage rescans the complete 17-analysis switch-pair universe and requires both
transcripts in a pair to provide at least 100 nt of intronic opportunity. The primary
call is an exact YCAY cluster (at least three `CCAC`, `CCAT`, `TCAC`, or `TCAT`
instances within 30 nt). Module enrichment uses genes as the statistical units and
adjusts for total intronic opportunity, eligible-pair count, and transcript count;
candidate output remains pair-level. All analysis keys include cohort tree and region
so the BrainSEQ and GTEx hippocampus analyses cannot collide.

The frozen run scored 79,763 transcripts and 374,296 eligible region-specific pairs.
Two adjusted regulons passed global BH q <= 0.05: BrainSEQ hippocampus M001 and M002,
both GO-visible. They contain 339 switched genes represented by 1,993 unique transcript
pairs. These rows are labeled `NOVA_FAMILY`, never NOVA1 or NOVA2. They are a testable
rescue for orthogonal NOVA-family CLIP validation, not evidence for the GO-invisible
disease-module claim. Exact-RBP attribution requires public NOVA2 cTag-/HITS-CLIP or
the optional controlled-access neuronal eCLIP extension.

Public Tier-2 NOVA2 cTag-CLIP validation is run with:

```bash
sbatch 07_rbp_regulation/_h/13.nova2_ctag_clip.sh
```

The configuration pins the GSE103315 processed archive and reciprocal UCSC
`hg38ToMm10`/`mm10ToHg38` chain files by byte count and SHA-256. The stage rebuilds
50-, 100-, and 250-nt strict within-gene case/control windows from the frozen 1,993
NOVA-family pairs before reading CLIP signal. A human window is mouse-testable only
when it maps completely within one alignment block, has a unique best forward and
reverse chain mapping, recovers chromosome and strand, and overlaps its original
human interval by at least 95% after reciprocal mapping.

GSE103315 provides one study-processed unique-tag coverage bedGraph for each selected
cell type: Emx1-positive cortical excitatory neurons, Gad2-positive cortical inhibitory
neurons, and Pcp2-positive Purkinje cells. Each map pools three biological replicates;
replicate-resolved processed peaks and matched input are not supplied. Consequently:

- exact NOVA2 occupancy is scored as overlap with positive pooled tag coverage, not
  as overlap with a called replicate-consensus peak;
- gene-level case/control effects, module effects, and the Emx1-by-Gad2 interaction
  are labeled descriptive pooled-signal tests regardless of p-value;
- Pcp2 is reported separately as a cerebellar context;
- zero pooled coverage is never interpreted as proof of no binding; and
- upgrading to replicate-consensus evidence requires a separately configured raw-SRA
  reprocessing stage.

The CLI writes each context atomically under
`07_rbp_regulation/_m/neuronal_clip/nova2_ctag/context_calls/` as soon as its scan finishes,
then writes reciprocal windows, window-, candidate-, and gene-level support, overall
and module effects, the cortical interaction, source/input/output hashes, and
`NOVA2_CTAG_VALIDATION.md`. Discovery rows remain labeled `NOVA_FAMILY`; `NOVA2` labels
only this exact-RBP orthogonal evidence.

**Frozen stage-26 result (2026-08-01; SLURM 42902680).** Of 42,718 matched
case/control window rows across widths, 1,272 (3.0%) passed the strict reciprocal
mapping rule. At the prespecified 100-nt width, only 9 genes retained jointly callable
case/control windows; neither side had positive pooled NOVA2 coverage in Emx1, Gad2,
or Pcp2. The corrected Emx1-by-Gad2 interaction therefore had 9 jointly callable genes
and 0 informative differences. The 50-nt Emx1 sensitivity retained 76 genes and was
sparse and nonsignificant (4 case-bound versus 2 control-bound; Haldane OR 1.8, exact
p=0.6875); at the matched-window level it was exactly balanced (2 case-only versus 2
control-only; OR 1.0, p=1.0). This is an informative null caused primarily by
cross-species interval attrition and does not upgrade the frozen discovery label from
`NOVA_FAMILY` to NOVA2.

Raw replicate reprocessing is not the next primary implementation. Although the nine
selected SRA-lite runs total only about 474 MB, replicate-level peak recovery cannot
repair the dominant reciprocal-mapping bottleneck and the current designated
environment lacks the required pinned SRA/alignment toolchain. The next orthogonal
stage is event-level mapping of the public GSE103314 Nova2-cKO splicing tables, kept as
a perturbation tier separate from cTag occupancy. Raw reprocessing remains an optional
sensitivity if perturbation-supported loci justify it.

Callable-window overlap is run only after that gate succeeds:

```bash
sbatch 07_rbp_regulation/_h/11.neuronal_clip_overlap.sh
```

The overlap CLI reads only candidates marked `confirmatory_eligible`. For TDP-43,
callability requires that the matched pair's gene occur in each replicate's
study-processed CLIP universe; binding requires strand-correct overlap with a q <=
0.05 enriched window in both configured replicates. For PTBP2, the two
context-labeled consensus bigBeds remain sensitivity-only; the primary context maps
are reconstructed as replicate-paired, library-normalized IP/input signal from the
strand-specific bigWigs. PTBP2 results are explicitly labeled signal support rather
than reconstructed CLAM peaks. Contexts failing the frozen context-QC manifest remain
`not estimable`, and controlled EGA contexts remain optional.

Heavy preprocessing, raw-read reanalysis, reciprocal liftover at scale, and full
overlap testing must run on SLURM under account `bio260021p` with 2000M per CPU.
Login-node work is limited to manifests, small aggregation, QC, and plotting.

Use:

- Python:
  `/ocean/projects/bio260021p/shared/opt/envs/isograph/bin/python`;
- figure R:
  `/ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript`;
- deterministic seed 13 for matching, bootstrap, and permutation;
- `Path.iterdir()` rather than `glob()` on Ocean FS.

### Required manifests and outputs

Before analysis, write:

- `dataset_manifest.tsv`: accession, dataset and sample IDs, RBP, species, tissue or
  cell type, condition, assay, evidence tier, assembly, replicate and input IDs,
  peak-calling source, access state, retrieval date, URL, local path, and checksum;
- `context_qc.tsv`: replicate and file completeness, enriched-window yield,
  confirmatory eligibility, controlled-access state, and the prespecified reason for
  every context-level exclusion;
- `candidate_manifest.parquet`: region, module, GO class, phenotype class, RBP,
  motif scope, nomination q-value, gene, transcript pair, and source-file hashes;
- `window_manifest.parquet`: case/control coordinates, structural class, strand,
  match ID, length, GC, mappability, input coverage, orthology status, and exclusion
  reason.

Analysis outputs:

- `overlap_calls.parquet`;
- `primary_effects.parquet`;
- `module_effects.parquet`;
- `perturbation_support.parquet`;
- `sensitivity_effects.parquet`;
- `neuronal_clip_validation.json`;
- `NEURONAL_CLIP_VALIDATION.md`;
- figure-ready CSV plus PDF and PNG.

Raw and processed CLIP caches belong under `inputs/raw/neuronal_clip/` and must not be
committed. Large `_m/` outputs and region output directories must not be committed.
The concise specification, configuration, CLI modules, SLURM wrapper, small manifests,
and final summary are the reproducibility record, subject to the repository's normal
selective-staging policy.

### Staged execution

1. Freeze and hash the candidate and dataset manifests.
2. Acquire public Tier-1 files; record EGA as an optional, non-blocking extension.
3. Build and QC GENCODE v47 switch windows and deterministic matches.
4. Run the human Tier-1 primary analysis.
5. Run required sensitivities and negative controls.
6. Add Tier-2/3 mouse evidence and perturbation support without pooling tiers.
7. Generate the summary, null-aware decision statement, and conditional figure.

The implemented stages are configuration-driven and SLURM-wrapped. Public NOVA2
acquisition and pooled Tier-2 overlap are stage 26. Replicate-level raw reprocessing,
additional Tier-3 overlap, and perturbation analysis remain separate stages and must
preserve the evidence-tier boundaries above.
