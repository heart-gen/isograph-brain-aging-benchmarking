# Nova2-cKO perturbation validation specification

## 1. Biological question

Do frozen IsoGraph `NOVA_FAMILY` switch-specific intronic windows preferentially
localize significant splicing responses to Nova2 loss in genetically defined mouse
neuronal populations?

## 2. Manuscript claim supported

This is orthogonal functional support for a human co-switch candidate set. A positive
result would show that exact-RBP perturbation-responsive splice neighborhoods are
preferentially localized to the switch-specific member of matched transcript pairs.
It would not establish adult-human occupancy, prove that NOVA2 is the sole causal
family member, or convert the discovery label from `NOVA_FAMILY` to NOVA2.

## 3. Hypothesis and expected direction

Within the frozen matched switch-window pairs, the candidate-bearing case window is
more likely than its sequence-matched control window to overlap a significant
Nova2-cKO splicing-event flank. The expected matched odds ratio is greater than one.

## 4. Statistical unit

Genes are the primary inferential unit. Multiple candidate transcript pairs and
multiple event flanks within a gene are collapsed with an any-support rule before
testing. Matched windows are a secondary diagnostic unit; modules are a secondary
stratified unit.

## 5. Comparison groups

For every frozen candidate transcript pair, the case window is the intronic flank
specific to the motif-bearing transcript and the control is its pre-matched
motif-negative flank. Each pair shares `match_id`, width, splice-flank class, and the
window-stage sequence matching constraints. Emx1 and Gad2 are tested separately as
primary cortical contexts. Pcp2 is a secondary cerebellar context. Contexts and
evidence tiers are never pooled.

## 6. Covariates and confound control

The primary matched design controls the window-stage matching factors: window length,
GC content, ambiguous sequence, splice-flank class, and source-gene context. Reciprocal
mm10-to-hg38-to-mm10 mapping is required to reduce cross-species coordinate artifacts,
with chromosome, strand, and at least 95% source-interval recovery. No regression-based
gene-set enrichment is attempted because the published GSE103314 supplements contain
significant events rather than the complete tested-event universe.

## 7. Primary analysis

The study-processed Quantas event table and BED12 coordinates are joined by event ID.
Only events with FDR at most 0.05 and absolute delta inclusion of at least 0.10 are
eligible. Intronic splice flanks of 100 nt are derived from adjacent BED12 blocks,
reciprocally mapped to hg38, and intersected with frozen human switch windows using
chromosome, strand, and splice-flank class. Event support is collapsed to case/control
indicators per gene. The exact two-sided McNemar sign test is applied to discordant
genes, with a Haldane-corrected matched odds ratio and 95% confidence interval.

## 8. Secondary and sensitivity analyses

- Repeat localization at 50 nt and 250 nt.
- Report matched-window effects as a diagnostic of the gene-level collapse.
- Estimate module-specific effects and apply Benjamini-Hochberg correction within the
  primary-width module family.
- Report cortical context replication without pooling Emx1 and Gad2.
- Intersect perturbation-localized windows with the independent stage-26 NOVA2 cTag-CLIP
  calls to identify exact direct-plus-responsive evidence. This remains descriptive
  because the processed occupancy maps are pooled and lack matched input.

## 9. Decision thresholds

An overall context is formally estimable only with at least 10 discordant genes. A
positive primary result requires matched odds ratio greater than one and exact
two-sided p-value at most 0.05 at 100 nt. Module findings require q-value at most 0.05
and the same direction. Sparse results are labeled descriptive or not estimable rather
than treated as negative binding evidence.

## 10. Missingness and exclusions

Events are excluded if identifiers do not reconcile between the Quantas and BED12
files, coordinates are malformed, adjacent BED blocks do not define a positive intron,
or reciprocal mapping fails. A window without overlap to the significant-only event
catalog is treated as not localized, not as proof of biological non-response. Ambiguous
multi-gene labels do not determine genomic overlap and are retained only as source
annotations; IsoGraph membership supplies the inferential gene identity.

## 11. Outputs and provenance

The analysis writes parsed events, derived mouse flanks, reciprocal mappings, per-window
localization calls, matched-pair, candidate, gene, module, context, replication, and
direct-plus-responsive tables; a JSON provenance record; and a Markdown report under
`07_rbp_regulation/_m/neuronal_clip/nova2_perturbation/`. The configuration pins source URLs,
byte sizes, SHA-256 checksums, thresholds, assemblies, frozen candidate identity, source
hashes, and output hashes. Heavy execution uses the stage-27 SLURM wrapper with seed 13.

## 12. Interpretation guardrails and reviewer risks

The strongest risk is ascertainment: the public tables are significant-event catalogs,
not full tested universes. Therefore the analysis tests spatial localization within
matched IsoGraph windows and makes no genome-wide enrichment claim. Cross-species loss
of mappability limits power and is reported explicitly. Event-flank overlap supports a
splice neighborhood, not necessarily the exact causal NOVA2 nucleotide. Directional
concordance is not tested because the IsoGraph/NOVA-family nomination has no independently
defined signed regulatory direction. A null result is an informative power- and
localization-limited outcome and does not negate the discovery nomination.
