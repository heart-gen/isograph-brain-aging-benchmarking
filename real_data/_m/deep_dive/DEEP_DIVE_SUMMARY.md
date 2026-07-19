## Per-gene mechanistic deep-dive of colocalized disease genes

### Purpose
Move beyond the aggregate co-switch-gene review to a per-gene mechanistic integration: for every
gene whose GWAS locus colocalizes with a brain QTL, ask whether the disease variant acts through
an IsoGraph isoform switch (splicing-led, DTU-without-DGE) or through gene-level abundance
(expression-led), and assemble the supporting layers into a single vignette per gene.

### Inputs
- `real_data/coloc/_m/coloc_isoform_events_combined.parquet` (141 colocalized isoform events,
  68 genes; eCAVIAR CLPP, sQTL/eQTL, tissue, junction→transcript→switch-pair mapping,
  `structural_consequence`, GTEx concordance, GO-invisibility, BrainSeq replication).
- `real_data/coloc/_m/coloc_direction_combined.parquet` (signed risk-allele direction).
- `real_data/_m/rbp/{rbp_switch_calls,rbp_regulon}.parquet` (RBP motifs switched per gene and
  module-level enrichment).
- `real_data/*/*/_m/isograph_vae/clinical_consequence/{gene_constraint,exon_clinvar}.parquet`
  (gnomAD LOEUF/missense o/e; per-exon ClinVar P/LP and CDS-overlap flags).
- GENCODE v47 exon annotation for the SNCA transcript-pair schematic.

### Methods Text
For each colocalized gene we joined, keyed by Ensembl gene identifier, its colocalized isoform
events, signed risk-allele QTL direction, switched RBP motifs (retaining RBPs that were both
called switched in the gene and enriched in its IsoGraph module at BH q < 0.05), and gnomAD
loss-of-function constraint. Each gene was classified as splicing-led when a colocalizing splicing
QTL mapped onto a GTEx-concordant IsoGraph switch pair, expression-led when only gene-level eQTL
colocalizations were present, or splicing-unresolved when a colocalizing sQTL did not map onto the
tissue's switch pair. The analysis is a set of deterministic joins over existing result tables
(`gene_deep_dive.py`); no new statistical model or random component is introduced.

### Results Text
Across 68 colocalized genes, 12 were splicing-led (a colocalizing sQTL resolving onto a concordant
IsoGraph switch pair), 23 splicing-unresolved, and 33 expression-led. **All 12 splicing-led genes
reside in GO-invisible switch modules** — every genetically anchored isoform switch sits in a
module that gene-level pathway enrichment would miss. Four splicing-led genes carried multiple
resolved colocalizations: PPP6R2 (ALS and SCZ, three events), SNCA (LBD and PD), and CDIP1 and
DLG1 (two events each within SCZ). eCAVIAR CLPP values were modest (max 0.39 for CTSH in AD
hippocampus; most < 0.1), so individual colocalizations are suggestive rather than definitive; the
strength of the splicing-led set is its coherence — cross-disease concordance and GO-invisibility —
rather than any single high-posterior locus.

**SNCA cross-disease vignette (candidate main figure).** In both dementia with Lewy bodies (risk
allele A at rs7680557, cortex, CLPP 0.038) and Parkinson's disease (risk allele C at rs1471483,
frontal cortex BA9, CLPP 0.025), the risk allele increases usage of the same junction
chr4:89,835,692–89,836,127, which maps onto one IsoGraph switch pair (ENST00000508895 /
ENST00000618500) in a GO-invisible module. The junction defines an **alternative first exon**: the
canonical isoform (ENST00000336904) uses a distal first exon (chr4:89,838,252–89,838,315), whereas
the switch isoforms use the alternative first exon (chr4:89,836,127–89,836,213), consistent with
the "no annotated structural change" coding-consequence call.

**Clinical consequence of the SNCA switch.** The alternative first exon and the other 5′ switched
exons are non-coding (`cds_overlap = False`) and carry no ClinVar pathogenic/likely-pathogenic
(P/LP) variants (0 P/LP; the alternative first exon has 0 ClinVar records) in all five regions
where SNCA switches. In contrast, SNCA's P/LP burden falls in constitutive, isoform-shared coding
exons — chr4:89,828,143–89,828,184 (4 P/LP) and chr4:89,835,547–89,835,692 (1 P/LP) — outside the
switched region. The clinical consequence of the switch is therefore **regulatory (5′UTR
remodeling), not a protein-coding change**: the common-variant risk mechanism (isoform choice at a
dosage-sensitive gene; SNCA LOEUF = 0.40) is distinct from, and complementary to, the rare coding
pathogenic variants that cause Mendelian synucleinopathy. This case instantiates the manuscript
thesis twice over — the switch is invisible to GO/pathway enrichment (GO-invisible) *and* to
coding-variant/ClinVar analysis (non-coding) — yet is genetically anchored across two
synucleinopathies. It is consistent with the aggregate switch-consequence result (productive UTR
remodeling, not decay) and the aggregate clinical-consequence result (switched exons less
P/LP-dense than constitutive exons; ratio < 1).

### Figure and Table Notes
- Potential main figure: `real_data/_m/figures/figGeneticAnchoring.{pdf,png}`
  - Rationale: primary payoff of the genetic-anchoring layer; integrates variant → switch →
    heritability → resolved-gene set.
  - Key message: disease variants resolve to GO-invisible isoform switches; SNCA converges across
    LBD and PD on one alternative-first-exon switch.
  - Required legend details: panel A is a coordinate-accurate transcript schematic (minus strand,
    5′→3′ left→right); CLPP is eCAVIAR; panel B is single-annot S-LDSC enrichment of the aging
    switch layer; asterisks are enrichment p (`* <0.05, ** <0.01, *** <0.001`).
- Potential supplementary table: `real_data/_m/deep_dive/deep_dive_panel.parquet`
  (+ `DEEP_DIVE_PANEL.md`) — one row per colocalized gene with verdict, CLPP, LOEUF, concordant
  traits, multi-locus and GO-invisible flags.
- Per-gene vignettes: `real_data/_m/deep_dive/<GENE>.md` (68 files).

### Reproducibility Information
- Analysis directory: `real_data/_m/deep_dive/`
- Primary scripts: `isograph_benchmark/real_data/gene_deep_dive.py`;
  figure `real_data/_h/genetic_anchoring_figure.R`.
- Execution command: `python -m isograph_benchmark.real_data.gene_deep_dive`;
  `Rscript real_data/_h/genetic_anchoring_figure.R`.
- Execution date: 2026-07-19.
- Compute environment: PSC Bridges-2; conda env
  `/ocean/projects/bio260021p/shared/opt/envs/isograph` (Python; igraph 1.0.0) for the join,
  `/ocean/projects/bio260021p/shared/opt/envs/rnaseq` (R; arrow/ggplot2/patchwork) for the figure.
- Random seed: not applicable (deterministic joins).
- Git commit: uncommitted at time of writing (branch `clinical-consequence-light`); wrapper fixes
  committed at b9664bd.
- Key package versions: not recorded in a runtime log for this run; extract from the env before
  manuscript submission.
- Missing reproducibility information: package version pins for the deep-dive run.

### Limitations and Integration Notes
Colocalization posteriors are modest and bulk-tissue-derived, so per-gene claims are suggestive;
the set-level pattern (splicing-led ≡ GO-invisible; cross-disease concordance) is the defensible
claim. RBP motif calls are motif presence, not measured binding, and do not yet test whether the
lead QTL variant sits within a switched-exon motif.

A curated **literature layer** (deep-dive layer 6) is now attached to each resolved
splicing-led gene as vignette Section 6 and as `deep_dive_literature.{parquet,tsv}` (Table S12):
four genes have established disease isoform biology matching their resolved switch — SNCA
[@doi:10.3389/fgene.2019.00584; @doi:10.3390/genes9020063], DLG1/SAP97 [@doi:10.1038/tp.2015.154],
CTSH [@doi:10.1038/s41386-023-01542-2], and ARVCF [@doi:10.1038/sj.mp.4001586] — and the other
eight (PPP6R2, GGNBP2, PGS1, CDIP1, PRRC2B, RTEL1, TBC1D15, TPCN1) are flagged `novel_candidate`
(no established disease-specific isoform literature), the under-characterized GO-invisible
switching the method is designed to nominate.

The integrated genetic-anchoring Results section is now drafted at
`real_data/_m/GENETIC_ANCHORING_RESULTS.md`, folding this summary together with
`SWITCH_CONSEQUENCE_SUMMARY.md`, `RBP_REGULON_SUMMARY.md`, `CLINICAL_CONSEQUENCE_META.md`,
`LDSC_SUMMARY.md`, and the coloc `SIGNED_DIRECTION_ISOFORM_EVENTS_SUMMARY.md`.
