# Genetic anchoring of the IsoGraph switch layer — integrated Results section

Manubot-ready synthesis folding the genetic-anchoring and downstream-biology layers
(colocalization → partitioned heritability → coding consequence → RBP regulons → clinical
consequence → per-gene deep-dive + literature) into one Results narrative. It backs main
**Fig 3** (`figQtlSpecificity`) and **Fig 4** (`figGeneticAnchoring`) and supplements
**S-real-5/6/7**. Numbers are traceable to the per-analysis summaries cited at each subsection;
citekeys are Manubot-formatted, with `[citation needed: …]` where a standard tool/resource
citekey still needs pinning (never fabricated). North-star: IsoGraph is a **complementary
DTU-without-DGE layer**, and this section is the orthogonal-evidence payoff that the layer is
genetically real and disease-relevant.

---

## Disease variants colocalize onto GO-invisible isoform switches

Source: `05_genetic_anchoring/_m/coloc/SIGNED_DIRECTION_ISOFORM_EVENTS_SUMMARY.md`,
`05_genetic_anchoring/_m/deep_dive/DEEP_DIVE_SUMMARY.md`.

To ask whether the switch layer carries genuine disease genetics rather than technical
structure, we colocalized brain splicing- and expression-QTL credible sets (GTEx v11 sQTL/eQTL,
SuSiE) with GWAS credible sets for five traits — schizophrenia (SCZ), Alzheimer's disease (AD),
Parkinson's disease (PD), dementia with Lewy bodies (LBD), and amyotrophic lateral sclerosis
(ALS) — using eCAVIAR CLPP, and mapped each colocalized junction back onto the IsoGraph
switch-pair structure of the matching tissue [citation needed: eCAVIAR];
[citation needed: GTEx v11]; [citation needed: SuSiE]. Across 141 colocalized isoform events in
68 genes, 12 genes were **splicing-led** — a colocalizing sQTL resolving onto a
GTEx-concordant IsoGraph switch pair (the DTU-without-DGE class) — 23 were splicing-unresolved
(colocalizing sQTL not mapping onto the tissue's switch pair), and 33 were expression-led
(gene-level eQTL only). **All 12 splicing-led genes reside in GO-invisible switch modules**:
every genetically anchored switch sits in a module that gene-level pathway enrichment would
miss, exactly the complementary biology the method is built to surface. Nine of these were
GTEx-tissue-matched concordant aging events, again all GO-invisible, a conservative floor
(concordance requires the sQTL junction to fall within 2 bp of an annotated GENCODE junction).

The headline case is **SNCA** (α-synuclein), which converges across two synucleinopathies: in
LBD (risk allele A at rs7680557, cortex, CLPP 0.038) and in PD (risk allele C at rs1471483,
frontal cortex BA9, CLPP 0.025) the risk allele raises usage of the same junction
chr4:89,835,692–89,836,127, mapping onto one IsoGraph switch pair (ENST00000508895 /
ENST00000618500) in a GO-invisible module. The event is an **alternative first exon**: the
switch isoforms use chr4:89,836,127–89,836,213 in place of the canonical distal first exon
(ENST00000336904; chr4:89,838,252–89,838,315), a 5′-regulatory rather than protein-coding
change. eCAVIAR posteriors are modest (max 0.39, CTSH in AD hippocampus; most < 0.1), so the
strength of the splicing-led set is its **coherence** — cross-disease concordance and uniform
GO-invisibility — not any single high-posterior locus. Schizophrenia produced no GTEx-concordant
event but five BrainSeq-replicated splicing events (Fig 4C–D).

## The switch layer carries partitioned neurodegenerative heritability

Source: `05_genetic_anchoring/_m/ldsc/LDSC_SUMMARY.md` (Fig 4B).

Stratified LD-score regression on baselineLD v2.2 — robust to the gene-size confound that
inflates MAGMA on giant modules — confirmed that heritability concentrates in the switch layer's
cis-regulatory variants [citation needed: S-LDSC/LDSC]; [citation needed: baselineLD v2.2]. In
single-annotation models the aging switch layer's brain **sQTL** annotation was enriched for
heritability in every neurodegenerative trait: LBD 7.17× (enrichment p = 0.068), PD 3.80×
(p = 0.049), AD 3.37× (p = 0.0021), ALS 2.96× (p = 0.0022), and SCZ 1.74× (p = 0.0051); the
matched eQTL and total-cis annotations were comparably enriched. The one honest asymmetry is
**disease** SCZ (the BrainSeq-SCZD annotation), where the switch layer is expression-led — sQTL
enrichment 0.71× (n.s.) versus eQTL 2.04× (tau p = 0.016), and the joint sQTL-vs-eQTL model
assigns the signal to eQTL (tau p = 0.0064) — matching the coloc verdict that SCZ resolves
through abundance more than splicing. Because a gene's sQTL and eQTL SNPs overlap, the joint
model splits and understates each annotation, and bulk GTEx QTLs under-sample cell-type-specific
splicing, so these enrichments are a floor, not a ceiling.

## The switch axis is productive UTR/CDS remodeling, not decay (S-real-5)

Source: `06_switch_mechanism/_m/SWITCH_CONSEQUENCE_SUMMARY.md` (`figSwitchConsequence`).

Characterizing *what* the switches do structurally, against a within-gene permutation null that
controls transcript count and length, only the productive classes were enriched: UTR change in
10/10 regions (median 1.27×, 47% of switches, Fisher-combined p = 1.1×10⁻¹⁹) and CDS change /
combined coding consequence in 10/10 (median 1.04×, 76%, p = 2.1×10⁻¹⁹). Every degradative
class was never enriched — NMD-status switch (0.92, 0/10), coding-status loss (0.86, 0/10),
biotype switch (0.87, 0/10) — and first/last/internal-exon differences were near-neutral.
Critically the **GO-invisible modules carried the identical signature** (UTR 1.29× vs 1.26×
GO-visible; CDS 1.04× in both), so the pathway-invisible biology is structurally the same
productive remodeling, not an artifact class.

## Co-switch modules share candidate RBP regulons (S-real-6)

Source: `07_rbp_regulation/_m/rbp/RBP_REGULON_SUMMARY.md` (`figRbpRegulon`).

The 3′UTR/CDS remodeling above is the substrate read out by sequence-specific RNA-binding
proteins, motivating a shared-*trans*-factor test for module coordination. Scanning switch-pair
transcripts against ATtRACT human RBP PWMs with MOODS [@doi:10.1093/database/baw035;
@doi:10.1093/bioinformatics/btp554] and testing per-module over-representation of within-pair
binding-site switches, **829 module × RBP pairs were significant at BH q < 0.05 (245 in
GO-invisible modules), spanning 129 distinct RBPs** across 16,973 tests — coordination is broad,
not driven by one factor. The most reproducible regulons were neuronal 3′UTR/splicing factors:
KHDRBS1 (SAM68) in 8/10 regions; A1CF, KHDRBS3, RBMS3, PPIE, RNASEL, and U2AF2 in 7/10; and a
frontal-cortex module (M008) enriched for PPIE (3.22×, q = 1.4×10⁻⁴⁰) and the neuronal ELAV
proteins ELAVL4/3/2. The recurrence of APA/3′UTR regulators (NUDT21, CPEB2/4, ELAV, RBMS)
mirrors the 3′UTR remodeling of the previous section. Motif presence is a computational
prediction, so these are *candidate* regulons for experimental follow-up.

## The switch layer is constrained but non-coding in consequence (S-real-7)

Source: `06_switch_mechanism/_m/CLINICAL_CONSEQUENCE_META.md` (`figClinicalConsequence`).

Switch genes are more loss-of-function constrained than genome-wide in every region (median
gnomAD LOEUF 0.72 vs 0.94; Fisher-combined MWU p = 1.1×10⁻⁹³; GO-invisible 0.700, p = 1.8×10⁻⁹²)
[citation needed: gnomAD LOEUF]. Yet the switched exons themselves carry **lower** ClinVar
pathogenic/likely-pathogenic (P/LP) density than the same genes' constitutive exons in 10/10
regions (median ratio 0.18, robust to a coding-only CDS scope at 0.21) — the expected
alternative-exon biology and a clean statement that the switch consequence is regulatory, not
Mendelian-coding. SNCA instantiates this precisely: its alternative first exon and other 5′
switched exons are non-coding (`cds_overlap = False`) with zero ClinVar P/LP, whereas SNCA's
four P/LP variants sit in a constitutive, isoform-shared coding exon outside the switched region.
The common-variant switch mechanism (isoform choice at a dosage-sensitive gene, SNCA LOEUF 0.40)
is therefore distinct from and complementary to the rare coding variants of Mendelian
synucleinopathy.

## Per-gene deep-dive and literature (Fig 4C–D; Tables S8–S12)

Source: `05_genetic_anchoring/_m/deep_dive/DEEP_DIVE_PANEL.md`, `deep_dive_literature.parquet`.

A deterministic per-gene join (`gene_deep_dive.py`) assembles all six layers into one vignette
per colocalized gene, classifies the verdict (12 splicing-led / 23 splicing-unresolved /
33 expression-led), and attaches a curated literature layer. Four splicing-led genes have
established disease isoform biology that matches their resolved switch: **SNCA**'s multi-5′UTR /
exon-skipping program in synucleinopathy [@doi:10.3389/fgene.2019.00584;
@doi:10.3390/genes9020063]; **DLG1/SAP97**, whose alternatively-spliced synaptic isoforms
include a variant reported down-regulated in early-onset schizophrenia at the 3q29 locus
[@doi:10.1038/tp.2015.154]; **CTSH**, the protective AD locus whose coding change affects only a
subset of transcript isoforms [@doi:10.1038/s41386-023-01542-2]; and **ARVCF**, a 22q11.2 SCZ gene
that is itself a splicing modulator [@doi:10.1038/sj.mp.4001586]. The remaining eight
(PPP6R2, GGNBP2, PGS1, CDIP1, PRRC2B, RTEL1, TBC1D15, TPCN1) lack established disease-specific
isoform literature and stand as **novel splicing-led candidates** — the under-characterized,
GO-invisible switching the method is designed to nominate. Machine-readable per-gene tables
(events S9, RBP regulators S10, exon-clinical S11, literature S12) let a reader reconstruct any
colocalized gene's mechanistic vignette without the hand-written narrative.

---

## Integration & limitations

Colocalization posteriors are modest and bulk-tissue-derived, so per-gene claims are suggestive;
the defensible claim is the **set-level pattern** — splicing-led ≡ GO-invisible, cross-disease
concordance (SNCA), productive-not-degradative remodeling, shared candidate RBP regulons,
LoF-constrained genes with non-coding switch consequence. RBP motif calls are sequence
predictions, not measured binding, and do not yet test whether the lead QTL variant sits inside a
switched-exon motif (the outstanding genomic-intronic-scope extension). This section should be
read against the matched-baseline control of Fig 3 (the sQTL-sparing specificity is IsoGraph-only
on identical WGCNA feature matrices) and the bounding supplement S-real-1 (on per-module rates
IsoGraph is not globally superior): together they place the genetic anchoring as evidence that a
**complementary** layer is real, not that the method dominates.

**Do not merge the two GO-invisible statements in this section with Fig 3's.** The claim here is
about *content*: the 12 splicing-led colocalized genes happen to sit in modules that gene-level
pathway enrichment would miss (n = 12, descriptive, no test against a GO-visible comparator).
Fig 3's set-level sQTL/eQTL contrast does **not** localise to the GO-invisible modules — on the
2026-08-29 refresh GO-invisible (1.068, p = 0.077) and GO-visible (1.084, p = 0.050) are
indistinguishable, and the contrast is carried by the phenotype-associated set (1.111,
p = 3.6e-4). Earlier drafts asserted a genetic GO-invisible localisation; that is retracted, and
the descriptive coloc observation above must not be used to re-import it.
