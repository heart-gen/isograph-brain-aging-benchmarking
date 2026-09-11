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

## Signal-level colocalization nominates candidate loci, and an event audit bounds them

Source: `05_genetic_anchoring/_m/coloc_signal_susie/all_introns/COLOC_SIGNAL_SUSIE.md`,
`05_genetic_anchoring/_m/locus_event_audit/susie_all_introns/LOCUS_EVENT_AUDIT.md`
(regenerated 2026-09-11). *This section is the primary locus-nomination layer. The eCAVIAR (CLPP)
section that follows is retained as orthogonal sensitivity evidence and must not be quoted as
signal-level.*

eCAVIAR CLPP multiplies credible-set inclusion probabilities, so it cannot represent a locus
carrying more than one causal signal and yields individually modest posteriors. We therefore
re-fine-mapped every GTEx v11 QTL from the all-pairs release with `susie_rss` on the GWAS's own
1000 Genomes European LD and colocalized signal against signal with `coloc.susie`
[@doi:10.1371/journal.pgen.1009440], keeping a re-fitted QTL credible set only where it agrees with
GTEx's own in-sample credible set and falling back to `coloc.abf`
[@doi:10.1371/journal.pgen.1004383] where either side does not fine-map. The fallback is the
majority estimator and is reported as such: only 230 of 579 GWAS loci fine-map, and 2,933 of
42,130 (cell, QTL class) rows — 7.0% — are scored signal against signal. Testing every intron
phenotype of each gene (the all-introns arm), **42 gene × trait pairs (40 genes) reach
PP4_sQTL ≥ 0.8** at the prespecified prior $p_{12} = 10^{-5}$; 17 are headlined by `coloc.susie`
and 25 by `coloc.abf`. Across the $p_{12}$ sweep, 7 calls hold at $10^{-6}$, 20 from
$5\times10^{-6}$, and 15 only from the primary prior.

Fourteen of the 42 are **sQTL-preferential** (PP4_eQTL < 0.5). That count is a selection — a
strong sQTL posterior is required first — so it nominates loci rather than estimating splicing
specificity, and a low eQTL posterior can equally reflect no eQTL, low power, a different causal
architecture or multiple signals. The unbiased within-gene comparison, run on GTEx's
representative intron so each gene contributes one phenotype per QTL class, shows no splicing
preference (21 splicing-only vs 43 expression-only discordant genes, exact McNemar P = 0.008 in
the expression direction; conditional-posterior Wilcoxon P = 0.33).

**SNCA** is the strongest signal-level nomination. In LBD, PP4_sQTL is 0.974 against PP4_eQTL
0.113, reaching the call in 11 of 13 tissues, with SuSiE resolving two QTL signals of which only
one colocalizes — the multi-signal architecture `coloc.abf` cannot represent. The PD signal
colocalizes with the same intron (chr4:89,835,692–89,836,127; PP4_sQTL 0.974, 8 tissues). Both
calls sit just below 0.8 at $p_{12} = 10^{-6}$ (PP4 0.79).

**A posterior names a gene, not an event.** An event audit walked each nomination from GWAS signal
to sQTL signal to intron phenotype to driver transcript, and compared the colocalizing intron with
curated, coordinate-verified literature events. No locus recovered a known mechanism. The *UNC13A*
ALS signal (PP4_sQTL 0.961; cerebellum and cerebellar hemisphere) colocalizes at
chr19:17,630,750–17,632,782, ~9 kb from the intron harbouring the TDP-43-dependent cryptic exon
[@doi:10.1038/s41586-022-04424-7; @doi:10.1038/s41586-022-04436-3]. GTEx carries an sQTL phenotype
at that intron in exactly those two tissues and it did not colocalize, so the signal is
context-distinct from the known event rather than a near miss. *PICALM* (AD) colocalizes only
under `coloc.abf` on the primary grid (PP4_sQTL 0.818, cortex), because its locus exceeds the
12,000-SNP fine-mapping limit, and at an intron distinct from the transcriptome-wide-association
event reported at the locus [@doi:10.1038/s41588-018-0238-1]. In a prespecified 30,000-SNP
sensitivity arm it colocalizes at signal level on the same intron (PP4 0.813) with an LD-robust
GWAS fit, but the call falls to 0.30 at $p_{12} = 10^{-6}$. The remaining 40 nominations have no
curated event and stand as novel splice-led candidates.

**The signal-level anchored switches behave like switches in orthogonal long-read data.** For 64
events in 24 nominated genes, the colocalizing intron maps onto a transcript of the gene's
tissue-matched IsoGraph switch pair. Scored in ONT long-read DLPFC (n = 12)
[@doi:10.1038/s41587-024-02245-9] against IsoGraph switch pairs from other genes matched on
abundance decile, 0.565 of the 124 detected anchored pairs are switch-like against a matched-null
mean of 0.238 (95% null interval 0.169–0.315; empirical P = 5 × 10⁻⁴), and 0.672 against 0.298
among the 61 pairs whose anchored isoform is usably expressed. The mean usage correlation itself
does not separate from its null (−0.104 vs −0.131, P = 0.50), so the evidence is the rate of
switch-like pairs read against the compositional background, not the strength of
anti-correlation. Six genes, SNCA and CTSH among them, never reach usable abundance for the
anchored isoform, so the confirmation holds for the set, not locus by locus. UNC13A and PICALM
cannot enter this test: neither has an IsoGraph switch pair in the tissue where it colocalizes.

**SMR agrees with most nominations, as expected, and disagrees at one informative locus.** SMR
and HEIDI [@doi:10.1038/ng.3538] were run on each nomination's colocalizing intron in every tissue
where its sQTL call holds, using the same GWAS and GTEx summary statistics as colocalization, so
agreement is a consistency check rather than independent evidence. Of the 42 nominations, 29 have
at least one tissue in which the colocalizing intron is SMR-significant (Bonferroni within trait
and QTL class) without HEIDI rejecting a single shared variant ($p_{HEIDI} \geq 0.01$), including
SNCA in LBD (11 of 11 tissues), UNC13A (2 of 2) and PICALM (1 of 1); 9 have no cis-QTL instrument
at $p < 5 \times 10^{-8}$ in any tissue. The disagreement is SNCA in PD: the same intron colocalizes
in 8 tissues, yet SMR is not significant after correction in 7 (minimum $p_{SMR}$ 1.1 × 10⁻³) and
HEIDI rejects in the eighth (median $p_{HEIDI}$ 7 × 10⁻¹¹). A HEIDI rejection does not overrule a
signal-level colocalization, but the PD arm of the shared SNCA event is supported by
colocalization alone, and the disagreement is reported as such. The HEIDI cut matters: of the 113
primary tissue probes not rejected at 0.01, 20 would be rejected at 0.05. No `b_SMR` is read as
causal direction.

## Disease variants colocalize onto GO-invisible isoform switches (eCAVIAR CLPP layer)

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
change. **This exact contrast is confirmed orthogonally in BrainSEQ short-read**
(`06_switch_mechanism/_m/junction_coloc_confirm/`): the LIBD PSI event measuring the
anchored proximal exon against the canonical distal one gives a minor-form usage of
**0.189 in DLPFC (n = 222) and 0.234 in caudate (n = 238)** — both forms carry substantial
usage, so the switch is real and used in the regions the colocalization was found in.
An earlier ONT long-read check put the anchored isoform at 0.29% of gene output and failed
to confirm it; short-read junction data, which measures the junction the sQTL actually tags
rather than a whole-transcript proxy at ~20x the long-read n, places it two orders of
magnitude higher. The long-read failure is therefore an assay limitation, not a refutation.
**CTSH does not confirm**: its anchored acceptor reaches only 0.016–0.020 minor-form usage
in hippocampus, below the pre-registered 0.05 threshold, so CTSH stays off any main figure. eCAVIAR posteriors are modest (max 0.39, CTSH in AD hippocampus; most < 0.1), so the
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
binding-site switches, **714 module × RBP pairs were significant at BH q < 0.05 (149 in
GO-invisible modules), spanning 138 distinct RBPs of 160 tested** across 34,240 module × RBP
cells in 14 regions — coordination is broad, not driven by one factor. Under the
opportunity-adjusted GLM (transcript length, GC, UTR composition) 310 cells are significant,
but **only 43 overlap the hypergeometric set**, so the adjustment reshuffles the ranking rather
than thinning it; the adjusted arm is the one to quote for a regulon claim. The most
reproducible regulons were neuronal 3′UTR/splicing factors:
KHDRBS1 (SAM68), ELAVL4 and PPRC1 in 8/14 regions; CPEB4, RBMS3, RNASEL, ELAVL3, RBM14 and
ADAR in 7/14; and a
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

## Two negative controls bound the genetic claim

Source: `05_genetic_anchoring/_m/coloc_modality_contrast/` (all four arms) and
`05_genetic_anchoring/_m/module_coloc_convergence/MODULE_COLOC_CONVERGENCE.md` (Fig 4E, S-real-10).

The set-level contrast of Fig 3 says that, *as a set*, phenotype-associated switch modules are
spared at splicing QTL relative to expression QTL. Two pre-specified tests asked whether that
survives at finer resolution. Neither does, and we report both.

**At per-gene resolution there is no modality preference attributable to module membership.**
We ran `coloc.abf` on GTEx v11 all-pairs for every gene at each GWAS locus, once against the
gene's sQTL and once against its eQTL — 126,390 fits over 13 brain tissues [citation needed:
coloc]. Because the two modalities are compared *within* a gene, module membership cancels, so
the only channel by which a module-detection method can move the statistic is which genes it
selects. We therefore ran four gene pools: IsoGraph switch genes, all testable genes at the same
loci, and the two matched-feature WGCNA baselines. The splicing share of discordant genes is
**~31% in every pool** — IsoGraph switch 20 splicing-only vs 46 expression-only (0.303),
background 108/234 (0.316), WGCNA-switch 59/129 (0.314), WGCNA-multiplex 95/209 (0.313) — and
the apparent significance of the expression-favouring skew tracks pool size alone (exact McNemar
P = 1.9 × 10⁻³ to 5.5 × 10⁻¹¹). The decisive comparison is locus-matched: within the background
pool, where both groups sit at the *same* loci under the same estimator, switch genes give 20/46
and non-switch genes 88/188 (**Fisher exact P = 0.88**). On the power-corrected conditional
posterior PP4/(PP3+PP4) the switch arm is null (0.2415 vs 0.2383, Wilcoxon P = 0.39) while the
background is significantly splicing-leaning (0.2575 vs 0.2413, P = 9.6 × 10⁻¹¹) — that is,
non-switch genes look *more* splicing-skewed, not less. The expression-favouring skew is
therefore GTEx eQTL discovery power, not biology, and **this test does not corroborate the
set-level contrast in either direction**. Fig 3 and this test are different estimands — a
set-level ratio of odds ratios against matched baselines, versus a within-gene modality
comparison — and only the former supports the splicing-specificity claim.

**Colocalizing genes do not concentrate in particular modules, in any trait (Fig 4E).** Against
a size-matched permutation null that holds module sizes fixed, no (trait, cohort) cell is more
concentrated than chance (P = 0.19–1.00, 10/10 cells), and no individual module survives FDR
(min q = 0.23 over 46 testable rows). This retracts a convergence result we previously reported:
SCZ colocalized genes appeared enriched in MAGMA-anchored age-sensitive modules (15/31 = 48% vs
a 25% background, P = 0.0095 for GTEx caudate) only because the denominator was *all* module
genes. A gene can colocalize only if it was coloc-tested, and anchored modules are defined by
MAGMA enrichment for the same GWAS that decides which genes enter the test, so the tested pool
is already anchored-rich; against it the same counts are null (**P = 0.30**). The convergence
claim is withdrawn. The module-disruption and candidate-regulator layers of that analysis do not
depend on the denominator and stand.

Both nulls narrow the claim in the same direction: the genetic anchoring is a property of the
switch layer **as a set**, not of individual switch genes or of individual modules.

## Integration & limitations

Colocalization posteriors are modest and bulk-tissue-derived, so per-gene claims are suggestive;
the defensible claim is the **set-level pattern** — splicing-led ≡ GO-invisible, cross-disease
concordance (SNCA), productive-not-degradative remodeling, shared candidate RBP regulons,
LoF-constrained genes with non-coding switch consequence. RBP motif calls are sequence
predictions, not measured binding, and do not yet test whether the lead QTL variant sits inside a
switched-exon motif (the outstanding genomic-intronic-scope extension). This section should be
read against the matched-baseline control of Fig 3 (the sQTL-sparing specificity is IsoGraph-only
on identical WGCNA feature matrices), the two negative controls above (the per-gene modality test
and the module-convergence test are both null), and the bounding supplement S-real-1 (on
per-module rates IsoGraph is not globally superior): together they place the genetic anchoring as
evidence that a **complementary** layer is real, at set level, not that the method dominates or
that individual switch genes and modules are separately anchored.

**Do not merge the two GO-invisible statements in this section with Fig 3's.** The claim here is
about *content*: the 12 splicing-led colocalized genes happen to sit in modules that gene-level
pathway enrichment would miss (n = 12, descriptive, no test against a GO-visible comparator).
Fig 3's set-level sQTL/eQTL contrast does **not** localise to the GO-invisible modules — on the
2026-08-29 refresh GO-invisible (1.068, p = 0.077) and GO-visible (1.084, p = 0.050) are
indistinguishable, and the contrast is carried by the phenotype-associated set (1.111,
p = 3.6e-4). Earlier drafts asserted a genetic GO-invisible localisation; that is retracted, and
the descriptive coloc observation above must not be used to re-import it.
