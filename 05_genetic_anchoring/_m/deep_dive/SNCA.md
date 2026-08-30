# SNCA — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** LBD,PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cortex)
- **Constraint:** LOEUF 0.397  ·  missense o/e 0.71

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | sQTL | Brain_Cortex | 0.04 | yes | yes | yes | no annotated structural change | risk allele A increases usage of junction chr4:89835692-89836127(-) (ENST00000508895.5,ENS |
| pd | sQTL | Brain_Frontal_Cortex_BA9 | 0.03 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr4:89835692-89836127(-) (ENST00000508895.5,ENS |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** DHX9, ZC3H10, ADAR, RALY, G3BP1, RBMS3, RBM41, RBMS1, CPEB2, IFIH1
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ADAR(6), AGO1(6), AGO2(6), AKAP1(6), CELF6(6), CNOT4(6), CPEB2(6), G3BP1(6), DHX9(6), HNRNPA3(6), GRSF1(6), HNRNPLL(6), IFIH1(6), IGHMBP2(6), LIN28A(6), HNRNPAB(6), ZCRB1(6), YTHDC1(6), ZFP36(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for SNCA in LBD,PD — the same switch is genetically anchored across more than one trait. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
SNCA carries an extensively documented alternative-splicing program that is disease-relevant in synucleinopathy: at least four alternative 5'UTR first exons plus internal exon-3/exon-5 skipping generate transcripts that are differentially expressed across PD and dementia-with-Lewy-bodies brain regions, and the coding splice variants (SNCA-126/112/98) modulate alpha-synuclein aggregation kinetics. The IsoGraph-resolved event here is a 5'-end (alternative first exon) choice, matching the well-established 5'UTR/regulatory arm of this program rather than a coding change -- consistent with a dosage mechanism at a LoF-constrained gene (LOEUF 0.40).

_References:_ @doi:10.3389/fgene.2019.00584; @doi:10.3390/genes9020063
