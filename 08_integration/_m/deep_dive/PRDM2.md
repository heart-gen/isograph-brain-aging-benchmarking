# PRDM2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · replicates in BrainSeq

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.06 (Brain_Cortex)
- **Constraint:** LOEUF 0.185  ·  missense o/e 0.90

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cortex | 0.06 | no | — | no | — | risk allele A increases usage of junction chr1:13816570-13821649(+) (unmapped transcript;  |
| als | sQTL | Brain_Cortex | 0.06 | yes | yes | no | cds, first exon, internal exon, last exon, utr | risk allele A decreases usage of junction chr1:13816570-13823159(+) (ENST00000235372.11,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cortex | rs2744682 | A | 0.9042878746986389 | risk allele A increases intron usage of chr1:13816570-13821649(+) |
| als | Brain_Cortex | rs2744682 | A | -1.0998578071594238 | risk allele A decreases intron usage of chr1:13816570-13823159(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBM6, SNRPB2, RBMS3, PPRC1, G3BP2, ENOX1, CNOT4, ZC3H10, HNRNPA3, SNRNP70, PABPC3, RBM46
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), CPEB1(6), CNOT4(6), CPEB2(6), CPEB4(6), EIF4B(6), DDX19B(6), G3BP2(6), FXR2(6), EIF4A3(6), DHX58(6), ELAVL3(6), ESRP2(6), ESRP1(6), ENOX1(6), HNRNPLL(6), HNRNPD(6), HNRNPCL1(6), HNRNPC(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PRDM2 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage; the switch also replicates in an independent BrainSeq cohort.

## 6. Literature (known isoform biology)
PRDM2/RIZ is the clearest documented isoform program in the anchored set: it is transcribed from alternative promoters into RIZ1, which carries the PR (SET-like) methyltransferase domain, and RIZ2, which lacks it, and the RIZ1:RIZ2 balance -- not total PRDM2 -- is what changes in disease states where the gene has been studied. The colocalizing ALS event here is a 5'-end choice, which is the same class of event as that switch, so the variant plausibly acts by selecting between PR-positive and PR-negative products rather than by changing gene dosage. Brain-specific and ALS-specific isoform work is not established.

_Curation: isoform_documented._ A disease-relevant isoform program is established for this gene.

_References:_ [citation needed: RIZ1/RIZ2 alternative-promoter isoforms and the PR-domain balance]
