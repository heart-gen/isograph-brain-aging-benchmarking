# DLG1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.441  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | no | no | — | risk allele T increases usage of junction chr3:197194589-197225857(-) (ENST00000469073.2,E |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | no | no annotated structural change | risk allele T decreases usage of junction chr3:197194589-197282679(-) (ENST00000346964.6,E |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | no | no annotated structural change | risk allele T increases usage of junction chr3:197282845-197296346(-) (ENST00000346964.6,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs9843908 | T | 0.5421736836433411 | risk allele T increases intron usage of chr3:197194589-197225857(-) |
| scz | Brain_Cerebellar_Hemisphere | rs9843908 | T | -0.3680476248264313 | risk allele T decreases intron usage of chr3:197194589-197282679(-) |
| scz | Brain_Cerebellar_Hemisphere | rs9843908 | T | 0.31144991517066956 | risk allele T increases intron usage of chr3:197282845-197296346(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, SNRPB2, PPRC1, G3BP2, ENOX1, CNOT4, ZC3H10, HNRNPA3, CELF4, CELF5, FXR1, CELF6, PABPC5
- **All switched-motif RBPs (recurrence across regions):** ACO1(8), CELF6(8), CNOT4(8), DHX58(8), DDX58(8), PPRC1(8), TRA2A(8), NELFE(8), HNRNPU(8), HNRNPA3(8), FXR1(8), YTHDC1(7), G3BP2(7), ZC3H10(7), RBM41(7), SNRPB2(6), SUPV3L1(6), AGO2(6), RBFOX2(6), ZFP36L2(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for DLG1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
DLG1/SAP97 is a canonical alternatively-spliced synaptic scaffold: N-terminal alpha vs beta isoforms, an internal I3 insert and additional cassette exons tune its PDZ/GK synaptic function. A DLG1 splice variant is reported to be expressed at reduced cortical levels in early-onset schizophrenia, and DLG1 sits in the 3q29 schizophrenia locus, so an sQTL that shifts DLG1 isoform choice is a mechanistically plausible splicing-led route to SCZ risk.

_References:_ @doi:10.1038/tp.2015.154
