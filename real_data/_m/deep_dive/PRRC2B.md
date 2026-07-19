# PRRC2B — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF 0.338  ·  missense o/e 0.94

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.02 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr9:131483445-131484686(+) (ENST00000682501.1,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs113866326 | C | 2.670367956161499 | risk allele C increases intron usage of chr9:131483445-131484686(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), CELF6(3), DHX58(3), ENOX1(3), FXR2(3), G3BP2(3), IGF2BP2(3), PABPC3(3), PABPN1(3), PPRC1(3), SAMD4A(3), SRSF4(3), YBX2(3), ZC3H10(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PRRC2B in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
PRRC2B is LoF-constrained (LOEUF 0.34) and colocalizes as a splicing-led switch in schizophrenia with no established disease isoform literature -- a novel candidate.
