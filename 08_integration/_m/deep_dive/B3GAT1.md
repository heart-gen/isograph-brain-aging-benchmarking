# B3GAT1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.01 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.749  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Caudate_basal_ganglia | 0.01 | no | — | yes | — | risk allele A decreases B3GAT1 expression (gene-level; no intron) |
| scz | sQTL | Brain_Cortex | 0.01 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr11:134380747-134381924(-) (ENST00000312527.9, |
| scz | sQTL | Brain_Cortex | 0.01 | yes | yes | yes | no annotated structural change | risk allele G decreases usage of junction chr11:134381817-134381924(-) (ENST00000524765.1; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs1440480 | A | -0.19964297115802765 | risk allele A decreases expression of B3GAT1 |
| scz | Brain_Cortex | rs61908686 | G | 1.265340805053711 | risk allele G increases intron usage of chr11:134380747-134381924(-) |
| scz | Brain_Cortex | rs61908686 | G | -1.1389235258102417 | risk allele G decreases intron usage of chr11:134381817-134381924(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PABPN1, PUM1, YBX2, RNASEL, PUM2, CELF6, RBM25, PABPC4, RBMS3, IGF2BP2, ELAVL4, CELF4, KHDRBS1, ZNF638, CELF5
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), AGO2(4), AKAP1(4), CELF4(4), CELF5(4), CELF6(4), CPEB1(4), CPEB4(4), CSTF2(4), DAZAP1(4), DDX19B(4), DHX58(4), EIF4A3(4), EIF4B(4), ELAVL3(4), ELAVL4(4), ESRP1(4), ESRP2(4), G3BP1(4), HNRNPA0(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for B3GAT1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
