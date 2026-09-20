# INO80E — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.282  ·  missense o/e 1.29

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellum | 0.03 | no | — | yes | — | risk allele G increases INO80E expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr16:30001040-30001212(+) (ENST00000540562.1,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr16:30001040-30005221(+) (ENST00000304516.11;  |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | — | yes | — | risk allele G increases usage of junction chr16:30001530-30003436(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr16:30004657-30005221(+) (ENST00000540562.1,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs8054556 | G | 0.4421675503253937 | risk allele G increases expression of INO80E |
| scz | Brain_Cerebellar_Hemisphere | rs4788204 | G | 0.34582898020744324 | risk allele G increases intron usage of chr16:30001040-30001212(+) |
| scz | Brain_Cerebellar_Hemisphere | rs4788204 | G | -0.2865675091743469 | risk allele G decreases intron usage of chr16:30001040-30005221(+) |
| scz | Brain_Cerebellar_Hemisphere | rs4788204 | G | 0.6303684711456299 | risk allele G increases intron usage of chr16:30001530-30003436(+) |
| scz | Brain_Cerebellar_Hemisphere | rs4788204 | G | 0.46899572014808655 | risk allele G increases intron usage of chr16:30004657-30005221(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, PPRC1, ENOX1, CNOT4, CELF4, CELF5, PABPC3, EIF4B, A1CF, CPEB2, ZFP36L2, ZNF638, RBM24, ESRP2, RALY
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), CELF1(2), CELF2(2), CELF4(2), CELF5(2), CMTR1(2), CNOT4(2), CPEB1(2), CPEB2(2), CPEB4(2), CSTF2(2), DDX19B(2), DHX58(2), EIF4A3(2), EIF4B(2), ELAVL1(2), ELAVL3(2), ENOX1(2), ESRP1(2), ESRP2(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for INO80E in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
INO80E is a subunit of the INO80 chromatin-remodelling complex and lies in the 16p11.2 schizophrenia locus. Its genetics here colocalize on both modalities and are not splicing-specific; isoform biology is not established.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
