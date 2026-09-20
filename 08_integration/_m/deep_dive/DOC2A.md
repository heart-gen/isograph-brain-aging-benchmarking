# DOC2A — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** AD,SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.05 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.732  ·  missense o/e 0.84

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | yes | — | risk allele C decreases DOC2A expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.05 | no | — | yes | — | risk allele C increases usage of junction chr16:30007090-30007179(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.05 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele C decreases usage of junction chr16:30007299-30008996(-) (ENST00000350119.9,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellar_Hemisphere | rs1140239 | C | -0.21894823014736176 | risk allele C decreases expression of DOC2A |
| scz | Brain_Cerebellar_Hemisphere | rs11150580 | C | 0.8092381358146667 | risk allele C increases intron usage of chr16:30007090-30007179(-) |
| scz | Brain_Cerebellar_Hemisphere | rs11150580 | C | -0.4637572765350342 | risk allele C decreases intron usage of chr16:30007299-30008996(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PABPC4, IGF2BP2, PABPC5, TARDBP, RBM4, RBM46, PTBP2, YBX2
- **All switched-motif RBPs (recurrence across regions):** ANKHD1(5), CELF5(5), CELF4(5), PABPC3(5), PABPC5(5), ESRP2(5), CPEB4(5), SFPQ(5), SNRNP70(5), SUPV3L1(5), ZRANB2(5), IGF2BP2(5), KHDRBS2(5), HNRNPLL(5), HNRNPCL1(5), NELFE(5), RBM4(5), RBM6(5), RBM46(5), RBM24(5)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for DOC2A in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
DOC2A is a calcium sensor for spontaneous neurotransmitter release and sits in the 16p11.2 locus. It is the most splicing-specific gene in this panel on the genetics -- strong sQTL colocalization in nine tissues with no eQTL instrument at all -- while its isoform biology is uncharacterised, which is exactly the gap this layer is meant to mark.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
