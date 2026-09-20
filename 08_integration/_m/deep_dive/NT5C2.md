# NT5C2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.128  ·  missense o/e 0.65

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele T decreases NT5C2 expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellum | 0.02 | yes | yes | no | biotype switch, cds, first exon, internal exon,  | risk allele G increases usage of junction chr10:103174982-103181185(-) (ENST00000404739.8, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs11191580 | T | -0.4148642420768738 | risk allele T decreases expression of NT5C2 |
| scz | Brain_Cerebellum | rs12412038 | G | 1.4276000261306763 | risk allele G increases intron usage of chr10:103174982-103181185(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM6, PPRC1, RBM46, RBMS1, IGF2BP2, PABPC4, DHX9, PABPN1, ADAR, SRSF11, RBM8A, CSTF2, ESRP1, DHX58
- **All switched-motif RBPs (recurrence across regions):** PABPN1(6), ESRP1(6), HNRNPA1L2(6), RBM8A(6), RBM6(6), RBM14(6), SRSF11(6), PPRC1(6), TARDBP(5), AGO2(5), DHX58(3), CSTF2(1), ADAR(1), IGF2BP2(1), DHX9(1), PABPC4(1), LIN28A(1), RBM46(1), RBMS1(1)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for NT5C2 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
NT5C2 is a cytosolic 5'-nucleotidase, an established schizophrenia GWAS gene, and the cause of hereditary spastic paraplegia SPG45 when lost. Isoform choice in brain is not characterised.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
