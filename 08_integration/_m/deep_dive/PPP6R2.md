# PPP6R2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.06 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.766  ·  missense o/e 0.90

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.06 | yes | yes | no | no annotated structural change | risk allele G decreases usage of junction chr22:50419462-50422254(+) (ENST00000216061.9,EN |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.06 | no | — | no | — | risk allele G increases usage of junction chr22:50419462-50422259(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs76300267 | G | -1.8367607593536377 | risk allele G decreases intron usage of chr22:50419462-50422254(+) |
| scz | Brain_Frontal_Cortex_BA9 | rs76300267 | G | 1.7728859186172485 | risk allele G increases intron usage of chr22:50419462-50422259(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM6, G3BP2, ENOX1, HNRNPA3, CELF4, CELF5, CELF6, MATR3, PABPC5, EIF4B, A1CF, SAMD4A, ZNF638, IGF2BP2, RBM42
- **All switched-motif RBPs (recurrence across regions):** A1CF(7), AGO1(7), CELF4(7), CELF5(7), CELF6(7), CSTF2(7), DHX58(7), EIF4B(7), ENOX1(7), ESRP1(7), FXR2(7), G3BP2(7), HNRNPA3(7), HNRNPM(7), IGF2BP2(7), LIN28A(7), MATR3(7), MSI1(7), NELFE(7), NUDT21(7)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PPP6R2 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
PPP6R2 (PP6 regulatory subunit) colocalizes as a splicing-led switch across both ALS and SCZ (three resolved events, the most in the panel), but disease-specific isoform biology is not established -- a novel cross-trait splicing-led candidate for follow-up.
