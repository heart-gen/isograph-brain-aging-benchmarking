# TPCN1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.558  ·  missense o/e 0.79

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele G decreases usage of junction chr12:113284771-113284869(+) (unmapped transcrip |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr12:113284771-113285889(+) (ENST00000335509.11 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM42, CPEB2, YBX2, RALY, CPEB1, ELAVL3, PPIE, DHX58, SYNCRIP, KHDRBS2, DDX19B, KHDRBS1, PPRC1, U2AF2, KHDRBS3
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CMTR1(3), CPEB2(3), DHX58(3), CSTF2(3), EIF4A3(3), HNRNPA1L2(3), FXR1(3), ESRP1(3), IGHMBP2(3), HNRNPLL(3), HNRNPD(3), HNRNPA0(3), HNRNPAB(3), RBM25(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TPCN1 in AD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
