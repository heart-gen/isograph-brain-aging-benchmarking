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
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, RBMS1, RBMS3, G3BP1, SNRPB2, HNRNPA1L2, PPRC1, HNRNPA3, RBM42, FXR1, CPEB1, ENOX1, PABPC3, DHX58, HNRNPD
- **All switched-motif RBPs (recurrence across regions):** AGO2(4), CELF4(4), AKAP1(4), CMTR1(4), CNOT4(4), CELF5(4), CELF6(4), CPEB2(4), EIF4A3(4), DDX19B(4), CSTF2(4), HNRNPA0(4), HNRNPA1L2(4), HNRNPA3(4), DHX58(4), ENOX1(4), G3BP1(4), FXR1(4), ESRP1(4), GRSF1(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TPCN1 in AD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
TPCN1 (endolysosomal two-pore Ca2+ channel; AD locus) colocalizes as a splicing-led switch with no established disease isoform biology -- a novel candidate.
