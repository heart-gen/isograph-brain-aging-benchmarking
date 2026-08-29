# MRPS10 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible · replicates in BrainSeq

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.336  ·  missense o/e 1.08

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr6:42214344-42217802(-) (ENST00000053468.4,ENS |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ESRP1(4), SNRNP70(4)

## 5. Interpretation
MRPS10 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
