# PCGF3 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** AD,PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.445  ·  missense o/e 0.51

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | no | — | risk allele A increases usage of junction chr4:705970-730630(+) (ENST00000362003.10,ENST00 |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | no | — | risk allele A decreases usage of junction chr4:725229-730630(+) (ENST00000433814.5; not in |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele A decreases usage of junction chr4:728628-730630(+) (unmapped transcript; not  |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | no | — | risk allele A increases usage of junction chr4:731110-733672(+) (ENST00000362003.10,ENST00 |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | no | — | risk allele A decreases usage of junction chr4:732498-733672(+) (ENST00000419774.5,ENST000 |
| pd | sQTL | Brain_Cerebellum | 0.01 | yes | yes | no | biotype switch, cds, first exon, internal exon,  | risk allele G increases usage of junction chr4:731110-733672(+) (ENST00000362003.10,ENST00 |
| pd | sQTL | Brain_Cerebellum | 0.01 | yes | yes | no | biotype switch, cds, first exon, internal exon,  | risk allele G increases usage of junction chr4:761416-764984(+) (ENST00000362003.10,ENST00 |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellar_Hemisphere | rs13116048 | A | 0.6466434597969055 | risk allele A increases intron usage of chr4:705970-730630(+) |
| ad | Brain_Cerebellar_Hemisphere | rs13116048 | A | -0.4211387038230896 | risk allele A decreases intron usage of chr4:725229-730630(+) |
| ad | Brain_Cerebellar_Hemisphere | rs13116048 | A | -0.4404391050338745 | risk allele A decreases intron usage of chr4:728628-730630(+) |
| ad | Brain_Cerebellar_Hemisphere | rs13116048 | A | 0.44158050417900085 | risk allele A increases intron usage of chr4:731110-733672(+) |
| ad | Brain_Cerebellar_Hemisphere | rs13116048 | A | -0.45222926139831543 | risk allele A decreases intron usage of chr4:732498-733672(+) |
| pd | Brain_Cerebellum | rs2242237 | G | 0.4140298068523407 | risk allele G increases intron usage of chr4:731110-733672(+) |
| pd | Brain_Cerebellum | rs2242237 | G | 0.42370930314064026 | risk allele G increases intron usage of chr4:761416-764984(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, RBM14, ZNF638, RBMS3, SNRPB2, PABPC4, G3BP2, ENOX1, ZC3H10, YBX2, RBMY1A1, HNRNPLL, AKAP1, TIA1, CELF5
- **All switched-motif RBPs (recurrence across regions):** ACO1(8), AGO1(8), AGO2(8), AKAP1(8), CELF4(8), CELF5(8), CELF6(8), CNOT4(8), CPEB2(8), ESRP1(8), CSTF2(8), ZRANB2(8), YBX2(8), ENOX1(8), HNRNPA3(8), HNRNPA0(8), G3BP2(8), FXR2(8), FXR1(8), HNRNPCL1(8)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PCGF3 in PD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
PCGF3 is a Polycomb RING-finger subunit of a non-canonical PRC1 complex; no disease-specific isoform biology is established, and its usage range here touches the detection floor.

_Curation: novel_candidate._ No established disease-specific isoform biology was found. That is a statement about the literature, not about the evidence here: an under-characterised switch is what this method is built to surface, and this row is a nomination rather than a null result.
