# TMEM175 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** AD,ALS,PD  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.48 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.295  ·  missense o/e 1.04

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases TMEM175 expression (gene-level; no intron) |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.48 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr4:932540-947709(+) (ENST00000264771.9,ENST000 |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.48 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr4:932540-950421(+) (ENST00000508204.5,ENST000 |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.48 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr4:947892-948116(+) (ENST00000264771.9,ENST000 |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.48 | no | no | yes | — | risk allele C increases usage of junction chr4:948615-950424(+) (ENST00000504180.5,ENST000 |
| pd | sQTL | Brain_Cortex | 0.05 | yes | yes | yes | no annotated structural change | risk allele A decreases usage of junction chr4:932540-947709(+) (ENST00000264771.9,ENST000 |
| pd | sQTL | Brain_Cortex | 0.05 | no | no | yes | — | risk allele A increases usage of junction chr4:948615-950424(+) (ENST00000504180.5,ENST000 |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | no | yes | — | risk allele C decreases usage of junction chr4:932540-945999(+) (ENST00000504850.1,ENST000 |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | no | yes | — | risk allele C decreases usage of junction chr4:946126-947709(+) (ENST00000504850.1,ENST000 |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | no | yes | — | risk allele C increases usage of junction chr4:951258-952367(+) (ENST00000509508.5; not in |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | — | yes | — | risk allele C decreases usage of junction chr4:951258-953190(+) (unmapped transcript; not  |
| ad | sQTL | Brain_Cerebellum | 0.03 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr4:951717-952367(+) (ENST00000264771.9,ENST000 |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | — | yes | — | risk allele C decreases usage of junction chr4:951717-953190(+) (unmapped transcript; not  |
| ad | sQTL | Brain_Cerebellum | 0.03 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr4:952450-953190(+) (ENST00000264771.9,ENST000 |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | — | yes | — | risk allele C increases usage of junction chr4:952539-953190(+) (unmapped transcript; not  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellum | rs11552301 | C | -0.7777791023254395 | risk allele C decreases intron usage of chr4:932540-945999(+) |
| ad | Brain_Cerebellum | rs11552301 | C | -0.6660060882568359 | risk allele C decreases intron usage of chr4:946126-947709(+) |
| ad | Brain_Cerebellum | rs11552301 | C | 0.8619697093963623 | risk allele C increases intron usage of chr4:951258-952367(+) |
| ad | Brain_Cerebellum | rs11552301 | C | -0.8472946882247925 | risk allele C decreases intron usage of chr4:951258-953190(+) |
| ad | Brain_Cerebellum | rs11552301 | C | 1.2017526626586914 | risk allele C increases intron usage of chr4:951717-952367(+) |
| ad | Brain_Cerebellum | rs11552301 | C | -1.3391225337982178 | risk allele C decreases intron usage of chr4:951717-953190(+) |
| ad | Brain_Cerebellum | rs11552301 | C | -0.5487470626831055 | risk allele C decreases intron usage of chr4:952450-953190(+) |
| ad | Brain_Cerebellum | rs11552301 | C | 1.2641254663467407 | risk allele C increases intron usage of chr4:952539-953190(+) |
| als | Brain_Cerebellar_Hemisphere | rs873786 | C | -0.7011362314224243 | risk allele C decreases intron usage of chr4:932540-947709(+) |
| als | Brain_Cerebellar_Hemisphere | rs873786 | C | 0.6464655995368958 | risk allele C increases intron usage of chr4:932540-950421(+) |
| als | Brain_Cerebellar_Hemisphere | rs873786 | C | -0.8334869742393494 | risk allele C decreases intron usage of chr4:947892-948116(+) |
| als | Brain_Cerebellar_Hemisphere | rs873786 | C | 0.7598555684089661 | risk allele C increases intron usage of chr4:948615-950424(+) |
| pd | Brain_Cerebellum | rs73211813 | T | -0.3880670964717865 | risk allele T decreases expression of TMEM175 |
| pd | Brain_Cortex | rs77060135 | A | -1.1121370792388916 | risk allele A decreases intron usage of chr4:932540-947709(+) |
| pd | Brain_Cortex | rs77060135 | A | 0.7959812879562378 | risk allele A increases intron usage of chr4:948615-950424(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPA0, ELAVL3, KHDRBS1, HNRNPD, ELAVL4, CPEB1, CPEB4, TIAL1, U2AF2, PPIE, DAZAP1, AGO1, SF1, ZRANB2
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO1(3), CPEB1(3), CPEB2(3), CPEB4(3), DAZAP1(3), DDX19B(3), EIF4A3(3), EIF4B(3), ELAVL3(3), ELAVL4(3), G3BP1(3), HNRNPA0(3), HNRNPA2B1(3), HNRNPCL1(3), HNRNPD(3), HNRNPK(3), HNRNPU(3), KHDRBS1(3), KHDRBS2(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TMEM175 in AD,ALS,PD — the same switch is genetically anchored across more than one trait. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
