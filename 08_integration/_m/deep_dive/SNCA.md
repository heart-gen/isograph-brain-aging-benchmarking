# SNCA — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** LBD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cortex)
- **Constraint:** LOEUF 0.397  ·  missense o/e 0.71

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | sQTL | Brain_Cortex | 0.04 | no | no | yes | — | risk allele A increases usage of junction chr4:89835692-89836127(-) (ENST00000508895.5,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| lbd | Brain_Cortex | rs7680557 | A | 0.7326975464820862 | risk allele A increases intron usage of chr4:89835692-89836127(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ADAR(1), AGO1(1), AGO2(1), AKAP1(1), CELF4(1), CELF6(1), CNOT4(1), CPEB1(1), CPEB2(1), CPEB4(1), CSTF2(1), DDX19B(1), DHX9(1), EIF4A3(1), ELAVL1(1), ELAVL3(1), ELAVL4(1), FXR2(1), G3BP1(1)

## 5. Interpretation
SNCA has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.

## 6. Literature (known isoform biology)
SNCA carries an extensively documented alternative-splicing program that is disease-relevant in synucleinopathy: at least four alternative 5'UTR first exons plus internal exon-3/exon-5 skipping generate transcripts differentially expressed across PD and dementia-with-Lewy-bodies brain regions, and the coding splice variants (SNCA-126/112/98) modulate alpha-synuclein aggregation kinetics. The event IsoGraph resolved here is a 5'-end (alternative first exon) choice, matching the regulatory arm of that program rather than a coding change. Not in the current anchored set: SNCA is no longer concordant at resolution 2.0.

_Curation: isoform_documented._ A disease-relevant isoform program is established for this gene.

_References:_ @doi:10.3389/fgene.2019.00584; @doi:10.3390/genes9020063
