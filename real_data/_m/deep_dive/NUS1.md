# NUS1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF n/a  ·  missense o/e n/a

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C increases NUS1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs80196932 | C | 0.24534142017364502 | risk allele C increases expression of NUS1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** —

## 5. Interpretation
NUS1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
