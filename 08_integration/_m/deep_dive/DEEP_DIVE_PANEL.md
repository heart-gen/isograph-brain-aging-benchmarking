# Per-gene mechanistic deep-dive — panel summary

Every colocalized disease gene, one row each (ranked splicing-led first, multi-locus above single). `resolved events` = colocalizing sQTLs that map onto a concordant IsoGraph switch pair (the splicing-led, DTU-without-DGE class); `multi-locus` flags genes with >=2 such events (multiple significant colocalizations, incl. cross-trait). Main-figure framing is reserved for the resolved cross-disease headliners (SNCA, CTSH); the remainder are supporting vignettes. See `<GENE>.md` for each.

| gene | traits | kinds | max CLPP | LOEUF | resolved events | multi-locus | concordant traits | BrainSeq rep | GO-inv | verdict |
|------|--------|-------|----------|-------|-----------------|-------------|-------------------|--------------|--------|---------|
| PPP6R2 | ALS,SCZ | sQTL | 0.07 | 0.77 | 3 | yes | ALS,SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| SNCA | LBD,PD | sQTL | 0.04 | 0.40 | 2 | yes | LBD,PD | no | yes | splicing-led (IsoGraph-resolved switch) |
| CDIP1 | SCZ | sQTL | 0.03 | 1.13 | 2 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| DLG1 | SCZ | sQTL | 0.03 | 0.44 | 2 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| CTSH | AD | eQTL,sQTL | 0.39 | 1.18 | 1 | — | AD | no | yes | splicing-led (IsoGraph-resolved switch) |
| PGS1 | ALS | sQTL | 0.09 | 1.00 | 1 | — | ALS | no | yes | splicing-led (IsoGraph-resolved switch) |
| GGNBP2 | ALS,SCZ | eQTL,sQTL | 0.06 | 0.20 | 1 | — | ALS | no | yes | splicing-led (IsoGraph-resolved switch) |
| TBC1D15 | PD | sQTL | 0.02 | 0.50 | 1 | — | PD | no | yes | splicing-led (IsoGraph-resolved switch) |
| PRRC2B | SCZ | sQTL | 0.02 | 0.34 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| RTEL1 | AD,SCZ | sQTL | 0.02 | 0.64 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| TPCN1 | AD | sQTL | 0.01 | 0.56 | 1 | — | AD | no | yes | splicing-led (IsoGraph-resolved switch) |
| ARVCF | SCZ | sQTL | 0.01 | 1.00 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| EDEM3 | SCZ | eQTL | 0.94 | 0.67 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| TPP1 | ALS | sQTL | 0.48 | 0.75 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| GPM6A | SCZ | eQTL,sQTL | 0.21 | 0.36 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| CR1 | AD | eQTL | 0.17 | 0.77 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PAK6 | SCZ | eQTL | 0.17 | 0.71 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| FDFT1 | PD | eQTL | 0.07 | 1.70 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| GALNT2 | SCZ | eQTL | 0.07 | 0.52 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| TMED4 | SCZ | eQTL,sQTL | 0.06 | 1.29 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ITPRIP | SCZ | eQTL | 0.06 | 1.14 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| DOC2A | AD,SCZ | eQTL,sQTL | 0.05 | 0.73 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| RNASEH2C | SCZ | eQTL,sQTL | 0.05 | 1.69 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| SH3GL2 | PD | eQTL | 0.05 | 0.94 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| ADAM10 | SCZ | eQTL | 0.05 | 0.28 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| GABBR2 | SCZ | sQTL | 0.04 | 0.38 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| STK39 | PD | eQTL | 0.04 | 0.38 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SPG7 | SCZ | sQTL | 0.04 | 1.43 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| PPIL2 | SCZ | eQTL,sQTL | 0.03 | 0.79 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| FCER1G | AD | eQTL | 0.03 | 0.70 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| TYK2 | LBD | eQTL | 0.03 | 0.62 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| BLNK | AD | eQTL | 0.03 | 0.51 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| TXNDC15 | ALS | sQTL | 0.03 | 0.90 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| MARK2 | SCZ | eQTL | 0.03 | 0.20 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| NDRG4 | SCZ | sQTL | 0.03 | 0.76 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| POGZ | SCZ | eQTL | 0.03 | 0.17 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PCGF3 | AD,PD | sQTL | 0.03 | 0.45 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ARHGAP44 | SCZ | eQTL | 0.02 | 0.59 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PITPNC1 | SCZ | eQTL | 0.02 | 0.42 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PTPRN | ALS | sQTL | 0.02 | 0.73 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SOWAHA | ALS | eQTL | 0.02 | n/a | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| WBP2 | SCZ | eQTL | 0.02 | 0.77 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PKD1 | PD | sQTL | 0.02 | 0.36 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| SCFD1 | ALS | eQTL | 0.02 | 0.70 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PPP1R13B | SCZ | eQTL | 0.02 | 0.46 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| TRAF3 | SCZ | eQTL | 0.02 | 0.31 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| NT5C2 | SCZ | eQTL,sQTL | 0.02 | 1.13 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| DYRK1A | PD | eQTL | 0.02 | 0.17 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| MEF2C | ALS | eQTL | 0.02 | 0.21 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PBX1 | SCZ | eQTL | 0.02 | 0.13 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| EIF3E | SCZ | sQTL | 0.02 | 0.76 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| MYO15A | AD | eQTL | 0.02 | 0.89 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| MYOM2 | SCZ | sQTL | 0.02 | 1.46 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| FAM221A | SCZ | sQTL | 0.02 | 1.55 | 0 | — | — | yes | yes | splicing (sQTL not resolved to switch pair) |
| RBFA | SCZ | sQTL | 0.01 | 0.98 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| CHRNB2 | SCZ | eQTL | 0.01 | 1.12 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| BAIAP3 | ALS | sQTL | 0.01 | 1.38 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| RWDD2A | SCZ | eQTL | 0.01 | n/a | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| MED15 | SCZ | eQTL | 0.01 | 0.41 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| MYO18A | ALS,SCZ | sQTL | 0.01 | 0.40 | 0 | — | — | yes | yes | splicing (sQTL not resolved to switch pair) |
| DPYSL5 | SCZ | eQTL | 0.01 | 0.27 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| TTC19 | PD | sQTL | 0.01 | 1.05 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| LRRC73 | SCZ | eQTL | 0.01 | n/a | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| G2E3 | ALS | eQTL | 0.01 | 0.60 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| NUS1 | SCZ | eQTL | 0.01 | n/a | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| GRAMD1B | AD | sQTL | 0.01 | 0.68 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| GRM4 | SCZ | eQTL | 0.01 | 0.46 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| MRPS10 | AD | sQTL | 0.01 | 1.34 | 0 | — | — | yes | yes | splicing (sQTL not resolved to switch pair) |
