# Per-gene mechanistic deep-dive — panel summary

Every colocalized disease gene, one row each (ranked splicing-led first, multi-locus above single). `resolved events` = colocalizing sQTLs that map onto a concordant IsoGraph switch pair (the splicing-led, DTU-without-DGE class); `multi-locus` flags genes with >=2 such events (multiple significant colocalizations, incl. cross-trait). Main-figure framing is a PI decision taken from the anchored gene summary (`08_integration/_m/anchored_gene_summary/`), not from this table; as of 2026-09-20 Fig 4A provisionally draws PRDM2. The remainder are supporting vignettes. See `<GENE>.md` for each.

| gene | traits | kinds | max CLPP | LOEUF | resolved events | multi-locus | concordant traits | BrainSeq rep | GO-inv | verdict |
|------|--------|-------|----------|-------|-----------------|-------------|-------------------|--------------|--------|---------|
| SPG7 | SCZ | sQTL | 0.04 | 1.43 | 10 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| NADSYN1 | SCZ | sQTL | 0.02 | 0.93 | 8 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| TMEM175 | AD,ALS,PD | eQTL,sQTL | 0.48 | 1.29 | 6 | yes | AD,ALS,PD | no | yes | splicing-led (IsoGraph-resolved switch) |
| FLCN | AD,SCZ | sQTL | 0.01 | 0.49 | 6 | yes | AD,SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| CRELD2 | SCZ | sQTL | 0.06 | 1.09 | 5 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| GSTO2 | SCZ | sQTL | 0.02 | 1.24 | 4 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| INO80E | SCZ | eQTL,sQTL | 0.03 | 1.28 | 3 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| PIGQ | ALS | sQTL | 0.01 | 1.27 | 3 | yes | ALS | no | yes | splicing-led (IsoGraph-resolved switch) |
| TARBP1 | SCZ | sQTL | 0.01 | 0.95 | 3 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| CTC1 | ALS | eQTL,sQTL | 0.05 | 0.78 | 2 | yes | ALS | no | no | splicing-led (IsoGraph-resolved switch) |
| IFNAR2 | AD | sQTL | 0.03 | 0.76 | 2 | yes | AD | no | no | splicing-led (IsoGraph-resolved switch) |
| DLG1 | SCZ | sQTL | 0.03 | 0.44 | 2 | yes | SCZ | no | no | splicing-led (IsoGraph-resolved switch) |
| DNAJA3 | SCZ | sQTL | 0.03 | 0.73 | 2 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| PCGF3 | AD,PD | sQTL | 0.02 | 0.45 | 2 | yes | PD | no | no | splicing-led (IsoGraph-resolved switch) |
| RPAIN | SCZ | sQTL | 0.02 | 1.32 | 2 | yes | SCZ | no | no | splicing-led (IsoGraph-resolved switch) |
| B3GAT1 | SCZ | eQTL,sQTL | 0.01 | 0.75 | 2 | yes | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| TPP1 | ALS | sQTL | 0.48 | 0.75 | 1 | — | ALS | no | yes | splicing-led (IsoGraph-resolved switch) |
| PPP6R2 | SCZ | sQTL | 0.06 | 0.77 | 1 | — | SCZ | no | no | splicing-led (IsoGraph-resolved switch) |
| PRDM2 | ALS | sQTL | 0.06 | 0.18 | 1 | — | ALS | yes | no | splicing-led (IsoGraph-resolved switch) |
| DOC2A | AD,SCZ | eQTL,sQTL | 0.05 | 0.73 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| SNAP91 | SCZ | sQTL | 0.05 | 0.36 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| THAP3 | SCZ | sQTL | 0.02 | 1.04 | 1 | — | SCZ | no | no | splicing-led (IsoGraph-resolved switch) |
| PTPRN | ALS | sQTL | 0.02 | 0.73 | 1 | — | ALS | no | yes | splicing-led (IsoGraph-resolved switch) |
| TBC1D15 | PD | sQTL | 0.02 | 0.50 | 1 | — | PD | no | no | splicing-led (IsoGraph-resolved switch) |
| NT5C2 | SCZ | eQTL,sQTL | 0.02 | 1.13 | 1 | — | SCZ | no | no | splicing-led (IsoGraph-resolved switch) |
| VAMP2 | ALS | sQTL | 0.01 | 0.18 | 1 | — | ALS | no | no | splicing-led (IsoGraph-resolved switch) |
| RBFA | SCZ | sQTL | 0.01 | 0.98 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| BAIAP3 | ALS | sQTL | 0.01 | 1.38 | 1 | — | ALS | no | no | splicing-led (IsoGraph-resolved switch) |
| CTSB | PD | sQTL | 0.01 | 1.71 | 1 | — | PD | no | no | splicing-led (IsoGraph-resolved switch) |
| GPR135 | SCZ | sQTL | 0.01 | 1.94 | 1 | — | SCZ | no | yes | splicing-led (IsoGraph-resolved switch) |
| DGKQ | ALS,LBD,PD | eQTL | 1.00 | 1.16 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| EDEM3 | SCZ | eQTL | 0.94 | 0.67 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| ACTR1B | SCZ | eQTL,sQTL | 0.56 | 0.97 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| CA14 | SCZ | eQTL | 0.43 | 1.16 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| GPM6A | SCZ | eQTL,sQTL | 0.25 | 0.36 | 0 | — | — | yes | no | splicing (sQTL not resolved to switch pair) |
| CD38 | PD | eQTL | 0.22 | 1.29 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| GALNT6 | AD | eQTL | 0.18 | 1.02 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SLC4A8 | AD | eQTL | 0.18 | 0.46 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PRSS36 | AD | eQTL | 0.13 | 1.22 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| CASS4 | AD | eQTL | 0.13 | 0.52 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| GSTO1 | SCZ | eQTL,sQTL | 0.12 | 1.24 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ECE2 | PD | sQTL | 0.10 | 0.93 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| RUNDC3B | SCZ | eQTL | 0.10 | 0.66 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| IDH3B | SCZ | eQTL,sQTL | 0.10 | 1.04 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| PABPC1 | AD | eQTL | 0.09 | 0.12 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| RPS6KL1 | ALS | eQTL,sQTL | 0.09 | 0.98 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| RCBTB1 | SCZ | eQTL | 0.09 | 0.75 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| TMEM163 | PD | eQTL | 0.08 | 1.02 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| FOXN2 | SCZ | eQTL,sQTL | 0.08 | 0.73 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| ARHGAP28 | SCZ | eQTL | 0.07 | 0.77 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| BCKDK | AD | sQTL | 0.07 | 0.44 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| FDFT1 | PD | eQTL | 0.07 | 1.70 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| RAD51C | AD,SCZ | sQTL | 0.06 | 1.13 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| TMED4 | SCZ | eQTL,sQTL | 0.06 | 1.29 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ITGB1BP1 | AD | sQTL | 0.06 | 0.91 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| GGNBP2 | ALS | eQTL,sQTL | 0.06 | 0.20 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| GLB1L2 | SCZ | eQTL | 0.05 | 1.02 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| ANKRD36B | SCZ | sQTL | 0.05 | 0.57 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SH3GL2 | PD | eQTL | 0.05 | 0.94 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| ADAM10 | SCZ | eQTL | 0.04 | 0.28 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| GABBR2 | SCZ | sQTL | 0.04 | 0.38 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| STK39 | PD | eQTL | 0.04 | 0.38 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| FAM120AOS | SCZ | sQTL | 0.04 | 1.18 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| NSMAF | ALS | eQTL,sQTL | 0.04 | 0.77 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SLC27A3 | SCZ | eQTL | 0.04 | 1.26 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| ZFYVE21 | SCZ | sQTL | 0.04 | 0.96 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| MYO19 | ALS | eQTL | 0.04 | 1.23 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| RAI1 | SCZ | sQTL | 0.04 | 0.11 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SNCA | LBD | sQTL | 0.04 | 0.40 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| RHPN1 | SCZ | sQTL | 0.04 | 1.60 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| PRMT7 | SCZ | sQTL | 0.03 | 0.94 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| TYK2 | LBD | eQTL | 0.03 | 0.62 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| TXNDC15 | ALS | sQTL | 0.03 | 0.90 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| LY6H | SCZ | eQTL,sQTL | 0.03 | 1.08 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| RANBP10 | ALS | eQTL | 0.03 | 0.48 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| RBM26 | SCZ | eQTL | 0.03 | 0.38 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| ADAM15 | LBD | sQTL | 0.03 | 0.75 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| CDIP1 | SCZ | sQTL | 0.03 | 1.13 | 0 | — | — | yes | yes | splicing (sQTL not resolved to switch pair) |
| ARL14EP | SCZ | eQTL | 0.03 | 0.82 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PLEKHA1 | AD | eQTL | 0.03 | 0.92 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| NDRG4 | SCZ | sQTL | 0.03 | 0.76 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| CCDC122 | SCZ | eQTL,sQTL | 0.03 | 1.53 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| POC1B | SCZ | eQTL | 0.03 | n/a | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| ITSN1 | SCZ | eQTL | 0.03 | 0.30 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| NMRAL1 | SCZ | sQTL | 0.03 | n/a | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| INTS8 | AD | sQTL | 0.03 | 0.45 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ARHGAP44 | SCZ | eQTL | 0.03 | 0.59 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| KCNN3 | SCZ | eQTL | 0.03 | 0.38 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| MAP2K5 | SCZ | sQTL | 0.02 | 1.20 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| FBXL15 | SCZ | eQTL | 0.02 | 1.05 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| NPEPL1 | SCZ | sQTL | 0.02 | 1.23 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| RNH1 | SCZ | sQTL | 0.02 | 1.12 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| AXIN1 | ALS | eQTL | 0.02 | 0.44 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| RERE | SCZ | eQTL | 0.02 | 0.35 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| NUP85 | ALS | sQTL | 0.02 | 0.33 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| TNFSF13 | ALS | eQTL,sQTL | 0.02 | n/a | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SOWAHA | ALS | eQTL | 0.02 | n/a | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| ZSWIM7 | PD | sQTL | 0.02 | 1.60 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ZNF365 | SCZ | sQTL | 0.02 | 0.85 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| PTPRU | SCZ | eQTL | 0.02 | 0.52 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| LMF1 | ALS | sQTL | 0.02 | 1.16 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| SCFD1 | ALS | eQTL | 0.02 | 0.70 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| PLXNB2 | SCZ | sQTL | 0.02 | 0.55 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| WBP2 | SCZ | eQTL | 0.02 | 0.77 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| NAP1L1 | SCZ | sQTL | 0.02 | 0.22 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| DDX56 | SCZ | eQTL | 0.02 | 0.89 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| L3HYPDH | SCZ | sQTL | 0.02 | 1.51 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| FNBP1 | ALS | sQTL | 0.02 | 0.46 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| CNNM2 | SCZ | eQTL | 0.02 | 0.29 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| PPP1R13B | SCZ | eQTL | 0.02 | 0.46 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| TMEM219 | SCZ | eQTL | 0.02 | 1.43 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| DUS2 | ALS | eQTL | 0.02 | 0.86 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| TRAF3 | SCZ | eQTL | 0.02 | 0.31 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SF3B1 | SCZ | eQTL | 0.02 | 0.10 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| NMRK1 | SCZ | sQTL | 0.02 | 1.20 | 0 | — | — | yes | yes | splicing (sQTL not resolved to switch pair) |
| NME4 | ALS | sQTL | 0.02 | 1.52 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| MEF2C | ALS | eQTL | 0.02 | 0.21 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| MYO15A | AD | eQTL | 0.02 | 0.89 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| ANKRD54 | AD | eQTL | 0.02 | 0.89 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SMAP1 | PD | sQTL | 0.02 | 0.45 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| ALG12 | SCZ | eQTL | 0.01 | 0.96 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| SETD6 | SCZ | eQTL | 0.01 | 1.62 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| CACNA1E | AD | eQTL | 0.01 | 0.15 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| DECR2 | ALS | sQTL | 0.01 | n/a | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| OLA1 | SCZ | eQTL | 0.01 | 0.43 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| MVD | SCZ | eQTL | 0.01 | 1.60 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SLCO1A2 | ALS | sQTL | 0.01 | 1.21 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| TMEM107 | ALS | sQTL | 0.01 | 1.47 | 0 | — | — | yes | no | splicing (sQTL not resolved to switch pair) |
| TMEM63A | SCZ | sQTL | 0.01 | 0.84 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| FANCL | SCZ | sQTL | 0.01 | 1.35 | 0 | — | — | yes | yes | splicing (sQTL not resolved to switch pair) |
| MYOM2 | SCZ | sQTL | 0.01 | 1.46 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SEC61A2 | SCZ | sQTL | 0.01 | 0.66 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| LRRC63 | ALS | sQTL | 0.01 | 1.30 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| RTEL1 | AD | sQTL | 0.01 | 0.64 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| CAMK1 | SCZ | sQTL | 0.01 | 0.97 | 0 | — | — | yes | no | splicing (sQTL not resolved to switch pair) |
| TNNT2 | SCZ | sQTL | 0.01 | n/a | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| STX16 | SCZ | sQTL | 0.01 | 0.96 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| RWDD2A | SCZ | eQTL | 0.01 | n/a | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| NSUN2 | SCZ | sQTL | 0.01 | 0.96 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| MYO18A | ALS | sQTL | 0.01 | 0.40 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| THOP1 | ALS | sQTL | 0.01 | 0.84 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| RAB40C | ALS | sQTL | 0.01 | 0.47 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| ACHE | AD | sQTL | 0.01 | 0.49 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| SLC6A7 | AD | sQTL | 0.01 | 0.98 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| MRPS33 | SCZ | sQTL | 0.01 | 1.55 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| COG7 | AD | eQTL | 0.01 | 0.65 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| IFT172 | SCZ | sQTL | 0.01 | 0.69 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| TTC19 | PD | sQTL | 0.01 | 1.05 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| G2E3 | ALS | eQTL | 0.01 | 0.60 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| COL9A3 | AD | sQTL | 0.01 | 1.03 | 0 | — | — | no | yes | splicing (sQTL not resolved to switch pair) |
| CHRNB2 | SCZ | eQTL | 0.01 | 1.12 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| EBPL | AD | eQTL | 0.01 | 1.41 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| CNIH3 | ALS | sQTL | 0.01 | 0.98 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| TMEM8B | ALS | sQTL | 0.01 | 0.82 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| CPNE7 | SCZ | eQTL | 0.01 | 1.52 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| PCCB | SCZ | eQTL | 0.01 | 0.84 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SPIRE2 | SCZ | eQTL | 0.01 | 1.35 | 0 | — | — | no | yes | expression-led (eQTL gene-level) |
| SLC36A1 | LBD | sQTL | 0.01 | 0.74 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
| GRM4 | SCZ | eQTL | 0.01 | 0.46 | 0 | — | — | no | no | expression-led (eQTL gene-level) |
| SNRNP35 | PD | sQTL | 0.01 | 0.97 | 0 | — | — | no | no | splicing (sQTL not resolved to switch pair) |
