# BrainSEQ QTL checks — `ea_only`

Two checks required before any BrainSEQ QTL result is read. Neither modifies the mapping outputs. Thresholds were fixed in `brainseq_qtl_checks.py` before the checks were run.

## Sign pin of the recomputed switch coordinate

`S_g` is PC1 of a within-gene CLR composition, oriented by `stable_sign` on the largest loading; recomputed on the expanded QTL cohort it can land on the opposite sign from the discovery fit. Per gene, the recomputed and discovery `S_g` are correlated over the libraries both contain (at least 20). A negative r is a flip; |r| < 0.5 marks a gene whose PC1 is a different axis on the expanded cohort, for which no sign makes a directional comparison meaningful. Pinned slopes are in `qtl/cis_qtl_switch.sign_pinned.parquet`; p-values, q-values and hidden factors are unaffected by construction.

**Pre-specified gate:** the recomputed `S_g` must reproduce the discovery coordinate at median |r| >= 0.99. Below that the arm mapped a different phenotype, and no sign pin can repair it: re-map before reading any swQTL result. (The first all_samples arm omitted the discovery transcript filter and sat at 0.31-0.35; recomputed with it, the same libraries reproduce at 1.000.)

| region | shared libraries | genes compared | median abs r | reproduces discovery | sign flipped | axis-unstable | swQTL genes | swQTL flipped | swQTL axis-unstable |
|---|---|---|---|---|---|---|---|---|---|
| caudate | 116 | 12,993 | 1.000 | **yes** | 1,174 (0.090) | 414 | 1,155 | 96 | 35 |
| dlpfc | 85 | 12,787 | 1.000 | **yes** | 1,222 (0.096) | 468 | 815 | 76 | 24 |
| hippocampus | 104 | 12,877 | 1.000 | **yes** | 1,196 (0.093) | 441 | 567 | 48 | 27 |

## Positive control: the `A_g` eQTL arm against GTEx v11 brain eGenes

Tissue-matched GTEx v11 eGenes (BH q < 0.05): caudate ↔ Caudate basal ganglia, DLPFC ↔ Frontal Cortex BA9, hippocampus ↔ Hippocampus. **Pre-specified pass rule**, per region: (i) Storey pi1 of the BrainSEQ `A_g` permutation p among GTEx eGenes >= 0.5, and (ii) direction concordance at GTEx's lead variant >= 0.9 among aligned pairs with BrainSEQ nominal p < 1e-05 (at least 20 pairs). A broken donor join destroys (i); an allele-coding error drives (ii) towards 0 or 0.5.

The counted allele is read from the data (BrainSEQ allele frequency against the panel's ALT frequency), not assumed. Palindromic SNPs are kept only when both panels agree on REF, since both code REF on the GRCh38 forward strand.

| region | GTEx eGenes tested | pi1 | q<0.05 | top-1000 recovery | pi1 all genes | lead pairs aligned | counted allele (af r) | concordance all | strong pairs | concordance strong | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| caudate | 9,651 | 0.763 | 0.571 | 0.937 | 0.576 | 6,131 | ALT (0.99) | 0.896 | 2,615 | 0.990 | **PASS** |
| dlpfc | 9,563 | 0.718 | 0.495 | 0.917 | 0.531 | 5,981 | ALT (0.99) | 0.890 | 2,208 | 0.992 | **PASS** |
| hippocampus | 6,709 | 0.677 | 0.442 | 0.883 | 0.444 | 4,267 | ALT (0.99) | 0.881 | 1,409 | 0.987 | **PASS** |

