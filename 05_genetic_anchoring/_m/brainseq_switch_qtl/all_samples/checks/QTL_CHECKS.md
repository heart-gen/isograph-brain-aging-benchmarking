# BrainSEQ QTL checks — `all_samples`

Two checks required before any BrainSEQ QTL result is read. Neither modifies the mapping outputs. Thresholds were fixed in `brainseq_qtl_checks.py` before the checks were run.

## Sign pin of the recomputed switch coordinate

`S_g` is PC1 of a within-gene CLR composition, oriented by `stable_sign` on the largest loading; recomputed on the expanded QTL cohort it can land on the opposite sign from the discovery fit. Per gene, the recomputed and discovery `S_g` are correlated over the libraries both contain (at least 20). A negative r is a flip; |r| < 0.5 marks a gene whose PC1 is a different axis on the expanded cohort, for which no sign makes a directional comparison meaningful. Pinned slopes are in `qtl/cis_qtl_switch.sign_pinned.parquet`; p-values, q-values and hidden factors are unaffected by construction.

| region | shared libraries | genes compared | sign flipped | median abs r | axis-unstable | swQTL genes | swQTL flipped | swQTL axis-unstable |
|---|---|---|---|---|---|---|---|---|
| caudate | 223 | 11,636 | 5,262 (0.452) | 0.345 | 6,808 | 1,583 | 443 | 554 |
| dlpfc | 198 | 11,541 | 5,232 (0.453) | 0.305 | 6,935 | 1,292 | 338 | 439 |
| hippocampus | 213 | 11,098 | 5,080 (0.458) | 0.330 | 6,597 | 956 | 269 | 338 |

## Positive control: the `A_g` eQTL arm against GTEx v11 brain eGenes

Tissue-matched GTEx v11 eGenes (BH q < 0.05): caudate ↔ Caudate basal ganglia, DLPFC ↔ Frontal Cortex BA9, hippocampus ↔ Hippocampus. **Pre-specified pass rule**, per region: (i) Storey pi1 of the BrainSEQ `A_g` permutation p among GTEx eGenes >= 0.5, and (ii) direction concordance at GTEx's lead variant >= 0.9 among aligned pairs with BrainSEQ nominal p < 1e-05 (at least 20 pairs). A broken donor join destroys (i); an allele-coding error drives (ii) towards 0 or 0.5.

The counted allele is read from the data (BrainSEQ allele frequency against the panel's ALT frequency), not assumed. Palindromic SNPs are kept only when both panels agree on REF, since both code REF on the GRCh38 forward strand.

| region | GTEx eGenes tested | pi1 | q<0.05 | top-1000 recovery | pi1 all genes | lead pairs aligned | counted allele (af r) | concordance all | strong pairs | concordance strong | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| caudate | 9,665 | 0.877 | 0.742 | 0.970 | 0.723 | 6,527 | ALT (0.96) | 0.887 | 3,230 | 0.989 | **PASS** |
| dlpfc | 9,575 | 0.873 | 0.732 | 0.969 | 0.705 | 6,360 | ALT (0.95) | 0.890 | 3,140 | 0.991 | **PASS** |
| hippocampus | 6,768 | 0.833 | 0.664 | 0.945 | 0.625 | 4,504 | ALT (0.95) | 0.881 | 2,030 | 0.987 | **PASS** |

