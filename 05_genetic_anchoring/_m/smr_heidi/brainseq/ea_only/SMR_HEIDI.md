# SMR + HEIDI on the signal-level colocalization nominations

QTL source: `brainseq` (`ea_only` arm). Instrument threshold `--peqtl-smr 5e-08`. HEIDI SNPs at p < 0.0015654, 3-20 SNPs, 2000 kb cis window. Significance: Bonferroni at 0.05 over the instrumented probes of the row's OWN family (primary confirmatory / secondary event-localization), never pooled across the two. HEIDI rejects a single shared variant at p < 0.01. LD reference: 1000G EUR.

## How these numbers may be read

- `b_SMR` is a signed ratio estimate relating genetically predicted phenotype to disease under a single-causal-variant, no-pleiotropy model. It does **not** establish causal direction and cannot distinguish causality from horizontal pleiotropy.
- BrainSEQ is the discovery cohort: agreement here is same-tissue genetic anchoring, not replication. Each gene has one `A_g` (abundance) and one `S_g` (switch) probe, compared with BrainSEQ's own coloc posterior for the same region and axis.
- `S_g` is PC1 of the gene's within-gene composition. Its `b_SMR` is sign-pinned to the discovery switch axis, so its sign is a direction along that axis, not "more splicing"; a probe with `sign_pinned = False` keeps its mapped orientation, and `axis_unstable` marks a gene whose recomputed axis differs from discovery. An `S_g` result names no intron or event.
- Failing to reject HEIDI is **not** evidence of a shared variant; at 169-229 donors per region HEIDI is weaker still than at GTEx. It reads "not rejected".
- A HEIDI rejection does **not** overrule a strong signal-level colocalization, and SMR significance does not promote a locus coloc did not support. Disagreements stay disagreements.

## What was testable, by family

Two families, corrected apart. The **primary confirmatory family** holds one pre-designated probe per gene; a gene's other introns form a **secondary event-localization family** with its own Bonferroni correction, so the primary threshold is not inflated by introns and the introns do not escape correction when they are discussed. `F` is the instrument strength `(b_eQTL/se_eQTL)^2` of the top cis-QTL SNP, reported on every row and never used to exclude one. It is summarized by its 5th percentile rather than a count below the conventional F < 10: an instrument that clears p < 5e-8 has |z| ≳ 5.4 and so F ≳ 30 (≈ 24 at the relaxed 1e-6 arm), so a weak-instrument count is zero by construction and carries no information; the per-row `weak_instrument` flag stays in `smr_results.parquet` for any run at a looser threshold.

| analysis | modality | family | probes | instrumented | threshold | F median [min-max] | F p5 | `no_instrument` | `instrumented_tested_null` | `smr_signal_heidi_unavailable` | `smr_signal_heidi_rejects` | `smr_heidi_supported` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | A_g | primary | 69 | 18 | 0.003 | 90 [32-825] | 34 | 51 | 11 | 0 | 1 | 6 |
| aging__ad | S_g | primary | 66 | 4 | 0.013 | 67 [31-99] | 35 | 62 | 2 | 0 | 0 | 2 |
| aging__als | A_g | primary | 60 | 11 | 0.005 | 58 [34-108] | 36 | 49 | 1 | 0 | 2 | 8 |
| aging__als | S_g | primary | 47 | 7 | 0.007 | 49 [31-103] | 32 | 40 | 2 | 0 | 3 | 2 |
| aging__lbd | A_g | primary | 3 | 0 | — | — | — | 3 | 0 | 0 | 0 | 0 |
| aging__lbd | S_g | primary | 3 | 0 | — | — | — | 3 | 0 | 0 | 0 | 0 |
| aging__pd | A_g | primary | 30 | 14 | 0.004 | 45 [31-553] | 34 | 16 | 4 | 0 | 1 | 9 |
| aging__pd | S_g | primary | 27 | 1 | 0.050 | 31 [31-31] | 31 | 26 | 0 | 0 | 0 | 1 |
| aging__scz | A_g | primary | 224 | 67 | 0.000746 | 52 [30-586] | 31 | 157 | 41 | 0 | 4 | 22 |
| aging__scz | S_g | primary | 203 | 15 | 0.003 | 46 [35-92] | 35 | 188 | 5 | 0 | 0 | 10 |
| brainseq-sczd__scz | A_g | primary | 39 | 8 | 0.006 | 65 [31-145] | 31 | 31 | 3 | 0 | 0 | 5 |
| brainseq-sczd__scz | S_g | primary | 37 | 1 | 0.050 | 35 [35-35] | 35 | 36 | 0 | 0 | 0 | 1 |

**`no_instrument` is not a negative result** — the probe was never tested, because no cis-QTL reached the instrument threshold. It is the largest cell in every arm here and must never be read as evidence against a locus.

## Agreement with coloc, primary confirmatory family

| analysis | modality | probes | instrumented | threshold | `coloc_and_smr_heidi_not_rejected` | `coloc_and_smr_heidi_rejected` | `coloc_and_smr_heidi_untestable` | `coloc_smr_not_significant` | `coloc_no_instrument` | `smr_without_coloc` | `smr_coloc_not_scored` | `neither_tested_null` | `neither_no_instrument` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | A_g | 69 | 18 | 0.003 | 3 | 0 | 0 | 0 | 2 | 4 | 0 | 11 | 49 |
| aging__ad | S_g | 66 | 4 | 0.013 | 1 | 0 | 0 | 0 | 1 | 1 | 0 | 2 | 61 |
| aging__als | A_g | 60 | 11 | 0.005 | 5 | 0 | 0 | 0 | 5 | 5 | 0 | 1 | 44 |
| aging__als | S_g | 47 | 7 | 0.007 | 1 | 2 | 0 | 1 | 1 | 2 | 0 | 1 | 39 |
| aging__lbd | A_g | 3 | 0 | — | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 2 |
| aging__lbd | S_g | 3 | 0 | — | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| aging__pd | A_g | 30 | 14 | 0.004 | 5 | 1 | 0 | 0 | 0 | 4 | 0 | 4 | 16 |
| aging__pd | S_g | 27 | 1 | 0.050 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 25 |
| aging__scz | A_g | 224 | 67 | 0.000746 | 14 | 4 | 0 | 2 | 11 | 8 | 0 | 39 | 146 |
| aging__scz | S_g | 203 | 15 | 0.003 | 9 | 0 | 0 | 2 | 8 | 1 | 0 | 3 | 180 |
| brainseq-sczd__scz | A_g | 39 | 8 | 0.006 | 2 | 0 | 0 | 0 | 3 | 3 | 0 | 3 | 28 |
| brainseq-sczd__scz | S_g | 37 | 1 | 0.050 | 1 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 34 |

## SNP attrition into HEIDI

HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so what survives is worth stating. Per SMR run; the per-probe end of the trace is `nsnp_HEIDI` in the results table, itself capped at 20.

**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends past its locus loses the remainder. That is not a QC failure and says nothing about the QTL data.

| step | median across runs |
|---|---|
| BESD SNPs in the probe's cis window | 23,302 |
| ... found in the LD panel by rsID | 100.0% |
| ... with matching alleles | 100.0% |
| ... also carried by the GWAS | 5,656 |
| dropped by SMR beyond that (frequency check, duplicates) | 0 |

Genuine harmonization loss — SNPs present in all three and still dropped — peaks at **6** SNPs in any run. The identifier and allele steps are the QC ones, and both sit at 100%.

**Do not compare BESD-SNP retention across QTL sources.** Its denominator is imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far larger set against the same GWAS loci. The ratio measures panel density, not quality.

## GTEx nominations on the BrainSEQ switch axis

Per GTEx signal-level nomination, the `S_g` agreement class in each region (`*` = the region's tissue-matched GTEx tissue carries the GTEx sQTL call).

| gene | trait | caudate | dlpfc | hippocampus |
|---|---|---|---|---|
| ACTR1B | scz | — | — | `coloc_no_instrument`* |
| AKT1 | ad | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| ASB3 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| AZI2 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| BCKDK | ad | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| C9orf72 | als | `coloc_and_smr_heidi_rejected`* | `coloc_and_smr_heidi_rejected`* | `smr_without_coloc`* |
| CATSPER2 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| CCDC122 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| CCDC62 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CCS | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CD46 | scz | `neither_tested_null` | `neither_no_instrument` | `neither_no_instrument` |
| CDHR3 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CDIP1 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument`* |
| COG7 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| COPA | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CRELD2 | scz | `coloc_and_smr_heidi_not_rejected`* | `coloc_and_smr_heidi_not_rejected`* | `coloc_and_smr_heidi_not_rejected`* |
| CTSB | pd | `neither_no_instrument` | `smr_without_coloc` | `neither_no_instrument` |
| DDRGK1 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| DGKZ | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| DNAJA3 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| DOC2A | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| DOC2A | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| EFHB | scz | `neither_tested_null` | `neither_no_instrument`* | `neither_no_instrument` |
| FAM120AOS | scz | `neither_no_instrument` | `neither_tested_null` | `neither_no_instrument` |
| FAM184A | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| FANCI | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| FGFR1 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| FNBP1 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| FOXN2 | scz | `coloc_and_smr_heidi_not_rejected` | `coloc_smr_not_significant` | `coloc_smr_not_significant` |
| G2E3 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| GABBR2 | scz | `neither_no_instrument` | — | — |
| GALNT15 | scz | `neither_no_instrument` | — | — |
| GLYCTK | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| GPM6A | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| GPR135 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| HMOX2 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| IDH3B | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| IFNAR2 | ad | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument`* |
| IKBIP | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| INO80E | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| INO80E | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| INTS8 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| IRF3 | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| ITGB1BP1 | ad | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| KLC1 | scz | `neither_no_instrument` | `neither_no_instrument` | `smr_without_coloc` |
| L3HYPDH | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| LPCAT4 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| MAD1L1 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| MAP2K5 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| MAP7D1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| MED19 | scz | — | — | `neither_no_instrument`* |
| MRPS33 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NCOR1 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NDUFAF7 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NDUFS3 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NEK4 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| NME4 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NMRAL1 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument`* |
| NSMAF | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NT5C2 | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument` |
| NUCB2 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NUP50 | scz | — | — | `neither_no_instrument`* |
| PAK6 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PAM16 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PBRM1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PCBP3 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PICALM | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PILRB | ad | `neither_no_instrument` | `neither_tested_null` | `neither_tested_null` |
| POLG | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PPIL2 | scz | `coloc_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument` |
| PPIP5K1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PRDM2 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| PRMT7 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PSMD6 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| PTPRN | als | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| RAD51C | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| RBM6 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| RCBTB1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| RCSD1 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| REEP2 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| RERE | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| RPS6KL1 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| SCFD1 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| SERPINB1 | ad | — | `neither_no_instrument` | `neither_no_instrument` |
| SETD6 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| SIRPA | ad | `neither_no_instrument`* | — | — |
| SLC39A13 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| SNAP91 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| SNCA | lbd | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument`* |
| SPAG9 | ad | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument`* |
| SPI1 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| SYT5 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| TAOK2 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TEAD4 | scz | `neither_no_instrument` | `neither_no_instrument` | — |
| THAP3 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TMED4 | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| TMEM175 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TNFSF13 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TPCN1 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TPP1 | als | `neither_no_instrument` | — | — |
| TSPAN31 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| TTC19 | pd | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| TUBGCP4 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TXNDC15 | als | — | — | `neither_no_instrument` |
| VWA5B2 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| WIPI2 | als | `smr_without_coloc` | `neither_no_instrument` | `neither_no_instrument` |
| YPEL3 | ad | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument` |
| YPEL3 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| YWHAB | scz | `coloc_and_smr_heidi_not_rejected`* | — | — |
| ZFYVE21 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| ZNF232 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| ZSWIM7 | pd | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |

## Primary probes

| gene | trait | region | axis | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | pin | agreement |
|---|---|---|---|---|---|---|---|---|---|
| IFNAR2 | ad | dlpfc | A_g | 0.896 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| IFNAR2 | ad | hippocampus | A_g | 0.873 (abf) | 0.138 (0.042) | 0.001 | 0.296 (6) |  | `coloc_and_smr_heidi_not_rejected` |
| INO80E | ad | caudate | A_g | 0.733 (susie) | 0.180 (0.039) | 3.44e-06 | 0.384 (20) |  | `smr_without_coloc` |
| INO80E | ad | dlpfc | A_g | 0.604 (susie) | 0.098 (0.020) | 8.62e-07 | 0.341 (20) |  | `smr_without_coloc` |
| INO80E | ad | hippocampus | A_g | 0.666 (susie) | 0.149 (0.029) | 1.83e-07 | 0.334 (20) |  | `smr_without_coloc` |
| NDUFS3 | ad | caudate | A_g | 0.000192 (abf) | -0.197 (0.054) | 0.000225 | 1.22e-05 (20) |  | `smr_without_coloc` |
| PILRB | ad | caudate | A_g | 2.78e-14 (susie) | 0.003 (0.010) | 0.722 | 0.074 (20) |  | `neither_tested_null` |
| PILRB | ad | dlpfc | A_g | 2.85e-14 (susie) | 0.004 (0.009) | 0.689 | 0.323 (20) |  | `neither_tested_null` |
| PILRB | ad | dlpfc | S_g | 2e-12 (susie) | -0.003 (0.009) | 0.715 | 0.351 (20) | pinned | `neither_tested_null` |
| PILRB | ad | hippocampus | A_g | 3.22e-14 (susie) | 0.006 (0.009) | 0.524 | 0.290 (20) |  | `neither_tested_null` |
| PILRB | ad | hippocampus | S_g | 3.11e-09 (susie) | -0.004 (0.011) | 0.722 | 0.392 (20) | pinned | `neither_tested_null` |
| RAD51C | ad | dlpfc | A_g | 0.000284 (susie) | 0.018 (0.011) | 0.110 | 0.199 (20) |  | `neither_tested_null` |
| RAD51C | ad | hippocampus | A_g | 0.000281 (susie) | 0.025 (0.018) | 0.172 | 0.506 (20) |  | `neither_tested_null` |
| SCAMP4 | ad | caudate | S_g | 0.878 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| SERPINB1 | ad | caudate | A_g | 0.012 (abf) | -0.043 (0.024) | 0.080 | 0.205 (20) |  | `neither_tested_null` |
| SIRPA | ad | caudate | A_g | 0.960 (abf) | 0.104 (0.025) | 4.02e-05 | 0.878 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | dlpfc | A_g | 0.961 (abf) | 0.166 (0.044) | 0.000143 | 0.835 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | hippocampus | A_g | 0.946 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| TMEM106B | ad | caudate | S_g | 0.703 (susie) | -0.060 (0.020) | 0.003 | 0.046 (20) | pinned | `smr_without_coloc` |
| TMEM106B | ad | dlpfc | S_g | 0.846 (susie) | -0.045 (0.014) | 0.002 | 0.242 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| VWA5B2 | ad | caudate | A_g | 0.078 (abf) | 0.072 (0.026) | 0.005 | 0.028 (20) |  | `neither_tested_null` |
| VWA5B2 | ad | dlpfc | A_g | 0.078 (abf) | 0.056 (0.020) | 0.005 | 0.050 (20) |  | `neither_tested_null` |
| VWA5B2 | ad | hippocampus | A_g | 0.093 (abf) | 0.105 (0.042) | 0.013 | 0.113 (19) |  | `neither_tested_null` |
| ZNF232 | ad | caudate | A_g | 0.046 (abf) | -0.045 (0.018) | 0.014 | 0.275 (6) |  | `neither_tested_null` |
| ZNF232 | ad | hippocampus | A_g | 1.71e-08 (abf) | -0.021 (0.015) | 0.153 | 0.238 (13) |  | `neither_tested_null` |
| C9orf72 | als | caudate | A_g | 6.86e-05 (abf) | 0.183 (0.041) | 6.6e-06 | 0.000147 (20) |  | `smr_without_coloc` |
| C9orf72 | als | caudate | S_g | 0.893 (abf) | 0.137 (0.023) | 1.99e-09 | 0.002 (20) | pinned | `coloc_and_smr_heidi_rejected` |
| C9orf72 | als | dlpfc | A_g | 4.35e-13 (abf) | 0.140 (0.031) | 7.51e-06 | 9.38e-05 (20) |  | `smr_without_coloc` |
| C9orf72 | als | dlpfc | S_g | 0.983 (abf) | 0.172 (0.021) | 3.57e-16 | 0.002 (20) | pinned | `coloc_and_smr_heidi_rejected` |
| C9orf72 | als | hippocampus | S_g | 0.139 (abf) | 0.124 (0.020) | 9.32e-10 | 0.005 (20) | pinned | `smr_without_coloc` |
| G2E3 | als | dlpfc | A_g | 0.822 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| GGNBP2 | als | caudate | A_g | 0.947 (susie) | 0.252 (0.064) | 8.04e-05 | 0.559 (8) |  | `coloc_and_smr_heidi_not_rejected` |
| GGNBP2 | als | hippocampus | A_g | 0.829 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| MYO19 | als | caudate | A_g | 0.953 (susie) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| MYO19 | als | dlpfc | A_g | 0.570 (susie) | -0.179 (0.051) | 0.000415 | 0.041 (10) |  | `smr_without_coloc` |
| MYO19 | als | dlpfc | S_g | 0.930 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| MYO19 | als | hippocampus | A_g | 0.528 (susie) | -0.178 (0.049) | 0.00032 | 0.118 (8) |  | `smr_without_coloc` |
| PRDM2 | als | dlpfc | A_g | 0.003 (abf) | -0.006 (0.065) | 0.932 | 0.353 (17) |  | `neither_tested_null` |
| PTPRN | als | caudate | A_g | 0.820 (susie) | -0.156 (0.036) | 1.69e-05 | 0.247 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SARM1 | als | caudate | S_g | 0.393 (susie) | 0.072 (0.029) | 0.012 | 0.170 (20) | pinned | `neither_tested_null` |
| SARM1 | als | dlpfc | S_g | 0.960 (susie) | 0.086 (0.033) | 0.009 | 0.325 (14) | pinned | `coloc_smr_not_significant` |
| SARM1 | als | hippocampus | A_g | 0.979 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| SARM1 | als | hippocampus | S_g | 0.999 (susie) | 0.126 (0.029) | 2.05e-05 | 0.281 (15) | pinned | `coloc_and_smr_heidi_not_rejected` |
| SCFD1 | als | caudate | A_g | 0.946 (susie) | 0.349 (0.070) | 6.96e-07 | 0.296 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SCFD1 | als | dlpfc | A_g | 0.937 (susie) | 0.197 (0.032) | 6.89e-10 | 0.444 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SCFD1 | als | hippocampus | A_g | 0.931 (susie) | 0.422 (0.079) | 8.61e-08 | 0.089 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| TNFSF13 | als | caudate | A_g | 0.873 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| TXNDC15 | als | hippocampus | A_g | 0.317 (abf) | -0.177 (0.057) | 0.002 | 0.189 (20) |  | `smr_without_coloc` |
| WIPI2 | als | caudate | S_g | 0.088 (abf) | -0.090 (0.031) | 0.004 | 0.016 (20) | pinned | `smr_without_coloc` |
| SNCA | lbd | caudate | A_g | 0.857 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| CDHR3 | pd | caudate | A_g | 0.008 (abf) | -0.054 (0.084) | 0.519 | 0.667 (13) |  | `neither_tested_null` |
| CDHR3 | pd | dlpfc | A_g | 0.009 (abf) | -0.045 (0.070) | 0.520 | 0.735 (20) |  | `neither_tested_null` |
| CTSB | pd | caudate | A_g | 0.891 (susie) | -0.161 (0.043) | 0.000165 | 0.000167 (20) |  | `coloc_and_smr_heidi_rejected` |
| CTSB | pd | dlpfc | A_g | 0.870 (susie) | -0.190 (0.055) | 0.000499 | 0.017 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CTSB | pd | dlpfc | S_g | 0.218 (abf) | 0.104 (0.038) | 0.006 | 0.079 (20) | pinned | `smr_without_coloc` |
| CTSB | pd | hippocampus | A_g | 0.672 (susie) | -0.281 (0.089) | 0.002 | 0.312 (20) |  | `smr_without_coloc` |
| DDRGK1 | pd | caudate | A_g | 0.081 (abf) | -0.044 (0.072) | 0.543 | 0.020 (20) |  | `neither_tested_null` |
| DDRGK1 | pd | dlpfc | A_g | 0.985 (abf) | 0.315 (0.083) | 0.000148 | 0.045 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| NCOR1 | pd | caudate | A_g | 0.892 (susie) | 0.444 (0.104) | 2.14e-05 | 0.058 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | caudate | A_g | 0.001 (abf) | 0.169 (0.091) | 0.064 | 0.359 (20) |  | `neither_tested_null` |
| SH3GL2 | pd | dlpfc | A_g | 0.964 (abf) | 0.198 (0.065) | 0.002 | 0.171 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| STX4 | pd | caudate | S_g | 0.942 (susie) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| TTC19 | pd | caudate | A_g | 0.732 (susie) | -0.326 (0.081) | 5.9e-05 | 0.222 (18) |  | `smr_without_coloc` |
| ZSWIM7 | pd | caudate | A_g | 0.775 (susie) | 0.081 (0.016) | 5.77e-07 | 0.089 (20) |  | `smr_without_coloc` |
| ZSWIM7 | pd | dlpfc | A_g | 0.767 (susie) | 0.087 (0.018) | 8.62e-07 | 0.147 (20) |  | `smr_without_coloc` |
| ZSWIM7 | pd | hippocampus | A_g | 0.873 (susie) | 0.177 (0.038) | 4.35e-06 | 0.378 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | caudate | A_g | 4.57e-05 (abf) | 0.024 (0.037) | 0.511 | 0.388 (15) |  | `neither_tested_null` |
| ACTR1B | scz | caudate | A_g | 4.57e-05 (abf) | 0.024 (0.037) | 0.511 | 0.388 (15) |  | `neither_tested_null` |
| ACTR1B | scz | dlpfc | A_g | 0.902 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| ACTR1B | scz | dlpfc | A_g | 0.902 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| ACTR1B | scz | hippocampus | S_g | 0.858 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| ACTR1B | scz | hippocampus | S_g | 0.858 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| ARL14EP | scz | caudate | S_g | 0.834 (abf) | -0.101 (0.031) | 0.001 | 0.090 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| ARL14EP | scz | caudate | S_g | 0.834 (abf) | -0.101 (0.031) | 0.001 | 0.090 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| ARL14EP | scz | dlpfc | A_g | 0.794 (abf) | 0.093 (0.025) | 0.000183 | 0.198 (20) |  | `smr_without_coloc` |
| ARL14EP | scz | dlpfc | A_g | 0.794 (abf) | 0.093 (0.025) | 0.000183 | 0.198 (20) |  | `smr_without_coloc` |
| ATE1 | scz | caudate | S_g | 0.912 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| ATE1 | scz | caudate | S_g | 0.912 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| ATE1 | scz | hippocampus | A_g | 0.716 (abf) | -0.158 (0.049) | 0.001 | 0.038 (19) |  | `neither_tested_null` |
| ATE1 | scz | hippocampus | A_g | 0.716 (abf) | -0.158 (0.049) | 0.001 | 0.038 (19) |  | `smr_without_coloc` |
| CCDC122 | scz | caudate | A_g | 0.956 (abf) | -0.067 (0.015) | 8.12e-06 | 0.190 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | dlpfc | A_g | 0.161 (abf) | -0.028 (0.017) | 0.093 | 0.016 (20) |  | `neither_tested_null` |
| CCDC122 | scz | hippocampus | A_g | 0.953 (abf) | -0.078 (0.018) | 9.86e-06 | 0.022 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CD46 | scz | caudate | A_g | 0.846 (abf) | 0.142 (0.039) | 0.000302 | 0.039 (17) |  | `coloc_and_smr_heidi_not_rejected` |
| CD46 | scz | caudate | S_g | 0.095 (abf) | 0.035 (0.017) | 0.044 | 0.053 (14) | pinned | `neither_tested_null` |
| CD46 | scz | dlpfc | A_g | 0.842 (abf) | 0.079 (0.020) | 0.000127 | 0.506 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CD46 | scz | hippocampus | A_g | 0.852 (abf) | 0.118 (0.031) | 0.000136 | 0.370 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CNOT7 | scz | caudate | A_g | 0.146 (abf) | 0.095 (0.044) | 0.029 | 0.012 (20) |  | `neither_tested_null` |
| CNOT7 | scz | caudate | S_g | 0.970 (abf) | 0.110 (0.027) | 3.81e-05 | 0.102 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| CNOT7 | scz | dlpfc | A_g | 0.002 (abf) | 0.048 (0.021) | 0.022 | 0.009 (20) |  | `neither_tested_null` |
| CNOT7 | scz | dlpfc | S_g | 0.993 (abf) | 0.091 (0.021) | 1.49e-05 | 0.048 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| CNOT7 | scz | hippocampus | A_g | 0.948 (abf) | 0.163 (0.039) | 3.33e-05 | 0.232 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CNOT7 | scz | hippocampus | S_g | 0.973 (abf) | 0.096 (0.023) | 4.24e-05 | 0.173 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | caudate | A_g | 0.200 (abf) | -0.104 (0.039) | 0.007 | 0.052 (18) |  | `neither_tested_null` |
| COPA | scz | dlpfc | A_g | 0.199 (abf) | -0.076 (0.025) | 0.002 | 0.037 (20) |  | `neither_tested_null` |
| COPA | scz | hippocampus | A_g | 0.243 (abf) | -0.143 (0.058) | 0.014 | 0.062 (20) |  | `neither_tested_null` |
| CRELD2 | scz | caudate | S_g | 0.924 (susie) | -0.061 (0.013) | 2.79e-06 | 0.369 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | dlpfc | S_g | 0.917 (susie) | -0.066 (0.018) | 0.000302 | 0.344 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | hippocampus | A_g | 0.702 (abf) | 0.134 (0.038) | 0.000405 | 0.309 (20) |  | `smr_without_coloc` |
| CRELD2 | scz | hippocampus | S_g | 0.815 (susie) | -0.062 (0.016) | 7.04e-05 | 0.419 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| EFHB | scz | caudate | A_g | 0.001 (abf) | 0.000393 (0.012) | 0.973 | 0.001 (20) |  | `neither_tested_null` |
| EFHB | scz | caudate | S_g | 0.001 (abf) | -0.000479 (0.014) | 0.973 | 0.270 (11) | pinned | `neither_tested_null` |
| EFHB | scz | dlpfc | A_g | 0.001 (abf) | 0.000247 (0.007) | 0.973 | 0.009 (20) |  | `neither_tested_null` |
| EFHB | scz | hippocampus | A_g | 0.001 (abf) | 0.000361 (0.011) | 0.973 | 0.004 (20) |  | `neither_tested_null` |
| FAM120AOS | scz | caudate | A_g | 0.586 (abf) | 0.094 (0.030) | 0.002 | 0.000845 (20) |  | `neither_tested_null` |
| FAM120AOS | scz | dlpfc | A_g | 0.463 (abf) | -0.002 (0.018) | 0.926 | 0.001 (20) |  | `neither_tested_null` |
| FAM120AOS | scz | dlpfc | S_g | 0.123 (abf) | -0.010 (0.015) | 0.501 | 0.119 (18) | pinned, axis unstable | `neither_tested_null` |
| FAM120AOS | scz | hippocampus | A_g | 0.938 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| FAM184A | scz | dlpfc | A_g | 0.012 (abf) | 0.026 (0.013) | 0.042 | 0.025 (20) |  | `neither_tested_null` |
| FAM184A | scz | hippocampus | A_g | 0.017 (abf) | 0.060 (0.030) | 0.047 | 0.058 (20) |  | `neither_tested_null` |
| FOXN2 | scz | caudate | A_g | 0.968 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| FOXN2 | scz | caudate | S_g | 0.978 (abf) | 0.091 (0.027) | 0.000886 | 0.159 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| FOXN2 | scz | dlpfc | A_g | 0.920 (abf) | -0.239 (0.072) | 0.000946 | 0.229 (20) |  | `coloc_smr_not_significant` |
| FOXN2 | scz | dlpfc | S_g | 0.894 (abf) | 0.066 (0.027) | 0.016 | 0.109 (20) | pinned | `coloc_smr_not_significant` |
| FOXN2 | scz | hippocampus | A_g | 0.916 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| FOXN2 | scz | hippocampus | S_g | 0.887 (abf) | 0.048 (0.030) | 0.109 | 0.016 (20) | pinned | `coloc_smr_not_significant` |
| IKBIP | scz | caudate | A_g | 1.19e-05 (susie) | -0.002 (0.017) | 0.900 | 0.335 (20) |  | `neither_tested_null` |
| IKBIP | scz | dlpfc | A_g | 1.16e-05 (susie) | -0.002 (0.014) | 0.900 | 0.340 (20) |  | `neither_tested_null` |
| IKBIP | scz | hippocampus | A_g | 1.36e-05 (susie) | -0.005 (0.021) | 0.808 | 0.262 (20) |  | `neither_tested_null` |
| INO80E | scz | caudate | A_g | 0.970 (susie) | 0.235 (0.045) | 2.07e-07 | 0.185 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| INO80E | scz | dlpfc | A_g | 0.974 (susie) | 0.131 (0.023) | 1.27e-08 | 0.000172 (20) |  | `coloc_and_smr_heidi_rejected` |
| INO80E | scz | hippocampus | A_g | 0.970 (susie) | 0.195 (0.032) | 1.76e-09 | 0.352 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| KLC1 | scz | caudate | A_g | 4.46e-08 (abf) | -0.112 (0.036) | 0.002 | 6.32e-05 (20) |  | `neither_tested_null` |
| KLC1 | scz | dlpfc | A_g | 7.96e-05 (abf) | -0.086 (0.028) | 0.002 | 0.009 (20) |  | `neither_tested_null` |
| KLC1 | scz | hippocampus | A_g | 0.00022 (abf) | -0.141 (0.049) | 0.004 | 0.092 (16) |  | `neither_tested_null` |
| KLC1 | scz | hippocampus | S_g | 0.590 (abf) | 0.187 (0.043) | 1.27e-05 | 0.152 (20) | pinned | `smr_without_coloc` |
| L3HYPDH | scz | caudate | A_g | 0.001 (susie) | 0.063 (0.024) | 0.010 | 0.015 (20) |  | `neither_tested_null` |
| L3HYPDH | scz | dlpfc | A_g | 0.003 (susie) | 0.060 (0.024) | 0.012 | 0.078 (20) |  | `neither_tested_null` |
| LPCAT4 | scz | caudate | A_g | 9.21e-06 (abf) | -0.005 (0.039) | 0.899 | 0.272 (6) |  | `neither_tested_null` |
| LPCAT4 | scz | caudate | A_g | 0.002 (susie) | -0.005 (0.039) | 0.899 | 0.272 (6) |  | `neither_tested_null` |
| MAP2K5 | scz | caudate | A_g | 0.007 (abf) | 0.039 (0.022) | 0.074 | 0.004 (20) |  | `neither_tested_null` |
| MAP2K5 | scz | dlpfc | A_g | 0.007 (abf) | 0.039 (0.023) | 0.092 | 0.008 (20) |  | `neither_tested_null` |
| MYO19 | scz | caudate | A_g | 0.938 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| MYO19 | scz | dlpfc | A_g | 0.765 (abf) | -0.131 (0.039) | 0.000744 | 0.283 (10) |  | `smr_without_coloc` |
| MYO19 | scz | dlpfc | S_g | 0.880 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| MYO19 | scz | hippocampus | A_g | 0.726 (abf) | -0.126 (0.038) | 0.00076 | 0.357 (9) |  | `neither_tested_null` |
| NDUFAF7 | scz | caudate | A_g | 0.950 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| NEK4 | scz | caudate | A_g | 0.003 (abf) | 0.198 (0.042) | 2.36e-06 | 0.028 (20) |  | `smr_without_coloc` |
| NEK4 | scz | dlpfc | A_g | 0.368 (abf) | 0.142 (0.030) | 2.36e-06 | 0.030 (20) |  | `smr_without_coloc` |
| NEK4 | scz | hippocampus | A_g | 0.076 (abf) | 0.244 (0.053) | 4.58e-06 | 0.098 (20) |  | `smr_without_coloc` |
| NMRAL1 | scz | caudate | A_g | 0.000719 (susie) | -0.030 (0.028) | 0.292 | 0.651 (20) |  | `neither_tested_null` |
| NMRAL1 | scz | dlpfc | A_g | 0.000749 (susie) | -0.030 (0.028) | 0.294 | 0.828 (19) |  | `neither_tested_null` |
| NMRAL1 | scz | hippocampus | A_g | 0.687 (abf) | 0.134 (0.039) | 0.000644 | 0.254 (15) |  | `smr_without_coloc` |
| NT5C2 | scz | dlpfc | A_g | 0.843 (susie) | 0.027 (0.029) | 0.352 | 0.004 (20) |  | `coloc_smr_not_significant` |
| NUCB2 | scz | caudate | A_g | 0.000718 (abf) | -0.016 (0.034) | 0.635 | 0.028 (20) |  | `neither_tested_null` |
| NUP50 | scz | caudate | A_g | 0.687 (abf) | -0.138 (0.040) | 0.00052 | 0.441 (20) |  | `smr_without_coloc` |
| NUP50 | scz | dlpfc | A_g | 0.768 (abf) | -0.119 (0.038) | 0.001 | 0.446 (12) |  | `neither_tested_null` |
| PAK6 | scz | caudate | A_g | 0.935 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| PAK6 | scz | caudate | A_g | 0.935 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| PAK6 | scz | dlpfc | A_g | 0.0006 (abf) | -0.150 (0.055) | 0.006 | 0.222 (3) |  | `neither_tested_null` |
| PAK6 | scz | dlpfc | A_g | 0.0006 (abf) | -0.150 (0.055) | 0.006 | 0.222 (3) |  | `smr_without_coloc` |
| PLXNB2 | scz | dlpfc | S_g | 0.885 (susie) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| POLG | scz | dlpfc | A_g | 0.081 (abf) | -0.094 (0.029) | 0.001 | 0.109 (20) |  | `neither_tested_null` |
| PPDPF | scz | caudate | A_g | 0.000998 (susie) | 0.064 (0.046) | 0.157 | 2.15e-05 (7) |  | `neither_tested_null` |
| PPDPF | scz | caudate | S_g | 0.890 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| PPIL2 | scz | caudate | S_g | 0.858 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| PRMT7 | scz | caudate | A_g | 0.956 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| PRMT7 | scz | dlpfc | A_g | 0.985 (abf) | 0.152 (0.038) | 5.31e-05 | 0.640 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| PRMT7 | scz | hippocampus | A_g | 0.976 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| PSMD6 | scz | caudate | A_g | 0.650 (abf) | -0.094 (0.038) | 0.013 | 0.002 (20) |  | `neither_tested_null` |
| RBM6 | scz | caudate | A_g | 0.813 (abf) | 0.149 (0.038) | 8.85e-05 | 0.000811 (20) |  | `coloc_and_smr_heidi_rejected` |
| RBM6 | scz | dlpfc | A_g | 0.837 (abf) | 0.084 (0.020) | 3.76e-05 | 4.54e-05 (20) |  | `coloc_and_smr_heidi_rejected` |
| RBM6 | scz | hippocampus | A_g | 0.810 (abf) | 0.100 (0.024) | 4.79e-05 | 0.001 (20) |  | `coloc_and_smr_heidi_rejected` |
| RCBTB1 | scz | caudate | A_g | 0.977 (susie) | -0.103 (0.024) | 2.57e-05 | 0.312 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| RCBTB1 | scz | dlpfc | A_g | 0.975 (susie) | -0.100 (0.024) | 2.73e-05 | 0.630 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| RCBTB1 | scz | hippocampus | A_g | 0.973 (susie) | -0.154 (0.039) | 9.48e-05 | 0.316 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SETD6 | scz | dlpfc | A_g | 0.945 (susie) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| SNAP91 | scz | caudate | A_g | 0.963 (susie) | 0.307 (0.076) | 5.75e-05 | 0.119 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | dlpfc | A_g | 0.966 (susie) | 0.247 (0.063) | 7.71e-05 | 0.198 (15) |  | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | hippocampus | A_g | 0.959 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| SYT5 | scz | dlpfc | A_g | 0.004 (abf) | 0.009 (0.038) | 0.811 | 0.003 (9) |  | `neither_tested_null` |
| TAOK2 | scz | caudate | A_g | 0.933 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| THAP3 | scz | caudate | A_g | 0.006 (susie) | -0.038 (0.029) | 0.188 | 0.056 (5) |  | `neither_tested_null` |
| TMED4 | scz | caudate | A_g | 0.228 (abf) | -0.139 (0.048) | 0.004 | 0.100 (11) |  | `neither_tested_null` |
| TMED4 | scz | dlpfc | A_g | 0.178 (abf) | -0.108 (0.035) | 0.002 | 0.109 (20) |  | `neither_tested_null` |
| VPS29 | scz | dlpfc | S_g | 0.822 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| YPEL1 | scz | dlpfc | A_g | 0.601 (abf) | -0.078 (0.024) | 0.000892 | 0.925 (20) |  | `neither_tested_null` |
| YWHAB | scz | caudate | A_g | 0.845 (abf) | 0.231 (0.068) | 0.000633 | 0.074 (11) |  | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | caudate | S_g | 0.837 (abf) | 0.061 (0.016) | 0.000164 | 0.265 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | dlpfc | A_g | 0.856 (abf) | 0.130 (0.037) | 0.000473 | 0.136 (7) |  | `coloc_and_smr_heidi_not_rejected` |
| ZNF592 | scz | hippocampus | S_g | 0.888 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |

Probes with neither an SMR instrument nor a coloc call are in `smr_results.parquet` and not listed here.

