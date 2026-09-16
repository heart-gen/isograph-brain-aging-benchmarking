# SMR + HEIDI on the signal-level colocalization nominations

QTL source: `brainseq` (`ea_only` arm). Instrument threshold `--peqtl-smr 5e-08`. Multi-SNP SMR (`--smr-multi`, LD pruned at r2 0.1) -- **SENSITIVITY ARM**. HEIDI SNPs at p < 0.0015654, 3-20 SNPs, 2000 kb cis window. Significance: Bonferroni at 0.05 over the instrumented probes of the row's OWN family (primary confirmatory / secondary event-localization), never pooled across the two. HEIDI rejects a single shared variant at p < 0.01. LD reference: 1000G EUR.

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
| aging__ad | A_g | primary | 12 | 5 | 0.010 | 48 [41-91] | 41 | 7 | 2 | 0 | 1 | 2 |
| aging__ad | S_g | primary | 10 | 0 | — | — | — | 10 | 0 | 0 | 0 | 0 |
| aging__als | A_g | primary | 18 | 4 | 0.013 | 48 [34-103] | 35 | 14 | 1 | 0 | 0 | 3 |
| aging__als | S_g | primary | 10 | 1 | 0.050 | 35 [35-35] | 35 | 9 | 0 | 0 | 0 | 1 |
| aging__pd | A_g | primary | 18 | 9 | 0.006 | 43 [31-117] | 33 | 9 | 2 | 0 | 1 | 6 |
| aging__pd | S_g | primary | 15 | 1 | 0.050 | 31 [31-31] | 31 | 14 | 0 | 0 | 0 | 1 |
| aging__scz | A_g | primary | 75 | 21 | 0.002 | 41 [30-77] | 31 | 54 | 7 | 0 | 2 | 12 |
| aging__scz | S_g | primary | 59 | 2 | 0.025 | 51 [36-67] | 38 | 57 | 0 | 0 | 0 | 2 |
| brainseq-sczd__scz | A_g | primary | 12 | 2 | 0.025 | 80 [52-108] | 55 | 10 | 1 | 0 | 0 | 1 |
| brainseq-sczd__scz | S_g | primary | 10 | 1 | 0.050 | 35 [35-35] | 35 | 9 | 0 | 0 | 0 | 1 |

**`no_instrument` is not a negative result** — the probe was never tested, because no cis-QTL reached the instrument threshold. It is the largest cell in every arm here and must never be read as evidence against a locus.

## Agreement with coloc, primary confirmatory family

| analysis | modality | probes | instrumented | threshold | `coloc_and_smr_heidi_not_rejected` | `coloc_and_smr_heidi_rejected` | `coloc_and_smr_heidi_untestable` | `coloc_smr_not_significant` | `coloc_no_instrument` | `smr_without_coloc` | `smr_coloc_not_scored` | `neither_tested_null` | `neither_no_instrument` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | A_g | 12 | 5 | 0.010 | 2 | 0 | 0 | 0 | 1 | 1 | 0 | 2 | 6 |
| aging__ad | S_g | 10 | 0 | — | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 10 |
| aging__als | A_g | 18 | 4 | 0.013 | 2 | 0 | 0 | 0 | 1 | 1 | 0 | 1 | 13 |
| aging__als | S_g | 10 | 1 | 0.050 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 9 |
| aging__pd | A_g | 18 | 9 | 0.006 | 4 | 1 | 0 | 0 | 0 | 2 | 0 | 2 | 9 |
| aging__pd | S_g | 15 | 1 | 0.050 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 13 |
| aging__scz | A_g | 75 | 21 | 0.002 | 4 | 0 | 0 | 0 | 3 | 10 | 0 | 7 | 51 |
| aging__scz | S_g | 59 | 2 | 0.025 | 1 | 0 | 0 | 0 | 3 | 1 | 0 | 0 | 54 |
| brainseq-sczd__scz | A_g | 12 | 2 | 0.025 | 0 | 0 | 0 | 0 | 2 | 1 | 0 | 1 | 8 |
| brainseq-sczd__scz | S_g | 10 | 1 | 0.050 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 8 |

## SNP attrition into HEIDI

HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so what survives is worth stating. Per SMR run; the per-probe end of the trace is `nsnp_HEIDI` in the results table, itself capped at 20.

**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends past its locus loses the remainder. That is not a QC failure and says nothing about the QTL data.

| step | median across runs |
|---|---|
| BESD SNPs in the probe's cis window | 9,542 |
| ... found in the LD panel by rsID | 100.0% |
| ... with matching alleles | 100.0% |
| ... also carried by the GWAS | 4,223 |
| dropped by SMR beyond that (frequency check, duplicates) | 0 |

Genuine harmonization loss — SNPs present in all three and still dropped — peaks at **6** SNPs in any run. The identifier and allele steps are the QC ones, and both sit at 100%.

**Do not compare BESD-SNP retention across QTL sources.** Its denominator is imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far larger set against the same GWAS loci. The ratio measures panel density, not quality.

## GTEx nominations on the BrainSEQ switch axis

Per GTEx signal-level nomination, the `S_g` agreement class in each region (`*` = the region's tissue-matched GTEx tissue carries the GTEx sQTL call).

| gene | trait | caudate | dlpfc | hippocampus |
|---|---|---|---|---|
| ACTR1B | scz | — | — | `coloc_no_instrument`* |
| CDIP1 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument`* |
| COPA | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CTSB | pd | `neither_no_instrument` | `smr_without_coloc` | `neither_no_instrument` |
| DDRGK1 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| DGKZ | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| FGFR1 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| GABBR2 | scz | `neither_no_instrument` | — | — |
| GLYCTK | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| IRF3 | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| KLC1 | scz | `neither_no_instrument` | `neither_no_instrument` | `smr_without_coloc` |
| LPCAT4 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| MED19 | scz | — | — | `neither_no_instrument`* |
| MRPS33 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NCOR1 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NDUFAF7 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NDUFS3 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NEK4 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| NUP50 | scz | — | — | `neither_no_instrument`* |
| POLG | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PPIL2 | scz | `coloc_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument` |
| PRDM2 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| PSMD6 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| PTPRN | als | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| RASA1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| RERE | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| SIRPA | ad | `neither_no_instrument`* | — | — |
| SNAP91 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| SYT5 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| TAOK2 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TPCN1 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TTC19 | pd | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| TXNDC15 | als | — | — | `neither_no_instrument` |
| WIPI2 | als | `smr_without_coloc` | `neither_no_instrument` | `neither_no_instrument` |
| YWHAB | scz | `coloc_and_smr_heidi_not_rejected`* | — | — |
| ZNF232 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |

## Primary probes

| gene | trait | region | axis | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | pin | agreement |
|---|---|---|---|---|---|---|---|---|---|
| NDUFS3 | ad | caudate | A_g | 7.03e-05 (susie) | -0.197 (0.054) | 0.000225 | 1.22e-05 (20) |  | `smr_without_coloc` |
| SIRPA | ad | caudate | A_g | 0.960 (abf) | 0.104 (0.025) | 4.02e-05 | 0.878 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | dlpfc | A_g | 0.961 (abf) | 0.166 (0.044) | 0.000143 | 0.835 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | hippocampus | A_g | 0.946 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| ZNF232 | ad | caudate | A_g | 0.046 (abf) | -0.045 (0.018) | 0.014 | 0.275 (6) |  | `neither_tested_null` |
| ZNF232 | ad | hippocampus | A_g | 1.71e-08 (abf) | -0.021 (0.015) | 0.153 | 0.238 (13) |  | `neither_tested_null` |
| GGNBP2 | als | caudate | A_g | 0.947 (susie) | 0.252 (0.064) | 8.04e-05 | 0.559 (8) |  | `coloc_and_smr_heidi_not_rejected` |
| GGNBP2 | als | hippocampus | A_g | 0.829 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| PRDM2 | als | dlpfc | A_g | 0.003 (abf) | -0.006 (0.065) | 0.932 | 0.353 (17) |  | `neither_tested_null` |
| PTPRN | als | caudate | A_g | 0.834 (susie) | -0.156 (0.036) | 1.69e-05 | 0.247 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| TXNDC15 | als | hippocampus | A_g | 0.317 (abf) | -0.177 (0.057) | 0.002 | 0.189 (20) |  | `smr_without_coloc` |
| WIPI2 | als | caudate | S_g | 0.088 (abf) | -0.090 (0.031) | 0.004 | 0.016 (20) | pinned | `smr_without_coloc` |
| CTSB | pd | caudate | A_g | 0.891 (susie) | -0.161 (0.043) | 0.000165 | 0.000167 (20) |  | `coloc_and_smr_heidi_rejected` |
| CTSB | pd | dlpfc | A_g | 0.870 (susie) | -0.190 (0.055) | 0.000499 | 0.017 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| CTSB | pd | dlpfc | S_g | 0.218 (abf) | 0.104 (0.038) | 0.006 | 0.079 (20) | pinned | `smr_without_coloc` |
| CTSB | pd | hippocampus | A_g | 0.672 (susie) | -0.281 (0.089) | 0.002 | 0.312 (20) |  | `smr_without_coloc` |
| DDRGK1 | pd | caudate | A_g | 0.081 (abf) | -0.044 (0.072) | 0.543 | 0.020 (20) |  | `neither_tested_null` |
| DDRGK1 | pd | dlpfc | A_g | 0.985 (abf) | 0.315 (0.083) | 0.000148 | 0.045 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| NCOR1 | pd | caudate | A_g | 0.891 (susie) | 0.444 (0.104) | 2.14e-05 | 0.058 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | caudate | A_g | 0.001 (abf) | 0.169 (0.091) | 0.064 | 0.359 (20) |  | `neither_tested_null` |
| SH3GL2 | pd | dlpfc | A_g | 0.964 (abf) | 0.198 (0.065) | 0.002 | 0.171 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| STX4 | pd | caudate | S_g | 0.942 (susie) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| TTC19 | pd | caudate | A_g | 0.737 (susie) | -0.326 (0.081) | 5.9e-05 | 0.222 (18) |  | `smr_without_coloc` |
| ACTR1B | scz | caudate | A_g | 4.57e-05 (abf) | 0.024 (0.037) | 0.511 | 0.388 (15) |  | `neither_tested_null` |
| ACTR1B | scz | caudate | A_g | 4.57e-05 (abf) | 0.024 (0.037) | 0.511 | 0.388 (15) |  | `neither_tested_null` |
| ACTR1B | scz | dlpfc | A_g | 0.902 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| ACTR1B | scz | dlpfc | A_g | 0.902 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| ACTR1B | scz | hippocampus | S_g | 0.858 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| ACTR1B | scz | hippocampus | S_g | 0.858 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| ARL14EP | scz | caudate | S_g | 0.834 (abf) | -0.101 (0.031) | 0.001 | 0.090 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| ARL14EP | scz | dlpfc | A_g | 0.794 (abf) | 0.093 (0.025) | 0.000183 | 0.198 (20) |  | `smr_without_coloc` |
| COPA | scz | caudate | A_g | 0.200 (abf) | -0.104 (0.039) | 0.007 | 0.052 (18) |  | `neither_tested_null` |
| COPA | scz | dlpfc | A_g | 0.199 (abf) | -0.076 (0.025) | 0.002 | 0.037 (20) |  | `smr_without_coloc` |
| COPA | scz | hippocampus | A_g | 0.243 (abf) | -0.143 (0.058) | 0.014 | 0.062 (20) |  | `neither_tested_null` |
| KLC1 | scz | caudate | A_g | 8.14e-08 (susie) | -0.112 (0.036) | 0.002 | 6.32e-05 (20) |  | `smr_without_coloc` |
| KLC1 | scz | dlpfc | A_g | 0.003 (susie) | -0.086 (0.028) | 0.002 | 0.009 (20) |  | `smr_without_coloc` |
| KLC1 | scz | hippocampus | A_g | 0.001 (susie) | -0.141 (0.049) | 0.004 | 0.092 (16) |  | `neither_tested_null` |
| KLC1 | scz | hippocampus | S_g | 0.721 (susie) | 0.187 (0.043) | 1.27e-05 | 0.152 (20) | pinned | `smr_without_coloc` |
| LPCAT4 | scz | caudate | A_g | 0.002 (susie) | -0.005 (0.039) | 0.899 | 0.272 (6) |  | `neither_tested_null` |
| NDUFAF7 | scz | caudate | A_g | 0.948 (susie) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| NEK4 | scz | caudate | A_g | 0.003 (abf) | 0.198 (0.042) | 2.36e-06 | 0.028 (20) |  | `smr_without_coloc` |
| NEK4 | scz | dlpfc | A_g | 0.368 (abf) | 0.142 (0.030) | 2.36e-06 | 0.030 (20) |  | `smr_without_coloc` |
| NEK4 | scz | hippocampus | A_g | 0.076 (abf) | 0.244 (0.053) | 4.58e-06 | 0.098 (20) |  | `smr_without_coloc` |
| NUP50 | scz | caudate | A_g | 0.687 (abf) | -0.138 (0.040) | 0.00052 | 0.441 (20) |  | `smr_without_coloc` |
| NUP50 | scz | dlpfc | A_g | 0.768 (abf) | -0.119 (0.038) | 0.001 | 0.446 (12) |  | `smr_without_coloc` |
| POLG | scz | dlpfc | A_g | 0.081 (abf) | -0.094 (0.029) | 0.001 | 0.109 (20) |  | `smr_without_coloc` |
| PPIL2 | scz | caudate | S_g | 0.858 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| PSMD6 | scz | caudate | A_g | 0.611 (susie) | -0.094 (0.038) | 0.013 | 0.002 (20) |  | `neither_tested_null` |
| SNAP91 | scz | caudate | A_g | 0.963 (susie) | 0.307 (0.076) | 5.75e-05 | 0.119 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | dlpfc | A_g | 0.966 (susie) | 0.247 (0.063) | 7.71e-05 | 0.198 (15) |  | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | hippocampus | A_g | 0.959 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| SYT5 | scz | dlpfc | A_g | 0.004 (abf) | 0.009 (0.038) | 0.811 | 0.003 (9) |  | `neither_tested_null` |
| TAOK2 | scz | caudate | A_g | 0.933 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| YPEL1 | scz | dlpfc | A_g | 0.601 (abf) | -0.078 (0.024) | 0.000892 | 0.925 (20) |  | `smr_without_coloc` |
| YWHAB | scz | caudate | A_g | 0.845 (abf) | 0.231 (0.068) | 0.000633 | 0.074 (11) |  | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | caudate | S_g | 0.837 (abf) | 0.061 (0.016) | 0.000164 | 0.265 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | dlpfc | A_g | 0.856 (abf) | 0.130 (0.037) | 0.000473 | 0.136 (7) |  | `coloc_and_smr_heidi_not_rejected` |
| ZNF592 | scz | hippocampus | S_g | 0.888 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |

## The multi-SNP arm

`--smr-multi` combines the cis SNPs surviving LD pruning at r2 0.1 instead of testing the top SNP alone. It is a robustness arm for the SMR estimate, not a second discovery pass: it answers whether a signal rests on one lead SNP or on the cis signal more broadly.

`p_SMR` and `smr_status` in this table are still the single-SNP quantities, so rows match the primary arm one for one. The multi-SNP verdict is `p_SMR_multi` / `smr_multi_status`, Bonferroni-corrected over the probes the multi-SNP test actually ran on (`n_multi_family`) -- a smaller denominator than the instrumented count, because SMR skips the test where too few cis SNPs survive pruning. Those probes are `multi_unavailable`, which is **not** a null result.

| smr_multi_status | probes |
|---|---|
| `no_instrument` | 193 |
| `multi_tested_null` | 14 |
| `multi_signal_heidi_rejects` | 3 |
| `multi_heidi_supported` | 29 |

Of the 33 probes significant on the single-SNP test, 32 are also significant under the multi-SNP test and 0 could not be tested. A probe that does not survive is a signal carried by its lead SNP alone; that is a caveat on the SMR estimate, not a refutation of the colocalization.

Probes with neither an SMR instrument nor a coloc call are in `smr_results.parquet` and not listed here.

