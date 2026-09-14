# SMR + HEIDI on the signal-level colocalization nominations

QTL source: `brainseq` (`ea_only` arm). Instrument threshold `--peqtl-smr 1e-06` -- **SENSITIVITY ARM**. HEIDI SNPs at p < 0.0015654, 3-20 SNPs, 2000 kb cis window. Significance: Bonferroni at 0.05 over the instrumented probes of the row's OWN family (primary confirmatory / secondary event-localization), never pooled across the two. HEIDI rejects a single shared variant at p < 0.01. LD reference: 1000G EUR.

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
| aging__ad | A_g | primary | 30 | 8 | 0.006 | 45 [24-91] | 24 | 22 | 4 | 0 | 1 | 3 |
| aging__ad | S_g | primary | 26 | 2 | 0.025 | 41 [33-49] | 34 | 24 | 1 | 0 | 0 | 1 |
| aging__als | A_g | primary | 27 | 10 | 0.005 | 58 [26-108] | 26 | 17 | 4 | 0 | 0 | 6 |
| aging__als | S_g | primary | 27 | 3 | 0.017 | 27 [27-33] | 27 | 24 | 1 | 0 | 0 | 2 |
| aging__lbd | A_g | primary | 3 | 0 | — | — | — | 3 | 0 | 0 | 0 | 0 |
| aging__lbd | S_g | primary | 3 | 0 | — | — | — | 3 | 0 | 0 | 0 | 0 |
| aging__pd | A_g | primary | 18 | 5 | 0.010 | 43 [29-47] | 31 | 13 | 2 | 0 | 0 | 3 |
| aging__pd | S_g | primary | 15 | 0 | — | — | — | 15 | 0 | 0 | 0 | 0 |
| aging__scz | A_g | primary | 60 | 20 | 0.003 | 38 [24-596] | 25 | 40 | 7 | 0 | 5 | 8 |
| aging__scz | S_g | primary | 59 | 4 | 0.013 | 39 [28-125] | 29 | 55 | 3 | 0 | 0 | 1 |

**`no_instrument` is not a negative result** — the probe was never tested, because no cis-QTL reached the instrument threshold. It is the largest cell in every arm here and must never be read as evidence against a locus.

## Agreement with coloc, primary confirmatory family

| analysis | modality | probes | instrumented | threshold | `coloc_and_smr_heidi_not_rejected` | `coloc_and_smr_heidi_rejected` | `coloc_and_smr_heidi_untestable` | `coloc_smr_not_significant` | `coloc_no_instrument` | `smr_without_coloc` | `smr_coloc_not_scored` | `neither_tested_null` | `neither_no_instrument` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | A_g | 30 | 8 | 0.006 | 3 | 0 | 0 | 0 | 1 | 1 | 0 | 4 | 21 |
| aging__ad | S_g | 26 | 2 | 0.025 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 24 |
| aging__als | A_g | 27 | 10 | 0.005 | 5 | 0 | 0 | 0 | 2 | 1 | 0 | 4 | 15 |
| aging__als | S_g | 27 | 3 | 0.017 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 24 |
| aging__lbd | A_g | 3 | 0 | — | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 2 |
| aging__lbd | S_g | 3 | 0 | — | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| aging__pd | A_g | 18 | 5 | 0.010 | 2 | 0 | 0 | 0 | 0 | 1 | 0 | 2 | 13 |
| aging__pd | S_g | 15 | 0 | — | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 15 |
| aging__scz | A_g | 60 | 20 | 0.003 | 2 | 0 | 0 | 1 | 0 | 11 | 0 | 6 | 40 |
| aging__scz | S_g | 59 | 4 | 0.013 | 0 | 0 | 0 | 1 | 3 | 1 | 0 | 2 | 52 |

## SNP attrition into HEIDI

HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so what survives is worth stating. Per SMR run; the per-probe end of the trace is `nsnp_HEIDI` in the results table, itself capped at 20.

**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends past its locus loses the remainder. That is not a QC failure and says nothing about the QTL data.

| step | median across runs |
|---|---|
| BESD SNPs in the probe's cis window | 10,046 |
| ... found in the LD panel by rsID | 100.0% |
| ... with matching alleles | 100.0% |
| ... also carried by the GWAS | 4,595 |
| dropped by SMR beyond that (frequency check, duplicates) | 0 |

Genuine harmonization loss — SNPs present in all three and still dropped — peaks at **4** SNPs in any run. The identifier and allele steps are the QC ones, and both sit at 100%.

**Do not compare BESD-SNP retention across QTL sources.** Its denominator is imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far larger set against the same GWAS loci. The ratio measures panel density, not quality.

## GTEx nominations on the BrainSEQ switch axis

Per GTEx signal-level nomination, the `S_g` agreement class in each region (`*` = the region's tissue-matched GTEx tissue carries the GTEx sQTL call).

| gene | trait | caudate | dlpfc | hippocampus |
|---|---|---|---|---|
| ASB3 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| AZI2 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CDIP1 | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument`* |
| COPA | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| CTSH | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument`* |
| DOC2A | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| DOC2A | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| FAM221A | scz | `smr_without_coloc` | `neither_no_instrument` | `neither_no_instrument` |
| G2E3 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| GABBR2 | scz | `neither_no_instrument` | `neither_no_instrument` | — |
| GGNBP2 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| GPM6A | scz | `neither_no_instrument`* | `neither_tested_null` | `neither_no_instrument` |
| HMOX2 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| KLC1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NCOR1 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NDUFS3 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| NEK4 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| NT5C2 | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument` |
| PGS1 | als | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| PICALM | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PITPNM2 | pd | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PLCB2 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PPIL2 | scz | `coloc_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument` |
| PPIP5K1 | scz | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| PTPRN | als | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| RNASEH2C | scz | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| SCFD1 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| SIRPA | ad | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| SNCA | lbd | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument`* |
| SNCA | pd | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument`* |
| SPAG9 | ad | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument`* |
| SYT5 | scz | `neither_no_instrument` | `neither_no_instrument`* | `neither_no_instrument` |
| TMED4 | scz | `neither_no_instrument`* | `neither_no_instrument`* | `neither_no_instrument`* |
| TPCN1 | ad | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TPP1 | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| TTC19 | pd | `neither_no_instrument`* | `neither_no_instrument` | `neither_no_instrument` |
| TXNDC15 | als | `neither_no_instrument`* | `smr_without_coloc` | `smr_without_coloc` |
| UNC13A | als | `neither_no_instrument` | `neither_no_instrument` | `neither_no_instrument` |
| WIPI2 | als | `neither_tested_null` | `neither_no_instrument` | `neither_no_instrument` |
| ZNF232 | ad | `neither_no_instrument` | `neither_no_instrument` | — |

## Primary probes

| gene | trait | region | axis | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | pin | agreement |
|---|---|---|---|---|---|---|---|---|---|
| CTSH | ad | caudate | A_g | 0.990 (abf) | 0.123 (0.032) | 9.3e-05 | 0.089 (7) |  | `coloc_and_smr_heidi_not_rejected` |
| DOC2A | ad | hippocampus | A_g | 0.164 (susie) | 0.034 (0.061) | 0.580 | 0.010 (4) |  | `neither_tested_null` |
| NDUFS3 | ad | caudate | A_g | 0.000144 (susie) | -0.197 (0.054) | 0.000225 | 1.22e-05 (20) |  | `smr_without_coloc` |
| SIRPA | ad | caudate | A_g | 0.960 (abf) | 0.104 (0.025) | 4.02e-05 | 0.878 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | dlpfc | A_g | 0.961 (abf) | 0.166 (0.044) | 0.000143 | 0.835 (10) |  | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | hippocampus | A_g | 0.946 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| TMEM106B | ad | caudate | A_g | 0.010 (abf) | 0.113 (0.073) | 0.124 | 0.438 (20) |  | `neither_tested_null` |
| TMEM106B | ad | caudate | S_g | 0.741 (susie) | -0.029 (0.018) | 0.107 | 0.012 (20) | pinned | `neither_tested_null` |
| TMEM106B | ad | dlpfc | S_g | 0.923 (susie) | 0.068 (0.017) | 4.9e-05 | 0.172 (20) | pinned | `coloc_and_smr_heidi_not_rejected` |
| ZNF232 | ad | caudate | A_g | 0.046 (abf) | -0.045 (0.018) | 0.014 | 0.275 (6) |  | `neither_tested_null` |
| ZNF232 | ad | hippocampus | A_g | 1.71e-08 (abf) | -0.021 (0.015) | 0.153 | 0.238 (13) |  | `neither_tested_null` |
| G2E3 | als | dlpfc | A_g | 0.822 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| GGNBP2 | als | caudate | A_g | 0.947 (susie) | 0.252 (0.064) | 8.04e-05 | 0.559 (8) |  | `coloc_and_smr_heidi_not_rejected` |
| GGNBP2 | als | hippocampus | A_g | 0.829 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| PGS1 | als | caudate | A_g | 0.006 (abf) | 0.054 (0.048) | 0.265 | 0.355 (18) |  | `neither_tested_null` |
| PGS1 | als | dlpfc | A_g | 0.026 (abf) | -0.060 (0.028) | 0.030 | 0.002 (15) |  | `neither_tested_null` |
| PGS1 | als | hippocampus | A_g | 0.087 (abf) | -0.125 (0.061) | 0.039 | 0.002 (16) |  | `neither_tested_null` |
| PTPRN | als | caudate | A_g | 0.812 (susie) | -0.156 (0.036) | 1.69e-05 | 0.247 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | dlpfc | A_g | 0.029 (abf) | -0.013 (0.075) | 0.865 | 0.052 (7) |  | `neither_tested_null` |
| SCFD1 | als | caudate | A_g | 0.946 (susie) | 0.349 (0.070) | 6.96e-07 | 0.296 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SCFD1 | als | dlpfc | A_g | 0.936 (susie) | 0.197 (0.032) | 6.89e-10 | 0.444 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SCFD1 | als | hippocampus | A_g | 0.931 (susie) | 0.422 (0.079) | 8.61e-08 | 0.089 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| TXNDC15 | als | dlpfc | S_g | 0.772 (abf) | 0.086 (0.029) | 0.003 | 0.164 (5) | pinned | `smr_without_coloc` |
| TXNDC15 | als | hippocampus | A_g | 0.317 (abf) | -0.177 (0.057) | 0.002 | 0.189 (20) |  | `smr_without_coloc` |
| TXNDC15 | als | hippocampus | S_g | 0.662 (abf) | -0.095 (0.033) | 0.004 | 0.224 (7) | pinned, axis unstable | `smr_without_coloc` |
| WIPI2 | als | caudate | S_g | 0.019 (abf) | -0.017 (0.050) | 0.737 | 0.002 (20) | pinned | `neither_tested_null` |
| SNCA | lbd | caudate | A_g | 0.857 (abf) | — (—) | — | — (—) |  | `coloc_no_instrument` |
| AZI2 | pd | dlpfc | A_g | 0.017 (abf) | -0.043 (0.101) | 0.671 | 0.200 (18) |  | `neither_tested_null` |
| NCOR1 | pd | caudate | A_g | 0.891 (susie) | 0.444 (0.104) | 2.14e-05 | 0.058 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | caudate | A_g | 0.001 (abf) | 0.169 (0.091) | 0.064 | 0.359 (20) |  | `neither_tested_null` |
| SH3GL2 | pd | dlpfc | A_g | 0.964 (abf) | 0.198 (0.065) | 0.002 | 0.171 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| TTC19 | pd | caudate | A_g | 0.737 (susie) | -0.326 (0.081) | 5.9e-05 | 0.222 (18) |  | `smr_without_coloc` |
| COPA | scz | caudate | A_g | 0.200 (abf) | -0.104 (0.039) | 0.007 | 0.052 (18) |  | `neither_tested_null` |
| COPA | scz | dlpfc | A_g | 0.199 (abf) | -0.076 (0.025) | 0.002 | 0.037 (20) |  | `smr_without_coloc` |
| COPA | scz | hippocampus | A_g | 0.243 (abf) | -0.143 (0.058) | 0.014 | 0.062 (20) |  | `neither_tested_null` |
| DOC2A | scz | hippocampus | A_g | 0.051 (susie) | 0.173 (0.076) | 0.023 | 0.107 (4) |  | `neither_tested_null` |
| FAM221A | scz | caudate | A_g | 0.010 (susie) | 0.029 (0.009) | 0.001 | 0.010 (20) |  | `smr_without_coloc` |
| FAM221A | scz | caudate | S_g | 0.722 (susie) | -0.053 (0.012) | 1.69e-05 | 0.094 (20) | pinned | `smr_without_coloc` |
| FAM221A | scz | dlpfc | A_g | 0.009 (susie) | 0.042 (0.014) | 0.002 | 0.003 (20) |  | `smr_without_coloc` |
| FAM221A | scz | hippocampus | A_g | 0.010 (susie) | 0.048 (0.016) | 0.002 | 0.027 (20) |  | `smr_without_coloc` |
| GPM6A | scz | dlpfc | S_g | 0.000606 (susie) | 0.022 (0.017) | 0.194 | 0.108 (20) | pinned | `neither_tested_null` |
| KLC1 | scz | caudate | A_g | 1.13e-06 (susie) | -0.112 (0.036) | 0.002 | 6.32e-05 (20) |  | `smr_without_coloc` |
| KLC1 | scz | dlpfc | A_g | 0.003 (susie) | -0.086 (0.028) | 0.002 | 0.009 (20) |  | `smr_without_coloc` |
| KLC1 | scz | hippocampus | A_g | 0.001 (susie) | -0.141 (0.049) | 0.004 | 0.092 (16) |  | `neither_tested_null` |
| NEK4 | scz | caudate | A_g | 0.003 (abf) | 0.198 (0.042) | 2.36e-06 | 0.028 (20) |  | `smr_without_coloc` |
| NEK4 | scz | dlpfc | A_g | 0.368 (abf) | 0.142 (0.030) | 2.36e-06 | 0.030 (20) |  | `smr_without_coloc` |
| NEK4 | scz | hippocampus | A_g | 0.076 (abf) | 0.244 (0.053) | 4.58e-06 | 0.098 (20) |  | `smr_without_coloc` |
| NGEF | scz | caudate | S_g | 0.886 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| NSF | scz | caudate | S_g | 0.094 (abf) | 0.033 (0.026) | 0.211 | 0.033 (12) | pinned | `neither_tested_null` |
| NSF | scz | dlpfc | S_g | 0.804 (abf) | 0.033 (0.019) | 0.088 | 0.008 (16) | pinned | `coloc_smr_not_significant` |
| NT5C2 | scz | dlpfc | A_g | 0.845 (susie) | 0.027 (0.029) | 0.352 | 0.004 (20) |  | `coloc_smr_not_significant` |
| NT5C2 | scz | hippocampus | A_g | 0.530 (abf) | -0.203 (0.066) | 0.002 | 0.001 (20) |  | `smr_without_coloc` |
| PPIL2 | scz | caudate | S_g | 0.811 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |
| RNASEH2C | scz | caudate | A_g | 0.981 (susie) | -0.170 (0.046) | 0.000213 | 0.597 (12) |  | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | dlpfc | A_g | 0.983 (susie) | -0.145 (0.038) | 0.000162 | 0.569 (20) |  | `coloc_and_smr_heidi_not_rejected` |
| SYT5 | scz | dlpfc | A_g | 0.004 (abf) | 0.009 (0.038) | 0.811 | 0.003 (9) |  | `neither_tested_null` |
| TMED4 | scz | caudate | A_g | 0.356 (susie) | -0.139 (0.048) | 0.004 | 0.100 (11) |  | `neither_tested_null` |
| TMED4 | scz | dlpfc | A_g | 0.190 (susie) | -0.108 (0.035) | 0.002 | 0.109 (20) |  | `smr_without_coloc` |
| ZNF592 | scz | hippocampus | S_g | 0.839 (abf) | — (—) | — | — (—) | pinned | `coloc_no_instrument` |

Probes with neither an SMR instrument nor a coloc call are in `smr_results.parquet` and not listed here.

