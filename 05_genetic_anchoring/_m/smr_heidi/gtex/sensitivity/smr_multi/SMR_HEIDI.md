# SMR + HEIDI on the signal-level colocalization nominations

QTL source: `gtex`. Instrument threshold `--peqtl-smr 5e-08`. Multi-SNP SMR (`--smr-multi`, LD pruned at r2 0.1) -- **SENSITIVITY ARM**. HEIDI SNPs at p < 0.0015654, 3-20 SNPs, 2000 kb cis window. Significance: Bonferroni at 0.05 over the instrumented probes of the row's OWN family (primary confirmatory / secondary event-localization), never pooled across the two. HEIDI rejects a single shared variant at p < 0.01. LD reference: 1000G EUR.

## How these numbers may be read

- `b_SMR` is a signed ratio estimate relating genetically predicted phenotype to disease under a single-causal-variant, no-pleiotropy model. It does **not** establish causal direction and cannot distinguish causality from horizontal pleiotropy.
- An sQTL probe is a LeafCutter intron-excision ratio, which is compositional within its cluster: introns sharing a splice site trade usage, so sibling probes carry opposite `b_SMR` signs by construction. A sign is read relative to its cluster, never alone.
- Failing to reject HEIDI is **not** evidence of a shared variant; HEIDI is underpowered at GTEx brain sample sizes. It reads "not rejected".
- A HEIDI rejection does **not** overrule a strong signal-level colocalization, and SMR significance does not promote a locus coloc did not support. Disagreements stay disagreements.

## What was testable, by family

Two families, corrected apart. The **primary confirmatory family** holds one pre-designated probe per gene; a gene's other introns form a **secondary event-localization family** with its own Bonferroni correction, so the primary threshold is not inflated by introns and the introns do not escape correction when they are discussed. `F` is the instrument strength `(b_eQTL/se_eQTL)^2` of the top cis-QTL SNP, reported on every row and never used to exclude one. It is summarized by its 5th percentile rather than a count below the conventional F < 10: an instrument that clears p < 5e-8 has |z| ≳ 5.4 and so F ≳ 30 (≈ 24 at the relaxed 1e-6 arm), so a weak-instrument count is zero by construction and carries no information; the per-row `weak_instrument` flag stays in `smr_results.parquet` for any run at a looser threshold.

| analysis | modality | family | probes | instrumented | threshold | F median [min-max] | F p5 | `no_instrument` | `instrumented_tested_null` | `smr_signal_heidi_unavailable` | `smr_signal_heidi_rejects` | `smr_heidi_supported` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | eQTL | primary | 63 | 14 | 0.004 | 55 [30-971] | 30 | 49 | 8 | 0 | 0 | 6 |
| aging__ad | sQTL | primary | 63 | 49 | 0.001 | 67 [30-531] | 32 | 14 | 4 | 0 | 1 | 44 |
| aging__ad | sQTL | secondary | 814 | 86 | 0.000581 | 67 [30-524] | 31 | 728 | 39 | 0 | 1 | 46 |
| aging__als | eQTL | primary | 45 | 18 | 0.003 | 94 [37-154] | 46 | 27 | 1 | 0 | 11 | 6 |
| aging__als | sQTL | primary | 45 | 39 | 0.001 | 50 [30-349] | 32 | 6 | 4 | 0 | 2 | 33 |
| aging__als | sQTL | secondary | 602 | 62 | 0.000806 | 68 [31-910] | 33 | 540 | 34 | 0 | 4 | 24 |
| aging__lbd | eQTL | primary | 11 | 0 | — | — | — | 11 | 0 | 0 | 0 | 0 |
| aging__lbd | sQTL | primary | 11 | 11 | 0.005 | 105 [39-221] | 57 | 0 | 0 | 0 | 0 | 11 |
| aging__lbd | sQTL | secondary | 81 | 21 | 0.002 | 66 [32-203] | 33 | 60 | 9 | 0 | 11 | 1 |
| aging__pd | eQTL | primary | 31 | 20 | 0.003 | 96 [32-630] | 33 | 11 | 2 | 0 | 0 | 18 |
| aging__pd | sQTL | primary | 31 | 27 | 0.002 | 69 [30-309] | 32 | 4 | 0 | 2 | 0 | 25 |
| aging__pd | sQTL | secondary | 384 | 50 | 0.001 | 150 [30-964] | 33 | 334 | 11 | 2 | 0 | 37 |
| aging__scz | eQTL | primary | 228 | 86 | 0.000581 | 61 [30-830] | 32 | 142 | 32 | 0 | 8 | 46 |
| aging__scz | sQTL | primary | 228 | 164 | 0.000305 | 60 [30-280] | 32 | 64 | 25 | 0 | 9 | 130 |
| aging__scz | sQTL | secondary | 2605 | 242 | 0.000207 | 54 [30-293] | 31 | 2363 | 137 | 0 | 2 | 103 |
| brainseq-sczd__scz | eQTL | primary | 32 | 16 | 0.003 | 45 [31-344] | 34 | 16 | 3 | 0 | 0 | 13 |
| brainseq-sczd__scz | sQTL | primary | 32 | 23 | 0.002 | 57 [30-225] | 32 | 9 | 2 | 0 | 1 | 20 |
| brainseq-sczd__scz | sQTL | secondary | 321 | 38 | 0.001 | 48 [30-177] | 30 | 283 | 13 | 0 | 0 | 25 |

**`no_instrument` is not a negative result** — the probe was never tested, because no cis-QTL reached the instrument threshold. It is the largest cell in every arm here and must never be read as evidence against a locus.

## Agreement with coloc, primary confirmatory family

| analysis | modality | probes | instrumented | threshold | `coloc_and_smr_heidi_not_rejected` | `coloc_and_smr_heidi_rejected` | `coloc_and_smr_heidi_untestable` | `coloc_smr_not_significant` | `coloc_no_instrument` | `smr_without_coloc` | `smr_coloc_not_scored` | `neither_tested_null` | `neither_no_instrument` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | eQTL | 63 | 14 | 0.004 | 2 | 0 | 0 | 0 | 8 | 4 | 0 | 8 | 41 |
| aging__ad | sQTL | 63 | 49 | 0.001 | 44 | 1 | 0 | 4 | 14 | 0 | 0 | 0 | 0 |
| aging__als | eQTL | 45 | 18 | 0.003 | 2 | 2 | 0 | 0 | 3 | 13 | 0 | 1 | 24 |
| aging__als | sQTL | 45 | 39 | 0.001 | 33 | 2 | 0 | 4 | 6 | 0 | 0 | 0 | 0 |
| aging__lbd | eQTL | 11 | 0 | — | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 11 |
| aging__lbd | sQTL | 11 | 11 | 0.005 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| aging__pd | eQTL | 31 | 20 | 0.003 | 10 | 0 | 0 | 0 | 1 | 8 | 0 | 2 | 10 |
| aging__pd | sQTL | 31 | 27 | 0.002 | 25 | 0 | 2 | 0 | 4 | 0 | 0 | 0 | 0 |
| aging__scz | eQTL | 228 | 86 | 0.000581 | 27 | 3 | 0 | 3 | 16 | 24 | 0 | 29 | 126 |
| aging__scz | sQTL | 228 | 164 | 0.000305 | 130 | 9 | 0 | 25 | 64 | 0 | 0 | 0 | 0 |
| brainseq-sczd__scz | eQTL | 32 | 16 | 0.003 | 11 | 0 | 0 | 0 | 5 | 2 | 0 | 3 | 11 |
| brainseq-sczd__scz | sQTL | 32 | 23 | 0.002 | 20 | 1 | 0 | 2 | 9 | 0 | 0 | 0 | 0 |

### Secondary family (event localization)

A gene's non-primary introns, corrected within their own family. These localize an event; they are not additional confirmatory evidence for a locus.

| analysis | modality | probes | instrumented | threshold | `smr_heidi_supported` | `smr_signal_heidi_rejects` |
|---|---|---|---|---|---|---|
| aging__ad | sQTL | 814 | 86 | 0.000581 | 46 | 1 |
| aging__als | sQTL | 602 | 62 | 0.000806 | 24 | 4 |
| aging__lbd | sQTL | 81 | 21 | 0.002 | 1 | 11 |
| aging__pd | sQTL | 384 | 50 | 0.001 | 37 | 0 |
| aging__scz | sQTL | 2605 | 242 | 0.000207 | 103 | 2 |
| brainseq-sczd__scz | sQTL | 321 | 38 | 0.001 | 25 | 0 |

## SNP attrition into HEIDI

HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so what survives is worth stating. Per SMR run; the per-probe end of the trace is `nsnp_HEIDI` in the results table, itself capped at 20.

**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends past its locus loses the remainder. That is not a QC failure and says nothing about the QTL data.

| step | median across runs |
|---|---|
| BESD SNPs in the probe's cis window | 8,987 |
| ... found in the LD panel by rsID | 100.0% |
| ... with matching alleles | 100.0% |
| ... also carried by the GWAS | 4,407 |
| dropped by SMR beyond that (frequency check, duplicates) | 0 |

Genuine harmonization loss — SNPs present in all three and still dropped — peaks at **7** SNPs in any run. The identifier and allele steps are the QC ones, and both sit at 100%.

**Do not compare BESD-SNP retention across QTL sources.** Its denominator is imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far larger set against the same GWAS loci. The ratio measures panel density, not quality.

## Primary probes

| gene | trait | tissue | modality | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | agreement |
|---|---|---|---|---|---|---|---|---|
| AKT1 | ad | Brain_Caudate_basal_ganglia | eQTL | 0.157 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| AKT1 | ad | Brain_Caudate_basal_ganglia | sQTL | 0.840 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| BCKDK | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.258 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| BCKDK | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.994 (susie) | -0.080 (0.017) | 2.01e-06 | 0.545 (20) | `coloc_and_smr_heidi_not_rejected` |
| BCKDK | ad | Brain_Cerebellum | eQTL | 0.117 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| BCKDK | ad | Brain_Cerebellum | sQTL | 0.995 (susie) | -0.083 (0.017) | 1.18e-06 | 0.343 (20) | `coloc_and_smr_heidi_not_rejected` |
| BCKDK | ad | Brain_Frontal_Cortex_BA9 | eQTL | 0.030 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| BCKDK | ad | Brain_Frontal_Cortex_BA9 | sQTL | 0.993 (susie) | -0.115 (0.028) | 3.6e-05 | 0.253 (14) | `coloc_and_smr_heidi_not_rejected` |
| BCKDK | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.050 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| BCKDK | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.995 (susie) | -0.103 (0.023) | 1.09e-05 | 0.330 (20) | `coloc_and_smr_heidi_not_rejected` |
| COG7 | ad | Brain_Amygdala | eQTL | 0.030 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| COG7 | ad | Brain_Amygdala | sQTL | 0.837 (abf) | -0.057 (0.018) | 0.001 | 0.513 (20) | `coloc_smr_not_significant` |
| DOC2A | ad | Brain_Amygdala | eQTL | 0.030 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | ad | Brain_Amygdala | sQTL | 0.881 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| IFNAR2 | ad | Brain_Caudate_basal_ganglia | eQTL | 0.029 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IFNAR2 | ad | Brain_Caudate_basal_ganglia | sQTL | 0.855 (abf) | -0.079 (0.024) | 0.001 | 0.055 (4) | `coloc_smr_not_significant` |
| IFNAR2 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.157 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IFNAR2 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.861 (abf) | 0.061 (0.017) | 0.000251 | 0.014 (10) | `coloc_and_smr_heidi_not_rejected` |
| IFNAR2 | ad | Brain_Cerebellum | eQTL | 0.626 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IFNAR2 | ad | Brain_Cerebellum | sQTL | 0.864 (abf) | -0.054 (0.015) | 0.000255 | 0.011 (14) | `coloc_and_smr_heidi_not_rejected` |
| IFNAR2 | ad | Brain_Hippocampus | eQTL | 0.043 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IFNAR2 | ad | Brain_Hippocampus | sQTL | 0.840 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| IFNAR2 | ad | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IFNAR2 | ad | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.870 (abf) | 0.076 (0.024) | 0.001 | 0.491 (6) | `coloc_smr_not_significant` |
| IFNAR2 | ad | Brain_Putamen_basal_ganglia | eQTL | 0.059 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IFNAR2 | ad | Brain_Putamen_basal_ganglia | sQTL | 0.843 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| IFNAR2 | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.829 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| IFNAR2 | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.865 (abf) | 0.079 (0.024) | 0.000946 | 0.004 (6) | `coloc_and_smr_heidi_rejected` |
| INO80E | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.679 (susie) | 0.125 (0.024) | 2e-07 | 0.135 (20) | `smr_without_coloc` |
| INO80E | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.815 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | ad | Brain_Cerebellum | eQTL | 0.568 (susie) | 0.111 (0.021) | 1.76e-07 | 0.048 (20) | `smr_without_coloc` |
| INO80E | ad | Brain_Cerebellum | sQTL | 0.903 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | ad | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.774 (susie) | — (—) | — | — (—) | `neither_no_instrument` |
| INO80E | ad | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.832 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| INTS8 | ad | Brain_Cerebellum | eQTL | 9.18e-06 (susie) | -0.036 (0.016) | 0.027 | 0.059 (20) | `neither_tested_null` |
| INTS8 | ad | Brain_Cerebellum | sQTL | 0.931 (susie) | 0.066 (0.014) | 2.18e-06 | 0.186 (20) | `coloc_and_smr_heidi_not_rejected` |
| INTS8 | ad | Brain_Hippocampus | eQTL | 0.086 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| INTS8 | ad | Brain_Hippocampus | sQTL | 0.925 (abf) | -0.069 (0.020) | 0.00079 | 0.081 (20) | `coloc_and_smr_heidi_not_rejected` |
| INTS8 | ad | Brain_Hypothalamus | eQTL | 0.012 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| INTS8 | ad | Brain_Hypothalamus | sQTL | 0.828 (susie) | 0.061 (0.018) | 0.000613 | 0.034 (20) | `coloc_and_smr_heidi_not_rejected` |
| INTS8 | ad | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.078 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| INTS8 | ad | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.940 (susie) | 0.061 (0.018) | 0.000508 | 0.073 (20) | `coloc_and_smr_heidi_not_rejected` |
| ITGB1BP1 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.215 (abf) | -0.099 (0.030) | 0.000765 | 0.073 (14) | `smr_without_coloc` |
| ITGB1BP1 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.979 (susie) | -0.131 (0.032) | 4e-05 | 0.273 (18) | `coloc_and_smr_heidi_not_rejected` |
| ITGB1BP1 | ad | Brain_Cerebellum | eQTL | 0.203 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ITGB1BP1 | ad | Brain_Cerebellum | sQTL | 0.977 (susie) | -0.113 (0.028) | 6.3e-05 | 0.682 (15) | `coloc_and_smr_heidi_not_rejected` |
| ITGB1BP1 | ad | Brain_Cortex | eQTL | 0.040 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ITGB1BP1 | ad | Brain_Cortex | sQTL | 0.955 (susie) | -0.098 (0.025) | 0.000124 | 0.287 (8) | `coloc_and_smr_heidi_not_rejected` |
| ITGB1BP1 | ad | Brain_Frontal_Cortex_BA9 | eQTL | 0.058 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ITGB1BP1 | ad | Brain_Frontal_Cortex_BA9 | sQTL | 0.951 (susie) | 0.090 (0.023) | 7.98e-05 | 0.036 (18) | `coloc_and_smr_heidi_not_rejected` |
| ITGB1BP1 | ad | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.079 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ITGB1BP1 | ad | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.947 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| NDUFS3 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.000281 (abf) | -0.162 (0.061) | 0.008 | 0.003 (20) | `neither_tested_null` |
| NDUFS3 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.864 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PICALM | ad | Brain_Cortex | eQTL | 0.117 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PICALM | ad | Brain_Cortex | sQTL | 0.818 (abf) | -0.309 (0.059) | 1.49e-07 | 0.170 (20) | `coloc_and_smr_heidi_not_rejected` |
| PILRB | ad | Brain_Putamen_basal_ganglia | eQTL | 2.79e-14 (abf) | 0.004 (0.011) | 0.715 | 0.086 (20) | `neither_tested_null` |
| PILRB | ad | Brain_Putamen_basal_ganglia | sQTL | 0.999 (susie) | -0.005 (0.014) | 0.722 | 0.002 (20) | `coloc_smr_not_significant` |
| RAD51C | ad | Brain_Cortex | eQTL | 0.00033 (susie) | 0.024 (0.013) | 0.074 | 0.001 (20) | `neither_tested_null` |
| RAD51C | ad | Brain_Cortex | sQTL | 0.992 (susie) | -0.087 (0.022) | 5.98e-05 | 0.721 (12) | `coloc_and_smr_heidi_not_rejected` |
| SERPINB1 | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.024 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SERPINB1 | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.917 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SIRPA | ad | Brain_Amygdala | eQTL | 0.106 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Amygdala | sQTL | 0.960 (abf) | 0.031 (0.007) | 1.09e-05 | 0.954 (13) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.914 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SIRPA | ad | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.965 (abf) | 0.029 (0.007) | 1.52e-05 | 0.929 (18) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Caudate_basal_ganglia | eQTL | 0.521 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Caudate_basal_ganglia | sQTL | 0.963 (abf) | 0.031 (0.007) | 7.81e-06 | 0.922 (18) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.905 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SIRPA | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.965 (abf) | 0.029 (0.007) | 8.57e-06 | 0.951 (16) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Cerebellum | eQTL | 0.947 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SIRPA | ad | Brain_Cerebellum | sQTL | 0.965 (abf) | 0.031 (0.007) | 8.42e-06 | 0.921 (16) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Cortex | eQTL | 0.936 (abf) | 0.201 (0.054) | 0.0002 | 0.978 (8) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Cortex | sQTL | 0.957 (abf) | 0.028 (0.007) | 1.93e-05 | 0.942 (16) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Frontal_Cortex_BA9 | eQTL | 0.939 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SIRPA | ad | Brain_Frontal_Cortex_BA9 | sQTL | 0.964 (abf) | 0.032 (0.007) | 8.35e-06 | 0.880 (16) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Hippocampus | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Hippocampus | sQTL | 0.950 (abf) | 0.030 (0.007) | 2.26e-05 | 0.910 (16) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Hypothalamus | eQTL | 0.172 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Hypothalamus | sQTL | 0.961 (abf) | -0.030 (0.007) | 9.01e-06 | 0.931 (14) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.507 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.936 (abf) | 0.029 (0.007) | 2.07e-05 | 0.925 (17) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Putamen_basal_ganglia | eQTL | 0.922 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SIRPA | ad | Brain_Putamen_basal_ganglia | sQTL | 0.965 (abf) | 0.031 (0.007) | 8.66e-06 | 0.917 (17) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.073 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.965 (abf) | 0.030 (0.007) | 1.02e-05 | 0.997 (13) | `coloc_and_smr_heidi_not_rejected` |
| SIRPA | ad | Brain_Substantia_nigra | eQTL | 0.042 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SIRPA | ad | Brain_Substantia_nigra | sQTL | 0.955 (abf) | 0.030 (0.007) | 2.81e-05 | 0.900 (13) | `coloc_and_smr_heidi_not_rejected` |
| SLC39A13 | ad | Brain_Hippocampus | eQTL | 0.007 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SLC39A13 | ad | Brain_Hippocampus | sQTL | 0.983 (abf) | -0.062 (0.012) | 5.74e-07 | 0.011 (20) | `coloc_and_smr_heidi_not_rejected` |
| SLC39A13 | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.042 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SLC39A13 | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.987 (abf) | -0.095 (0.017) | 2.07e-08 | 0.195 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Amygdala | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SPAG9 | ad | Brain_Amygdala | sQTL | 0.811 (abf) | 0.050 (0.014) | 0.000235 | 0.499 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.014 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SPAG9 | ad | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.823 (abf) | 0.080 (0.023) | 0.000563 | 0.432 (16) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Caudate_basal_ganglia | eQTL | 0.011 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SPAG9 | ad | Brain_Caudate_basal_ganglia | sQTL | 0.809 (abf) | 0.053 (0.014) | 0.000201 | 0.367 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Cortex | eQTL | 0.025 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SPAG9 | ad | Brain_Cortex | sQTL | 0.802 (abf) | -0.071 (0.021) | 0.000899 | 0.107 (12) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Hippocampus | eQTL | 0.023 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SPAG9 | ad | Brain_Hippocampus | sQTL | 0.823 (abf) | 0.039 (0.011) | 0.00026 | 0.258 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Hypothalamus | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SPAG9 | ad | Brain_Hypothalamus | sQTL | 0.839 (abf) | 0.054 (0.014) | 0.000144 | 0.557 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Putamen_basal_ganglia | eQTL | 0.005 (abf) | 0.080 (0.064) | 0.215 | 0.131 (9) | `neither_tested_null` |
| SPAG9 | ad | Brain_Putamen_basal_ganglia | sQTL | 0.827 (abf) | 0.046 (0.013) | 0.00038 | 0.104 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPAG9 | ad | Brain_Substantia_nigra | eQTL | 0.006 (abf) | 0.052 (0.042) | 0.215 | 0.080 (7) | `neither_tested_null` |
| SPAG9 | ad | Brain_Substantia_nigra | sQTL | 0.821 (abf) | 0.040 (0.011) | 0.000315 | 0.187 (20) | `coloc_and_smr_heidi_not_rejected` |
| SPI1 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.881 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SPI1 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.926 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TPCN1 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.819 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| TPCN1 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.861 (susie) | 0.051 (0.013) | 9.84e-05 | 0.343 (20) | `coloc_and_smr_heidi_not_rejected` |
| TPCN1 | ad | Brain_Cerebellum | eQTL | 0.243 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TPCN1 | ad | Brain_Cerebellum | sQTL | 0.836 (susie) | -0.062 (0.017) | 0.000207 | 0.745 (10) | `coloc_and_smr_heidi_not_rejected` |
| VWA5B2 | ad | Brain_Amygdala | eQTL | 0.084 (abf) | 0.034 (0.012) | 0.004 | 0.012 (20) | `neither_tested_null` |
| VWA5B2 | ad | Brain_Amygdala | sQTL | 0.864 (abf) | 0.058 (0.017) | 0.00071 | 0.156 (15) | `coloc_and_smr_heidi_not_rejected` |
| YPEL3 | ad | Brain_Caudate_basal_ganglia | eQTL | 0.556 (susie) | -0.293 (0.075) | 9.11e-05 | 0.174 (9) | `smr_without_coloc` |
| YPEL3 | ad | Brain_Caudate_basal_ganglia | sQTL | 0.822 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YPEL3 | ad | Brain_Frontal_Cortex_BA9 | eQTL | 0.932 (susie) | -0.293 (0.076) | 0.000119 | 0.266 (10) | `coloc_and_smr_heidi_not_rejected` |
| YPEL3 | ad | Brain_Frontal_Cortex_BA9 | sQTL | 0.845 (abf) | -0.110 (0.026) | 1.81e-05 | 0.262 (13) | `coloc_and_smr_heidi_not_rejected` |
| YPEL3 | ad | Brain_Putamen_basal_ganglia | eQTL | 0.765 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL3 | ad | Brain_Putamen_basal_ganglia | sQTL | 0.862 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YPEL3 | ad | Brain_Substantia_nigra | eQTL | 0.539 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL3 | ad | Brain_Substantia_nigra | sQTL | 0.935 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ZNF232 | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.009 (abf) | -0.039 (0.016) | 0.012 | 0.711 (9) | `neither_tested_null` |
| ZNF232 | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.899 (abf) | -0.073 (0.017) | 1.17e-05 | 0.020 (14) | `coloc_and_smr_heidi_not_rejected` |
| C9orf72 | als | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 7.42e-10 (abf) | 0.138 (0.032) | 1.84e-05 | 0.000162 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.910 (abf) | -0.263 (0.050) | 1.44e-07 | 0.069 (20) | `coloc_and_smr_heidi_not_rejected` |
| C9orf72 | als | Brain_Caudate_basal_ganglia | eQTL | 2.1e-24 (abf) | 0.136 (0.029) | 2.8e-06 | 1.76e-05 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Caudate_basal_ganglia | sQTL | 0.966 (abf) | -0.249 (0.040) | 3.22e-10 | 0.284 (20) | `coloc_and_smr_heidi_not_rejected` |
| C9orf72 | als | Brain_Cerebellar_Hemisphere | eQTL | 2.44e-33 (abf) | 0.109 (0.023) | 2.56e-06 | 1.36e-08 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Cerebellar_Hemisphere | sQTL | 0.988 (abf) | -0.157 (0.018) | 6.53e-19 | 0.156 (20) | `coloc_and_smr_heidi_not_rejected` |
| C9orf72 | als | Brain_Cerebellum | eQTL | 2.53e-33 (abf) | 0.094 (0.020) | 3.54e-06 | 1.44e-06 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Cerebellum | sQTL | 0.991 (abf) | -0.235 (0.028) | 1.48e-16 | 0.005 (20) | `coloc_and_smr_heidi_rejected` |
| C9orf72 | als | Brain_Cortex | eQTL | 8.54e-20 (abf) | 0.134 (0.029) | 4.44e-06 | 2.64e-05 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Cortex | sQTL | 0.929 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| C9orf72 | als | Brain_Frontal_Cortex_BA9 | eQTL | 1.65e-15 (abf) | 0.152 (0.035) | 1.35e-05 | 9.92e-06 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Frontal_Cortex_BA9 | sQTL | 0.831 (abf) | -0.176 (0.034) | 2.47e-07 | 0.008 (20) | `coloc_and_smr_heidi_rejected` |
| C9orf72 | als | Brain_Hippocampus | eQTL | 2.9e-16 (abf) | 0.168 (0.036) | 4.06e-06 | 1.93e-05 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Hippocampus | sQTL | 0.984 (abf) | -0.264 (0.045) | 5.48e-09 | 0.484 (20) | `coloc_and_smr_heidi_not_rejected` |
| C9orf72 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 3.45e-09 (abf) | 0.196 (0.028) | 5.8e-12 | 5.06e-06 (20) | `smr_without_coloc` |
| C9orf72 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.984 (abf) | -0.220 (0.032) | 4.68e-12 | 0.058 (20) | `coloc_and_smr_heidi_not_rejected` |
| FNBP1 | als | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.017 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FNBP1 | als | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.844 (abf) | -0.049 (0.012) | 5.92e-05 | 0.208 (20) | `coloc_and_smr_heidi_not_rejected` |
| FNBP1 | als | Brain_Cortex | eQTL | 0.034 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FNBP1 | als | Brain_Cortex | sQTL | 0.882 (abf) | -0.042 (0.010) | 3.86e-05 | 0.226 (20) | `coloc_and_smr_heidi_not_rejected` |
| FNBP1 | als | Brain_Hippocampus | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FNBP1 | als | Brain_Hippocampus | sQTL | 0.929 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| FNBP1 | als | Brain_Hypothalamus | eQTL | 0.020 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FNBP1 | als | Brain_Hypothalamus | sQTL | 0.821 (abf) | -0.047 (0.012) | 8.74e-05 | 0.398 (20) | `coloc_and_smr_heidi_not_rejected` |
| FNBP1 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.099 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FNBP1 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.882 (abf) | -0.047 (0.011) | 4.03e-05 | 0.228 (20) | `coloc_and_smr_heidi_not_rejected` |
| G2E3 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.012 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| G2E3 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.839 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GGNBP2 | als | Brain_Cerebellum | eQTL | 0.383 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GGNBP2 | als | Brain_Cerebellum | sQTL | 0.912 (susie) | 0.128 (0.035) | 0.000231 | 0.117 (13) | `coloc_and_smr_heidi_not_rejected` |
| NME4 | als | Brain_Cortex | eQTL | 0.289 (abf) | 0.195 (0.062) | 0.002 | 0.329 (20) | `smr_without_coloc` |
| NME4 | als | Brain_Cortex | sQTL | 0.866 (abf) | -0.116 (0.032) | 0.000285 | 0.343 (20) | `coloc_and_smr_heidi_not_rejected` |
| NSMAF | als | Brain_Cerebellar_Hemisphere | eQTL | 0.951 (abf) | -0.167 (0.047) | 0.000393 | 0.004 (13) | `coloc_and_smr_heidi_rejected` |
| NSMAF | als | Brain_Cerebellar_Hemisphere | sQTL | 0.867 (abf) | 0.097 (0.031) | 0.001 | 0.320 (6) | `coloc_smr_not_significant` |
| NSMAF | als | Brain_Cerebellum | eQTL | 0.939 (abf) | -0.154 (0.044) | 0.000417 | 0.003 (14) | `coloc_and_smr_heidi_rejected` |
| NSMAF | als | Brain_Cerebellum | sQTL | 0.940 (abf) | 0.091 (0.026) | 0.000452 | 0.127 (10) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.826 (abf) | -0.074 (0.021) | 0.000423 | 0.814 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Cortex | eQTL | 0.010 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Cortex | sQTL | 0.964 (abf) | -0.059 (0.015) | 9.19e-05 | 0.882 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Hippocampus | eQTL | 0.033 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Hippocampus | sQTL | 0.900 (abf) | -0.097 (0.029) | 0.000879 | 0.616 (18) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.027 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.945 (abf) | -0.059 (0.015) | 8.36e-05 | 0.887 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Putamen_basal_ganglia | eQTL | 0.038 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Putamen_basal_ganglia | sQTL | 0.856 (abf) | -0.124 (0.035) | 0.000415 | 0.914 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Substantia_nigra | eQTL | 0.153 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Substantia_nigra | sQTL | 0.947 (abf) | -0.097 (0.028) | 0.000463 | 0.891 (19) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cerebellar_Hemisphere | eQTL | 0.031 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cerebellar_Hemisphere | sQTL | 0.880 (susie) | 0.121 (0.034) | 0.000323 | 0.564 (8) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cerebellum | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cerebellum | sQTL | 0.881 (susie) | -0.099 (0.026) | 0.000113 | 0.962 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cortex | eQTL | 0.585 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cortex | sQTL | 0.879 (susie) | -0.120 (0.031) | 0.000118 | 0.663 (12) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Frontal_Cortex_BA9 | eQTL | 0.074 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Frontal_Cortex_BA9 | sQTL | 0.832 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PTPRN | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.964 (abf) | -0.172 (0.043) | 5.97e-05 | 0.444 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.918 (abf) | 0.193 (0.050) | 0.000125 | 0.618 (16) | `coloc_and_smr_heidi_not_rejected` |
| RCSD1 | als | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.013 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RCSD1 | als | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.873 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RPS6KL1 | als | Brain_Cerebellar_Hemisphere | eQTL | 0.932 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RPS6KL1 | als | Brain_Cerebellar_Hemisphere | sQTL | 0.858 (abf) | 0.038 (0.016) | 0.017 | 0.123 (20) | `coloc_smr_not_significant` |
| RPS6KL1 | als | Brain_Cerebellum | eQTL | 0.043 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RPS6KL1 | als | Brain_Cerebellum | sQTL | 0.912 (abf) | 0.058 (0.015) | 7.62e-05 | 0.302 (20) | `coloc_and_smr_heidi_not_rejected` |
| RPS6KL1 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.403 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RPS6KL1 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.829 (abf) | 0.070 (0.029) | 0.015 | 0.320 (20) | `coloc_smr_not_significant` |
| SCFD1 | als | Brain_Hypothalamus | eQTL | 0.929 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| SCFD1 | als | Brain_Hypothalamus | sQTL | 0.856 (abf) | -0.184 (0.042) | 9.15e-06 | 0.401 (17) | `coloc_and_smr_heidi_not_rejected` |
| TMEM175 | als | Brain_Cerebellar_Hemisphere | eQTL | 0.011 (susie) | -0.149 (0.046) | 0.001 | 0.080 (20) | `smr_without_coloc` |
| TMEM175 | als | Brain_Cerebellar_Hemisphere | sQTL | 0.999 (susie) | 0.150 (0.039) | 9.93e-05 | 0.808 (20) | `coloc_and_smr_heidi_not_rejected` |
| TMEM175 | als | Brain_Cerebellum | eQTL | 0.012 (susie) | -0.131 (0.044) | 0.003 | 0.155 (20) | `neither_tested_null` |
| TMEM175 | als | Brain_Cerebellum | sQTL | 0.997 (susie) | -0.149 (0.041) | 0.000314 | 0.362 (18) | `coloc_and_smr_heidi_not_rejected` |
| TMEM175 | als | Brain_Cortex | eQTL | 0.012 (susie) | -0.107 (0.033) | 0.001 | 0.005 (20) | `smr_without_coloc` |
| TMEM175 | als | Brain_Cortex | sQTL | 0.998 (susie) | -0.091 (0.022) | 4.88e-05 | 0.426 (20) | `coloc_and_smr_heidi_not_rejected` |
| TNFSF13 | als | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.901 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TNFSF13 | als | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.906 (abf) | 0.084 (0.025) | 0.00073 | 0.692 (11) | `coloc_and_smr_heidi_not_rejected` |
| TNFSF13 | als | Brain_Putamen_basal_ganglia | eQTL | 0.912 (abf) | -0.296 (0.084) | 0.000388 | 0.091 (9) | `coloc_and_smr_heidi_not_rejected` |
| TNFSF13 | als | Brain_Putamen_basal_ganglia | sQTL | 0.911 (abf) | 0.090 (0.026) | 0.000537 | 0.578 (20) | `coloc_and_smr_heidi_not_rejected` |
| TPP1 | als | Brain_Cerebellum | eQTL | 0.774 (abf) | 0.213 (0.064) | 0.000948 | 0.538 (7) | `smr_without_coloc` |
| TPP1 | als | Brain_Cerebellum | sQTL | 0.996 (abf) | -0.124 (0.030) | 2.99e-05 | 0.950 (11) | `coloc_and_smr_heidi_not_rejected` |
| TXNDC15 | als | Brain_Caudate_basal_ganglia | eQTL | 0.347 (abf) | -0.109 (0.034) | 0.001 | 0.066 (20) | `smr_without_coloc` |
| TXNDC15 | als | Brain_Caudate_basal_ganglia | sQTL | 0.969 (abf) | -0.113 (0.029) | 0.000119 | 0.346 (20) | `coloc_and_smr_heidi_not_rejected` |
| UNC13A | als | Brain_Cerebellar_Hemisphere | eQTL | 0.377 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| UNC13A | als | Brain_Cerebellar_Hemisphere | sQTL | 0.936 (susie) | -0.184 (0.028) | 9.5e-11 | 0.205 (6) | `coloc_and_smr_heidi_not_rejected` |
| UNC13A | als | Brain_Cerebellum | eQTL | 0.090 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| UNC13A | als | Brain_Cerebellum | sQTL | 0.959 (susie) | -0.165 (0.026) | 9.03e-11 | 0.166 (6) | `coloc_and_smr_heidi_not_rejected` |
| WHAMM | als | Brain_Hippocampus | eQTL | 0.043 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| WHAMM | als | Brain_Hippocampus | sQTL | 0.916 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| WIPI2 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.020 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| WIPI2 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.893 (abf) | -0.097 (0.031) | 0.002 | 0.014 (14) | `coloc_smr_not_significant` |
| SNCA | lbd | Brain_Amygdala | eQTL | 0.113 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Amygdala | sQTL | 0.955 (susie) | 0.245 (0.044) | 2.93e-08 | 0.866 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.103 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.951 (susie) | 0.218 (0.037) | 3.43e-09 | 0.177 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Cerebellar_Hemisphere | eQTL | 0.009 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Cerebellar_Hemisphere | sQTL | 0.952 (susie) | 0.267 (0.047) | 1.28e-08 | 0.479 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Cerebellum | eQTL | 0.067 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Cerebellum | sQTL | 0.940 (susie) | 0.288 (0.055) | 1.35e-07 | 0.627 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Cortex | eQTL | 0.009 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Cortex | sQTL | 0.974 (susie) | 0.324 (0.059) | 4.88e-08 | 0.294 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Frontal_Cortex_BA9 | eQTL | 0.040 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Frontal_Cortex_BA9 | sQTL | 0.956 (susie) | 0.275 (0.049) | 2.22e-08 | 0.101 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Hippocampus | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Hippocampus | sQTL | 0.954 (susie) | 0.307 (0.056) | 5.63e-08 | 0.206 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Hypothalamus | eQTL | 0.069 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Hypothalamus | sQTL | 0.970 (susie) | 0.259 (0.046) | 2.46e-08 | 0.500 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.045 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.956 (susie) | 0.743 (0.166) | 7.45e-06 | 0.094 (17) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.044 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.962 (susie) | 0.279 (0.052) | 7.17e-08 | 0.379 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Substantia_nigra | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Substantia_nigra | sQTL | 0.971 (susie) | 0.303 (0.059) | 2.78e-07 | 0.396 (20) | `coloc_and_smr_heidi_not_rejected` |
| AZI2 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| AZI2 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.803 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CCDC62 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.081 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CCDC62 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.946 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CDHR3 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.009 (abf) | -0.031 (0.060) | 0.603 | 0.844 (3) | `neither_tested_null` |
| CDHR3 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.995 (abf) | 0.164 (0.041) | 5.52e-05 | — (—) | `coloc_and_smr_heidi_untestable` |
| CDHR3 | pd | Brain_Cerebellum | eQTL | 0.008 (abf) | -0.022 (0.042) | 0.602 | 0.858 (12) | `neither_tested_null` |
| CDHR3 | pd | Brain_Cerebellum | sQTL | 0.922 (abf) | 0.181 (0.047) | 0.000126 | — (—) | `coloc_and_smr_heidi_untestable` |
| CTSB | pd | Brain_Amygdala | eQTL | 0.709 (susie) | -0.247 (0.080) | 0.002 | 0.071 (20) | `smr_without_coloc` |
| CTSB | pd | Brain_Amygdala | sQTL | 0.868 (susie) | 0.078 (0.021) | 0.000254 | 0.193 (20) | `coloc_and_smr_heidi_not_rejected` |
| DDRGK1 | pd | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.902 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DDRGK1 | pd | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.884 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NCOR1 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.032 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NCOR1 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.896 (susie) | -0.169 (0.043) | 8.14e-05 | 0.184 (20) | `coloc_and_smr_heidi_not_rejected` |
| NCOR1 | pd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.772 (susie) | 0.368 (0.096) | 0.000128 | 0.154 (16) | `smr_without_coloc` |
| NCOR1 | pd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.805 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SH3GL2 | pd | Brain_Amygdala | eQTL | 0.027 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SH3GL2 | pd | Brain_Amygdala | sQTL | 0.834 (abf) | -0.165 (0.040) | 3.66e-05 | 0.555 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.653 (abf) | 0.587 (0.193) | 0.002 | 0.264 (20) | `smr_without_coloc` |
| SH3GL2 | pd | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.834 (abf) | -0.114 (0.026) | 1.53e-05 | 0.315 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Caudate_basal_ganglia | eQTL | 0.005 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SH3GL2 | pd | Brain_Caudate_basal_ganglia | sQTL | 0.834 (abf) | -0.157 (0.037) | 2.52e-05 | 0.709 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Cortex | eQTL | 0.920 (abf) | 0.510 (0.144) | 0.000401 | 0.495 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Cortex | sQTL | 0.834 (abf) | -0.091 (0.020) | 6.44e-06 | 0.284 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Frontal_Cortex_BA9 | eQTL | 0.963 (abf) | 0.537 (0.128) | 2.62e-05 | 0.489 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Frontal_Cortex_BA9 | sQTL | 0.834 (abf) | -0.135 (0.030) | 6.81e-06 | 0.284 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Hypothalamus | eQTL | 0.017 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SH3GL2 | pd | Brain_Hypothalamus | sQTL | 0.834 (abf) | -0.125 (0.029) | 1.4e-05 | 0.386 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.005 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SH3GL2 | pd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.834 (abf) | -0.140 (0.033) | 1.63e-05 | 0.225 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.068 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SH3GL2 | pd | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.834 (abf) | -0.110 (0.026) | 2.2e-05 | 0.531 (20) | `coloc_and_smr_heidi_not_rejected` |
| SH3GL2 | pd | Brain_Substantia_nigra | eQTL | 0.020 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SH3GL2 | pd | Brain_Substantia_nigra | sQTL | 0.834 (abf) | -0.100 (0.023) | 1.59e-05 | 0.426 (20) | `coloc_and_smr_heidi_not_rejected` |
| TTC19 | pd | Brain_Caudate_basal_ganglia | eQTL | 0.711 (susie) | -0.418 (0.099) | 2.55e-05 | 0.048 (20) | `smr_without_coloc` |
| TTC19 | pd | Brain_Caudate_basal_ganglia | sQTL | 0.804 (susie) | -0.247 (0.065) | 0.000135 | 0.046 (20) | `coloc_and_smr_heidi_not_rejected` |
| TTC19 | pd | Brain_Putamen_basal_ganglia | eQTL | 0.726 (susie) | -0.346 (0.085) | 4.74e-05 | 0.457 (20) | `smr_without_coloc` |
| TTC19 | pd | Brain_Putamen_basal_ganglia | sQTL | 0.899 (susie) | -0.185 (0.045) | 4.13e-05 | 0.768 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Amygdala | eQTL | 0.805 (susie) | 0.153 (0.032) | 1.73e-06 | 0.021 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Amygdala | sQTL | 0.816 (susie) | -0.170 (0.040) | 1.88e-05 | 0.178 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.806 (susie) | 0.106 (0.021) | 7.3e-07 | 0.053 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.822 (susie) | -0.166 (0.038) | 1.28e-05 | 0.110 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Caudate_basal_ganglia | eQTL | 0.809 (susie) | 0.129 (0.026) | 4.83e-07 | 0.070 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Caudate_basal_ganglia | sQTL | 0.902 (susie) | -0.226 (0.059) | 0.000117 | 0.667 (16) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.820 (abf) | -0.180 (0.037) | 1.46e-06 | 0.067 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.827 (susie) | -0.165 (0.042) | 7.31e-05 | 0.176 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Cortex | eQTL | 0.806 (susie) | 0.106 (0.021) | 3.36e-07 | 0.090 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Cortex | sQTL | 0.885 (susie) | 0.190 (0.049) | 0.00012 | 0.129 (19) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Frontal_Cortex_BA9 | eQTL | 0.781 (susie) | 0.129 (0.026) | 6.03e-07 | 0.096 (20) | `smr_without_coloc` |
| ZSWIM7 | pd | Brain_Frontal_Cortex_BA9 | sQTL | 0.846 (susie) | -0.184 (0.042) | 1.38e-05 | 0.332 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Hippocampus | eQTL | 0.849 (susie) | 0.209 (0.045) | 4.37e-06 | 0.194 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Hippocampus | sQTL | 0.844 (susie) | -0.222 (0.059) | 0.000157 | 0.303 (19) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Hypothalamus | eQTL | 0.800 (susie) | 0.160 (0.033) | 1.03e-06 | 0.139 (20) | `smr_without_coloc` |
| ZSWIM7 | pd | Brain_Hypothalamus | sQTL | 0.803 (susie) | -0.160 (0.037) | 1.35e-05 | 0.141 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.780 (susie) | 0.114 (0.023) | 6.8e-07 | 0.041 (20) | `smr_without_coloc` |
| ZSWIM7 | pd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.906 (susie) | -0.166 (0.040) | 3.08e-05 | 0.104 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Putamen_basal_ganglia | eQTL | 0.805 (susie) | 0.164 (0.034) | 1.13e-06 | 0.061 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Putamen_basal_ganglia | sQTL | 0.904 (susie) | -0.220 (0.050) | 1.01e-05 | 0.024 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.752 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZSWIM7 | pd | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.823 (susie) | 0.161 (0.043) | 0.000167 | 0.180 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Substantia_nigra | eQTL | 0.839 (susie) | 0.193 (0.043) | 5.72e-06 | 0.265 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZSWIM7 | pd | Brain_Substantia_nigra | sQTL | 0.878 (susie) | -0.192 (0.047) | 3.79e-05 | 0.010 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Amygdala | eQTL | 0.032 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ACTR1B | scz | Brain_Amygdala | eQTL | 0.032 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ACTR1B | scz | Brain_Amygdala | sQTL | 0.811 (abf) | -0.084 (0.022) | 9.9e-05 | 0.000929 (12) | `coloc_and_smr_heidi_rejected` |
| ACTR1B | scz | Brain_Amygdala | sQTL | 0.811 (abf) | -0.084 (0.022) | 9.9e-05 | 0.000929 (12) | `coloc_and_smr_heidi_rejected` |
| ACTR1B | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.886 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.886 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.848 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.848 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Caudate_basal_ganglia | eQTL | 0.994 (abf) | -0.262 (0.067) | 8.35e-05 | 0.112 (10) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Caudate_basal_ganglia | eQTL | 0.994 (abf) | -0.262 (0.067) | 8.35e-05 | 0.112 (10) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Caudate_basal_ganglia | sQTL | 0.997 (abf) | -0.093 (0.021) | 6.39e-06 | 0.610 (11) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Caudate_basal_ganglia | sQTL | 0.997 (abf) | -0.093 (0.021) | 6.39e-06 | 0.610 (11) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.995 (abf) | -0.088 (0.017) | 1.47e-07 | 0.028 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.995 (abf) | -0.088 (0.017) | 1.47e-07 | 0.028 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.995 (abf) | -0.057 (0.011) | 2.25e-07 | 0.075 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.995 (abf) | -0.057 (0.011) | 2.25e-07 | 0.075 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellum | eQTL | 0.995 (abf) | -0.074 (0.014) | 1.12e-07 | 0.025 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellum | eQTL | 0.995 (abf) | -0.074 (0.014) | 1.12e-07 | 0.025 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellum | sQTL | 0.995 (abf) | -0.050 (0.010) | 2.04e-07 | 0.044 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cerebellum | sQTL | 0.995 (abf) | -0.050 (0.010) | 2.04e-07 | 0.044 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cortex | eQTL | 0.995 (abf) | -0.231 (0.053) | 1.27e-05 | 0.041 (9) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cortex | eQTL | 0.995 (abf) | -0.231 (0.053) | 1.27e-05 | 0.041 (9) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cortex | sQTL | 0.995 (abf) | -0.063 (0.013) | 1.02e-06 | 0.045 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Cortex | sQTL | 0.995 (abf) | -0.063 (0.013) | 1.02e-06 | 0.045 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.931 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.931 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.996 (abf) | -0.088 (0.020) | 9.19e-06 | 0.026 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.996 (abf) | -0.088 (0.020) | 9.19e-06 | 0.026 (20) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Hippocampus | eQTL | 0.962 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Hippocampus | eQTL | 0.962 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Hippocampus | sQTL | 0.978 (abf) | -0.088 (0.022) | 4.83e-05 | 0.163 (7) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Hippocampus | sQTL | 0.978 (abf) | -0.088 (0.022) | 4.83e-05 | 0.163 (7) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Hypothalamus | eQTL | 0.996 (abf) | -0.142 (0.032) | 6.41e-06 | 0.037 (10) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Hypothalamus | eQTL | 0.996 (abf) | -0.142 (0.032) | 6.41e-06 | 0.037 (10) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Hypothalamus | sQTL | 0.996 (abf) | -0.083 (0.019) | 7.06e-06 | 0.243 (12) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Hypothalamus | sQTL | 0.996 (abf) | -0.083 (0.019) | 7.06e-06 | 0.243 (12) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.958 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.958 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.995 (abf) | -0.078 (0.017) | 3.17e-06 | 0.053 (12) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.995 (abf) | -0.078 (0.017) | 3.17e-06 | 0.053 (12) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.897 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.897 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ACTR1B | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.989 (abf) | -0.072 (0.017) | 2.19e-05 | 0.100 (10) | `coloc_and_smr_heidi_not_rejected` |
| ACTR1B | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.989 (abf) | -0.072 (0.017) | 2.19e-05 | 0.100 (10) | `coloc_and_smr_heidi_not_rejected` |
| ASB3 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.150 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ASB3 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.859 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CATSPER2 | scz | Brain_Hippocampus | eQTL | 0.039 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CATSPER2 | scz | Brain_Hippocampus | eQTL | 0.039 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CATSPER2 | scz | Brain_Hippocampus | sQTL | 0.894 (abf) | -0.064 (0.013) | 3.6e-07 | 0.322 (20) | `coloc_and_smr_heidi_not_rejected` |
| CATSPER2 | scz | Brain_Hippocampus | sQTL | 0.894 (abf) | -0.064 (0.013) | 3.6e-07 | 0.322 (20) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Hippocampus | eQTL | 0.933 (abf) | -0.131 (0.035) | 0.000223 | 0.498 (20) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Hippocampus | sQTL | 0.913 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CCDC122 | scz | Brain_Hypothalamus | eQTL | 0.880 (abf) | -0.120 (0.031) | 0.0001 | 0.053 (20) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Hypothalamus | sQTL | 0.856 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CCDC122 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.925 (abf) | -0.127 (0.034) | 0.000161 | 0.061 (20) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.932 (abf) | 0.072 (0.019) | 0.000187 | 0.553 (16) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.933 (abf) | -0.096 (0.026) | 0.000243 | 0.564 (20) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.944 (abf) | 0.053 (0.013) | 5.29e-05 | 0.313 (20) | `coloc_and_smr_heidi_not_rejected` |
| CCDC122 | scz | Brain_Substantia_nigra | eQTL | 0.684 (abf) | -0.048 (0.021) | 0.025 | 0.052 (18) | `neither_tested_null` |
| CCDC122 | scz | Brain_Substantia_nigra | sQTL | 0.834 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CCS | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.008 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CCS | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.937 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CD46 | scz | Brain_Cerebellum | eQTL | 0.783 (abf) | 0.056 (0.014) | 7.91e-05 | 0.428 (20) | `smr_without_coloc` |
| CD46 | scz | Brain_Cerebellum | sQTL | 0.804 (abf) | -0.059 (0.016) | 0.000155 | 0.104 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Amygdala | eQTL | 0.624 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Amygdala | sQTL | 0.929 (susie) | 0.066 (0.017) | 0.000166 | 0.726 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.713 (susie) | 0.192 (0.055) | 0.000455 | 0.088 (20) | `smr_without_coloc` |
| CDIP1 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.951 (susie) | -0.099 (0.025) | 7.18e-05 | 0.563 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Cortex | eQTL | 0.802 (susie) | 0.108 (0.027) | 6.92e-05 | 0.548 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Cortex | sQTL | 0.832 (susie) | 0.065 (0.016) | 4.86e-05 | 0.448 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Hippocampus | eQTL | 0.021 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Hippocampus | sQTL | 0.801 (susie) | 0.085 (0.026) | 0.000961 | 0.378 (20) | `coloc_smr_not_significant` |
| CDIP1 | scz | Brain_Hypothalamus | eQTL | 0.037 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Hypothalamus | sQTL | 0.935 (susie) | -0.089 (0.023) | 0.000132 | 0.893 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Substantia_nigra | eQTL | 0.562 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Substantia_nigra | sQTL | 0.928 (abf) | 0.083 (0.022) | 0.000138 | 0.683 (20) | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | Brain_Cerebellum | eQTL | 0.822 (abf) | -0.095 (0.025) | 0.000127 | 0.587 (20) | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | Brain_Cerebellum | sQTL | 0.830 (abf) | -0.114 (0.035) | 0.001 | 0.237 (20) | `coloc_smr_not_significant` |
| CRELD2 | scz | Brain_Amygdala | eQTL | 0.413 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Amygdala | sQTL | 0.918 (susie) | -0.084 (0.021) | 5.73e-05 | 0.414 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.217 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.934 (susie) | 0.049 (0.010) | 3.93e-07 | 0.338 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.713 (susie) | 0.190 (0.052) | 0.000282 | 0.374 (20) | `smr_without_coloc` |
| CRELD2 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.933 (susie) | -0.052 (0.011) | 1e-06 | 0.757 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.017 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.943 (susie) | -0.052 (0.010) | 4.87e-07 | 0.370 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Cerebellum | eQTL | 0.050 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Cerebellum | sQTL | 0.949 (susie) | 0.046 (0.009) | 2.3e-07 | 0.509 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Cortex | eQTL | 0.122 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Cortex | sQTL | 0.947 (susie) | 0.048 (0.009) | 3.72e-07 | 0.753 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.943 (susie) | -0.052 (0.010) | 6.1e-07 | 0.666 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Hippocampus | eQTL | 0.763 (susie) | 0.141 (0.036) | 9.07e-05 | 0.323 (20) | `smr_without_coloc` |
| CRELD2 | scz | Brain_Hippocampus | sQTL | 0.919 (susie) | -0.064 (0.015) | 2.45e-05 | 0.198 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Hypothalamus | eQTL | 0.852 (susie) | 0.185 (0.045) | 3.7e-05 | 0.619 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Hypothalamus | sQTL | 0.932 (susie) | 0.043 (0.009) | 4.57e-07 | 0.671 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.629 (susie) | 0.136 (0.033) | 3.76e-05 | 0.147 (20) | `smr_without_coloc` |
| CRELD2 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.939 (susie) | -0.057 (0.011) | 8.29e-07 | 0.555 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.380 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.945 (susie) | 0.048 (0.009) | 4.59e-07 | 0.616 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.825 (susie) | 0.139 (0.038) | 0.000242 | 0.229 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.919 (susie) | 0.055 (0.011) | 9.47e-07 | 0.380 (20) | `coloc_and_smr_heidi_not_rejected` |
| CRELD2 | scz | Brain_Substantia_nigra | eQTL | 0.442 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CRELD2 | scz | Brain_Substantia_nigra | sQTL | 0.927 (susie) | 0.045 (0.009) | 7.41e-07 | 0.612 (20) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Amygdala | eQTL | 0.072 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Amygdala | sQTL | 0.878 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.849 (abf) | -0.139 (0.027) | 3.63e-07 | 0.034 (20) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Caudate_basal_ganglia | eQTL | 0.495 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Caudate_basal_ganglia | sQTL | 0.917 (abf) | -0.194 (0.043) | 6e-06 | 0.268 (15) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.858 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.924 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Cortex | eQTL | 0.038 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Cortex | sQTL | 0.939 (abf) | -0.145 (0.028) | 1.75e-07 | 0.356 (20) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.949 (abf) | -0.151 (0.028) | 7.34e-08 | 0.004 (20) | `coloc_and_smr_heidi_rejected` |
| DGKZ | scz | Brain_Hippocampus | eQTL | 0.035 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Hippocampus | sQTL | 0.838 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.132 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.928 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.035 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.892 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DNAJA3 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DNAJA3 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.963 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| DNAJA3 | scz | Brain_Cerebellum | eQTL | 0.282 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DNAJA3 | scz | Brain_Cerebellum | sQTL | 0.902 (susie) | -0.076 (0.022) | 0.000428 | 0.457 (20) | `coloc_smr_not_significant` |
| DNAJA3 | scz | Brain_Hypothalamus | eQTL | 0.028 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DNAJA3 | scz | Brain_Hypothalamus | sQTL | 0.900 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DOC2A | scz | Brain_Amygdala | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Amygdala | sQTL | 0.947 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DOC2A | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.009 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.959 (abf) | -0.140 (0.028) | 9.29e-07 | 0.061 (20) | `coloc_and_smr_heidi_not_rejected` |
| DOC2A | scz | Brain_Caudate_basal_ganglia | eQTL | 0.671 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Caudate_basal_ganglia | sQTL | 0.871 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DOC2A | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.034 (susie) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.969 (susie) | -0.141 (0.029) | 1.06e-06 | 0.704 (17) | `coloc_and_smr_heidi_not_rejected` |
| DOC2A | scz | Brain_Cortex | eQTL | 0.063 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Cortex | sQTL | 0.967 (susie) | 0.087 (0.014) | 5.39e-10 | 0.007 (20) | `coloc_and_smr_heidi_rejected` |
| DOC2A | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.970 (susie) | 0.143 (0.024) | 3.8e-09 | 0.006 (20) | `coloc_and_smr_heidi_rejected` |
| DOC2A | scz | Brain_Hippocampus | eQTL | 0.105 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Hippocampus | sQTL | 0.964 (abf) | -0.152 (0.032) | 2.6e-06 | 0.491 (19) | `coloc_and_smr_heidi_not_rejected` |
| DOC2A | scz | Brain_Hypothalamus | eQTL | 0.059 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Hypothalamus | sQTL | 0.949 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DOC2A | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.755 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.901 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| EFHB | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.000984 (abf) | 0.000261 (0.008) | 0.973 | 0.002 (20) | `neither_tested_null` |
| EFHB | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.965 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| FAM120AOS | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.000146 (abf) | -0.000713 (0.031) | 0.981 | 0.000133 (20) | `neither_tested_null` |
| FAM120AOS | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.889 (abf) | 0.082 (0.021) | 7.84e-05 | 0.109 (20) | `coloc_and_smr_heidi_not_rejected` |
| FAM120AOS | scz | Brain_Cerebellum | eQTL | 0.113 (abf) | -0.005 (0.027) | 0.862 | 0.000313 (20) | `neither_tested_null` |
| FAM120AOS | scz | Brain_Cerebellum | sQTL | 0.919 (abf) | -0.073 (0.020) | 0.000208 | 0.152 (20) | `coloc_and_smr_heidi_not_rejected` |
| FAM184A | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.057 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FAM184A | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.832 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| FANCI | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.166 (abf) | -0.098 (0.032) | 0.002 | 0.177 (20) | `neither_tested_null` |
| FANCI | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.879 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| FGFR1 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FGFR1 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.889 (abf) | 0.092 (0.023) | 6.81e-05 | 0.998 (14) | `coloc_and_smr_heidi_not_rejected` |
| FGFR1 | scz | Brain_Cerebellum | eQTL | 0.653 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FGFR1 | scz | Brain_Cerebellum | sQTL | 0.841 (abf) | -0.058 (0.012) | 3.18e-06 | 0.838 (20) | `coloc_and_smr_heidi_not_rejected` |
| FOXN2 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.997 (abf) | -0.063 (0.013) | 5.01e-07 | 0.859 (20) | `coloc_and_smr_heidi_not_rejected` |
| FOXN2 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.997 (abf) | 0.073 (0.017) | 1.97e-05 | 0.274 (20) | `coloc_and_smr_heidi_not_rejected` |
| FOXN2 | scz | Brain_Cerebellum | eQTL | 0.997 (abf) | -0.054 (0.011) | 5.12e-07 | 0.631 (20) | `coloc_and_smr_heidi_not_rejected` |
| FOXN2 | scz | Brain_Cerebellum | sQTL | 0.996 (abf) | 0.056 (0.012) | 3.14e-06 | 0.731 (20) | `coloc_and_smr_heidi_not_rejected` |
| GABBR2 | scz | Brain_Cerebellum | eQTL | 0.014 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GABBR2 | scz | Brain_Cerebellum | sQTL | 0.975 (susie) | -0.100 (0.026) | 9.97e-05 | 0.604 (5) | `coloc_and_smr_heidi_not_rejected` |
| GALNT15 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.008 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GALNT15 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.850 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GLYCTK | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.499 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GLYCTK | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.499 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GLYCTK | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.812 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GLYCTK | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.812 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GPM6A | scz | Brain_Amygdala | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Amygdala | sQTL | 0.841 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GPM6A | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.061 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.992 (susie) | 0.188 (0.043) | 1.18e-05 | 0.004 (13) | `coloc_and_smr_heidi_rejected` |
| GPM6A | scz | Brain_Caudate_basal_ganglia | eQTL | 0.044 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Caudate_basal_ganglia | sQTL | 0.993 (susie) | 0.167 (0.037) | 4.97e-06 | 0.036 (8) | `coloc_and_smr_heidi_not_rejected` |
| GPM6A | scz | Brain_Putamen_basal_ganglia | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Putamen_basal_ganglia | sQTL | 0.908 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GPR135 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.011 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPR135 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.933 (susie) | 0.031 (0.015) | 0.039 | 0.016 (20) | `coloc_smr_not_significant` |
| GPR135 | scz | Brain_Cerebellum | eQTL | 0.004 (susie) | 0.073 (0.025) | 0.003 | 0.006 (20) | `neither_tested_null` |
| GPR135 | scz | Brain_Cerebellum | sQTL | 0.933 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| HMOX2 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| HMOX2 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.829 (susie) | -0.089 (0.025) | 0.000371 | 0.138 (20) | `coloc_smr_not_significant` |
| IDH3B | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.676 (abf) | -0.178 (0.056) | 0.002 | 0.748 (15) | `neither_tested_null` |
| IDH3B | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.676 (abf) | -0.178 (0.056) | 0.002 | 0.748 (15) | `smr_without_coloc` |
| IDH3B | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.840 (abf) | 0.042 (0.011) | 0.000169 | 0.644 (18) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.840 (abf) | 0.042 (0.011) | 0.000169 | 0.644 (18) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Cerebellum | eQTL | 0.686 (abf) | -0.169 (0.056) | 0.003 | 0.882 (12) | `neither_tested_null` |
| IDH3B | scz | Brain_Cerebellum | eQTL | 0.686 (abf) | -0.169 (0.056) | 0.003 | 0.882 (12) | `smr_without_coloc` |
| IDH3B | scz | Brain_Cerebellum | sQTL | 0.860 (abf) | 0.041 (0.011) | 0.000145 | 0.659 (20) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Cerebellum | sQTL | 0.860 (abf) | 0.041 (0.011) | 0.000145 | 0.659 (20) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Cortex | eQTL | 0.811 (abf) | -0.214 (0.063) | 0.000695 | 0.636 (13) | `coloc_smr_not_significant` |
| IDH3B | scz | Brain_Cortex | eQTL | 0.811 (abf) | -0.214 (0.063) | 0.000695 | 0.636 (13) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Cortex | sQTL | 0.859 (abf) | 0.047 (0.013) | 0.000211 | 0.783 (15) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Cortex | sQTL | 0.859 (abf) | 0.047 (0.013) | 0.000211 | 0.783 (15) | `coloc_and_smr_heidi_not_rejected` |
| IDH3B | scz | Brain_Hippocampus | eQTL | 0.662 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IDH3B | scz | Brain_Hippocampus | eQTL | 0.662 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IDH3B | scz | Brain_Hippocampus | sQTL | 0.826 (abf) | 0.064 (0.018) | 0.000439 | 0.926 (13) | `coloc_smr_not_significant` |
| IDH3B | scz | Brain_Hippocampus | sQTL | 0.826 (abf) | 0.064 (0.018) | 0.000439 | 0.926 (13) | `coloc_and_smr_heidi_not_rejected` |
| IKBIP | scz | Brain_Frontal_Cortex_BA9 | eQTL | 1.27e-05 (susie) | -0.002 (0.017) | 0.900 | 0.345 (20) | `neither_tested_null` |
| IKBIP | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.853 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.967 (susie) | 0.165 (0.027) | 1.43e-09 | 0.009 (20) | `coloc_and_smr_heidi_rejected` |
| INO80E | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.958 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | scz | Brain_Cerebellum | eQTL | 0.970 (susie) | 0.149 (0.024) | 6.71e-10 | 0.006 (20) | `coloc_and_smr_heidi_rejected` |
| INO80E | scz | Brain_Cerebellum | sQTL | 0.953 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | scz | Brain_Hippocampus | eQTL | 0.493 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| INO80E | scz | Brain_Hippocampus | sQTL | 0.921 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.964 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.936 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| INO80E | scz | Brain_Putamen_basal_ganglia | eQTL | 0.962 (susie) | 0.207 (0.043) | 1.41e-06 | 0.010 (20) | `coloc_and_smr_heidi_not_rejected` |
| INO80E | scz | Brain_Putamen_basal_ganglia | sQTL | 0.943 (abf) | 0.027 (0.017) | 0.121 | 0.000387 (20) | `coloc_smr_not_significant` |
| IRF3 | scz | Brain_Amygdala | eQTL | 0.274 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Amygdala | sQTL | 0.971 (abf) | 0.077 (0.017) | 6.72e-06 | 0.096 (15) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.144 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.992 (abf) | 0.094 (0.020) | 1.57e-06 | 0.102 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.403 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.988 (abf) | 0.075 (0.013) | 2.5e-08 | 0.241 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Cortex | eQTL | 0.244 (abf) | -0.196 (0.060) | 0.001 | 0.002 (7) | `neither_tested_null` |
| IRF3 | scz | Brain_Cortex | sQTL | 0.994 (abf) | 0.078 (0.015) | 1.17e-07 | 0.070 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.074 (abf) | -0.151 (0.044) | 0.000522 | 0.010 (16) | `smr_without_coloc` |
| IRF3 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.983 (abf) | 0.087 (0.018) | 2.12e-06 | 0.166 (18) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Hippocampus | eQTL | 0.066 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Hippocampus | sQTL | 0.976 (abf) | 0.090 (0.018) | 6.69e-07 | 0.158 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Hypothalamus | eQTL | 0.148 (abf) | -0.142 (0.040) | 0.000376 | 0.004 (20) | `smr_without_coloc` |
| IRF3 | scz | Brain_Hypothalamus | sQTL | 0.990 (abf) | 0.089 (0.018) | 9.49e-07 | 0.248 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.561 (abf) | -0.204 (0.055) | 0.000205 | 0.071 (20) | `smr_without_coloc` |
| IRF3 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.990 (abf) | 0.094 (0.019) | 1.35e-06 | 0.072 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.614 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.991 (abf) | 0.084 (0.016) | 1.17e-07 | 0.193 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.064 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.991 (abf) | 0.090 (0.019) | 1.77e-06 | 0.103 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Substantia_nigra | eQTL | 0.074 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Substantia_nigra | sQTL | 0.990 (abf) | 0.082 (0.016) | 2.91e-07 | 0.179 (20) | `coloc_and_smr_heidi_not_rejected` |
| KLC1 | scz | Brain_Cortex | eQTL | 4.35e-05 (abf) | 0.150 (0.097) | 0.123 | 0.000978 (20) | `neither_tested_null` |
| KLC1 | scz | Brain_Cortex | sQTL | 0.855 (abf) | -0.293 (0.064) | 4.16e-06 | 0.582 (20) | `coloc_and_smr_heidi_not_rejected` |
| L3HYPDH | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.003 (abf) | 0.046 (0.015) | 0.003 | 0.036 (20) | `neither_tested_null` |
| L3HYPDH | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.944 (susie) | 0.110 (0.029) | 0.000143 | 0.092 (19) | `coloc_and_smr_heidi_not_rejected` |
| L3HYPDH | scz | Brain_Cerebellum | eQTL | 0.003 (susie) | 0.040 (0.013) | 0.003 | 0.055 (20) | `neither_tested_null` |
| L3HYPDH | scz | Brain_Cerebellum | sQTL | 0.937 (susie) | -0.032 (0.016) | 0.040 | 0.021 (20) | `coloc_smr_not_significant` |
| LPCAT4 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.007 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| LPCAT4 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.007 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| LPCAT4 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.884 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| LPCAT4 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.884 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| MAD1L1 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.026 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MAD1L1 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.992 (abf) | -0.152 (0.030) | 2.72e-07 | 0.077 (20) | `coloc_and_smr_heidi_not_rejected` |
| MAD1L1 | scz | Brain_Cortex | eQTL | 0.112 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MAD1L1 | scz | Brain_Cortex | sQTL | 0.933 (abf) | -0.141 (0.027) | 1.5e-07 | 0.083 (20) | `coloc_and_smr_heidi_not_rejected` |
| MAD1L1 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.059 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MAD1L1 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.966 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| MAP2K5 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.012 (abf) | 0.087 (0.042) | 0.041 | 0.004 (20) | `neither_tested_null` |
| MAP2K5 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.822 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| MAP2K5 | scz | Brain_Cerebellum | eQTL | 0.007 (abf) | 0.045 (0.024) | 0.060 | 0.507 (20) | `neither_tested_null` |
| MAP2K5 | scz | Brain_Cerebellum | sQTL | 0.852 (abf) | -0.074 (0.022) | 0.000754 | 0.054 (20) | `coloc_smr_not_significant` |
| MAP2K5 | scz | Brain_Cortex | eQTL | 0.006 (abf) | 0.084 (0.049) | 0.085 | 0.088 (19) | `neither_tested_null` |
| MAP2K5 | scz | Brain_Cortex | sQTL | 0.956 (abf) | -0.086 (0.025) | 0.000514 | 0.700 (20) | `coloc_smr_not_significant` |
| MAP7D1 | scz | Brain_Hippocampus | eQTL | 0.061 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MAP7D1 | scz | Brain_Hippocampus | sQTL | 0.860 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| MAP7D1 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.129 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MAP7D1 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.961 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| MAP7D1 | scz | Brain_Substantia_nigra | eQTL | 0.481 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MAP7D1 | scz | Brain_Substantia_nigra | sQTL | 0.911 (abf) | 0.080 (0.019) | 3.89e-05 | 0.004 (20) | `coloc_and_smr_heidi_rejected` |
| MED19 | scz | Brain_Hippocampus | eQTL | 0.036 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MED19 | scz | Brain_Hippocampus | sQTL | 0.810 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| MRPS33 | scz | Brain_Cerebellum | eQTL | 0.017 (abf) | 0.098 (0.041) | 0.018 | 0.465 (16) | `neither_tested_null` |
| MRPS33 | scz | Brain_Cerebellum | sQTL | 0.859 (abf) | -0.055 (0.021) | 0.007 | 0.210 (20) | `coloc_smr_not_significant` |
| MRPS33 | scz | Brain_Hypothalamus | eQTL | 0.021 (abf) | 0.128 (0.056) | 0.022 | 0.282 (20) | `neither_tested_null` |
| MRPS33 | scz | Brain_Hypothalamus | sQTL | 0.994 (abf) | -0.082 (0.020) | 3.04e-05 | 0.096 (20) | `coloc_and_smr_heidi_not_rejected` |
| MRPS33 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.580 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| MRPS33 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.848 (abf) | -0.044 (0.016) | 0.005 | 0.110 (20) | `coloc_smr_not_significant` |
| NDUFAF7 | scz | Brain_Amygdala | eQTL | 0.113 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NDUFAF7 | scz | Brain_Amygdala | sQTL | 0.859 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NDUFAF7 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.890 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NDUFAF7 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.925 (abf) | -0.055 (0.012) | 8.4e-06 | 0.171 (20) | `coloc_and_smr_heidi_not_rejected` |
| NEK4 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.341 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NEK4 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.820 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NMRAL1 | scz | Brain_Amygdala | eQTL | 0.595 (susie) | 0.089 (0.024) | 0.000187 | 0.310 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Amygdala | sQTL | 0.915 (susie) | -0.048 (0.012) | 7.84e-05 | 0.287 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.652 (susie) | 0.078 (0.021) | 0.000168 | 0.068 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.912 (susie) | 0.055 (0.014) | 8.58e-05 | 0.308 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.609 (susie) | 0.036 (0.009) | 7.42e-05 | 0.036 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.899 (susie) | -0.062 (0.016) | 0.000119 | 0.396 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Cerebellum | eQTL | 0.760 (susie) | 0.038 (0.009) | 1.64e-05 | 0.378 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Cerebellum | sQTL | 0.902 (susie) | -0.053 (0.014) | 0.000145 | 0.223 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Cortex | eQTL | 0.584 (susie) | 0.084 (0.022) | 0.000145 | 0.053 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Cortex | sQTL | 0.932 (susie) | -0.041 (0.010) | 2.46e-05 | 0.271 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.539 (susie) | 0.091 (0.024) | 0.000145 | 0.017 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.914 (susie) | 0.049 (0.012) | 2.94e-05 | 0.253 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Hippocampus | eQTL | 0.737 (susie) | 0.076 (0.019) | 6.53e-05 | 0.143 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Hippocampus | sQTL | 0.906 (susie) | -0.040 (0.010) | 5.4e-05 | 0.158 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Hypothalamus | eQTL | 0.640 (susie) | 0.093 (0.024) | 0.000137 | 0.352 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Hypothalamus | sQTL | 0.909 (susie) | 0.082 (0.023) | 0.000486 | 0.392 (20) | `coloc_smr_not_significant` |
| NMRAL1 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.628 (susie) | -0.036 (0.034) | 0.292 | 0.150 (20) | `neither_tested_null` |
| NMRAL1 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.904 (susie) | 0.068 (0.018) | 0.000178 | 0.243 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.692 (susie) | 0.093 (0.025) | 0.000177 | 0.074 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.931 (susie) | -0.047 (0.011) | 1.72e-05 | 0.346 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.604 (susie) | 0.074 (0.019) | 0.000124 | 0.040 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.925 (susie) | 0.041 (0.009) | 1.59e-05 | 0.390 (20) | `coloc_and_smr_heidi_not_rejected` |
| NMRAL1 | scz | Brain_Substantia_nigra | eQTL | 0.722 (susie) | 0.073 (0.018) | 6.46e-05 | 0.111 (20) | `smr_without_coloc` |
| NMRAL1 | scz | Brain_Substantia_nigra | sQTL | 0.859 (susie) | -0.037 (0.009) | 8.36e-05 | 0.085 (20) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Amygdala | eQTL | 0.225 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Amygdala | sQTL | 0.823 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NT5C2 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.243 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.901 (susie) | 0.155 (0.031) | 6.91e-07 | 0.243 (12) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Cerebellar_Hemisphere | eQTL | 1.84e-08 (susie) | 0.107 (0.051) | 0.035 | 0.001 (20) | `neither_tested_null` |
| NT5C2 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.944 (susie) | 0.113 (0.019) | 5.87e-09 | 0.645 (20) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Cerebellum | eQTL | 1.8e-08 (susie) | 0.079 (0.035) | 0.024 | 0.005 (20) | `neither_tested_null` |
| NT5C2 | scz | Brain_Cerebellum | sQTL | 0.959 (susie) | 0.088 (0.014) | 8.8e-10 | 0.774 (20) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Cortex | eQTL | 0.963 (susie) | 0.029 (0.030) | 0.328 | 0.097 (20) | `coloc_smr_not_significant` |
| NT5C2 | scz | Brain_Cortex | sQTL | 0.907 (susie) | 0.134 (0.029) | 3.03e-06 | 0.552 (17) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.004 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.822 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NT5C2 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.406 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.907 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NUCB2 | scz | Brain_Hypothalamus | eQTL | 0.023 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NUCB2 | scz | Brain_Hypothalamus | sQTL | 0.808 (abf) | -0.075 (0.024) | 0.002 | 0.380 (18) | `coloc_smr_not_significant` |
| NUP50 | scz | Brain_Amygdala | eQTL | 0.895 (abf) | -0.091 (0.024) | 0.000166 | 0.678 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Amygdala | sQTL | 0.835 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NUP50 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.925 (abf) | -0.100 (0.024) | 4.51e-05 | 0.978 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.918 (abf) | 0.064 (0.017) | 0.000159 | 0.576 (17) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.972 (abf) | -0.098 (0.023) | 1.64e-05 | 0.746 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.960 (abf) | -0.094 (0.026) | 0.000307 | 0.634 (19) | `coloc_smr_not_significant` |
| NUP50 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.948 (abf) | -0.163 (0.044) | 0.000237 | 0.802 (18) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.945 (abf) | -0.092 (0.025) | 0.000191 | 0.831 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cerebellum | eQTL | 0.760 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NUP50 | scz | Brain_Cerebellum | sQTL | 0.953 (abf) | -0.080 (0.021) | 0.000135 | 0.888 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cortex | eQTL | 0.906 (abf) | -0.106 (0.027) | 8.16e-05 | 0.409 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cortex | sQTL | 0.935 (abf) | 0.066 (0.018) | 0.000345 | 0.813 (20) | `coloc_smr_not_significant` |
| NUP50 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.943 (abf) | -0.127 (0.032) | 6.09e-05 | 0.385 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.927 (abf) | 0.058 (0.015) | 0.000167 | 0.951 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Hippocampus | eQTL | 0.932 (abf) | -0.126 (0.031) | 4.96e-05 | 0.774 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Hippocampus | sQTL | 0.910 (abf) | 0.069 (0.020) | 0.000583 | 0.817 (17) | `coloc_smr_not_significant` |
| NUP50 | scz | Brain_Hypothalamus | eQTL | 0.791 (abf) | -0.115 (0.030) | 0.000156 | 0.806 (20) | `smr_without_coloc` |
| NUP50 | scz | Brain_Hypothalamus | sQTL | 0.918 (abf) | 0.065 (0.018) | 0.000389 | 0.457 (15) | `coloc_smr_not_significant` |
| NUP50 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.839 (abf) | -0.088 (0.022) | 7.71e-05 | 0.660 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.956 (abf) | -0.069 (0.018) | 0.000101 | 0.848 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.733 (abf) | -0.073 (0.020) | 0.000229 | 0.657 (20) | `smr_without_coloc` |
| NUP50 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.805 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PAK6 | scz | Brain_Cerebellum | eQTL | 0.992 (abf) | 0.166 (0.038) | 1.13e-05 | 0.167 (8) | `coloc_and_smr_heidi_not_rejected` |
| PAK6 | scz | Brain_Cerebellum | eQTL | 0.992 (abf) | 0.166 (0.038) | 1.13e-05 | 0.167 (8) | `coloc_and_smr_heidi_not_rejected` |
| PAK6 | scz | Brain_Cerebellum | sQTL | 0.866 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PAK6 | scz | Brain_Cerebellum | sQTL | 0.866 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PAM16 | scz | Brain_Amygdala | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PAM16 | scz | Brain_Amygdala | sQTL | 0.908 (susie) | -0.095 (0.026) | 0.000287 | 0.342 (20) | `coloc_and_smr_heidi_not_rejected` |
| PBRM1 | scz | Brain_Cerebellum | eQTL | 0.523 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PBRM1 | scz | Brain_Cerebellum | sQTL | 0.985 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PCBP3 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.757 (abf) | 0.127 (0.034) | 0.00022 | 0.000985 (20) | `smr_without_coloc` |
| PCBP3 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.819 (abf) | -0.055 (0.014) | 0.000107 | 0.000458 (20) | `coloc_and_smr_heidi_rejected` |
| POLG | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.038 (abf) | -0.157 (0.054) | 0.004 | 0.020 (20) | `neither_tested_null` |
| POLG | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.982 (abf) | -0.048 (0.010) | 4.79e-07 | 0.770 (20) | `coloc_and_smr_heidi_not_rejected` |
| POLG | scz | Brain_Cerebellum | eQTL | 0.047 (abf) | -0.110 (0.034) | 0.001 | 0.009 (20) | `neither_tested_null` |
| POLG | scz | Brain_Cerebellum | sQTL | 0.972 (abf) | -0.049 (0.010) | 6e-07 | 0.709 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Amygdala | eQTL | 0.484 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIL2 | scz | Brain_Amygdala | sQTL | 0.943 (abf) | -0.040 (0.010) | 0.000113 | 0.037 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.299 (abf) | 0.120 (0.047) | 0.011 | 0.146 (20) | `neither_tested_null` |
| PPIL2 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.933 (abf) | -0.049 (0.013) | 0.000116 | 0.011 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.928 (abf) | 0.174 (0.051) | 0.000569 | 0.915 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.938 (abf) | -0.040 (0.010) | 5.36e-05 | 0.012 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Cerebellum | eQTL | 0.938 (abf) | 0.113 (0.030) | 0.00015 | 0.774 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Cerebellum | sQTL | 0.941 (abf) | -0.037 (0.009) | 5.19e-05 | 0.007 (20) | `coloc_and_smr_heidi_rejected` |
| PPIL2 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.299 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIL2 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.939 (abf) | -0.040 (0.010) | 6.53e-05 | 0.008 (20) | `coloc_and_smr_heidi_rejected` |
| PPIL2 | scz | Brain_Hypothalamus | eQTL | 0.159 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIL2 | scz | Brain_Hypothalamus | sQTL | 0.935 (abf) | -0.045 (0.012) | 0.000104 | 0.105 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.926 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PPIL2 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.855 (abf) | -0.040 (0.012) | 0.00052 | 0.032 (20) | `coloc_smr_not_significant` |
| PPIL2 | scz | Brain_Substantia_nigra | eQTL | 0.867 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PPIL2 | scz | Brain_Substantia_nigra | sQTL | 0.936 (abf) | -0.053 (0.014) | 0.000186 | 0.031 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIP5K1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIP5K1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.926 (abf) | 0.114 (0.026) | 1.34e-05 | 0.422 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIP5K1 | scz | Brain_Cerebellum | eQTL | 0.003 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIP5K1 | scz | Brain_Cerebellum | sQTL | 0.869 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PRMT7 | scz | Brain_Amygdala | eQTL | 0.023 (abf) | 0.063 (0.026) | 0.016 | 0.289 (20) | `neither_tested_null` |
| PRMT7 | scz | Brain_Amygdala | sQTL | 0.934 (abf) | 0.035 (0.011) | 0.001 | 0.459 (20) | `coloc_smr_not_significant` |
| PSMD6 | scz | Brain_Hippocampus | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PSMD6 | scz | Brain_Hippocampus | sQTL | 0.976 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RAI1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.717 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RAI1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.946 (susie) | -0.074 (0.017) | 1.07e-05 | 0.439 (20) | `coloc_and_smr_heidi_not_rejected` |
| RAI1 | scz | Brain_Cerebellum | eQTL | 0.308 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RAI1 | scz | Brain_Cerebellum | sQTL | 0.971 (abf) | 0.092 (0.022) | 2.03e-05 | 0.414 (20) | `coloc_and_smr_heidi_not_rejected` |
| RBM6 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.750 (abf) | 0.082 (0.020) | 6.03e-05 | 3.61e-05 (20) | `smr_without_coloc` |
| RBM6 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.869 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RCBTB1 | scz | Brain_Hippocampus | eQTL | 0.971 (susie) | -0.095 (0.022) | 2.16e-05 | 0.229 (20) | `coloc_and_smr_heidi_not_rejected` |
| RCBTB1 | scz | Brain_Hippocampus | sQTL | 0.886 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| REEP2 | scz | Brain_Cortex | eQTL | 0.006 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| REEP2 | scz | Brain_Cortex | sQTL | 0.913 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RERE | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.118 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RERE | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.818 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RERE | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.965 (susie) | 0.172 (0.040) | 2.03e-05 | 0.455 (20) | `coloc_and_smr_heidi_not_rejected` |
| RERE | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.856 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SETD6 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.951 (susie) | -0.182 (0.046) | 7.37e-05 | 0.005 (16) | `coloc_and_smr_heidi_rejected` |
| SETD6 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.926 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SNAP91 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.563 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNAP91 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.954 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| SNAP91 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.041 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNAP91 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.955 (susie) | 0.162 (0.039) | 2.96e-05 | 0.193 (15) | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | Brain_Cerebellum | eQTL | 0.030 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNAP91 | scz | Brain_Cerebellum | sQTL | 0.974 (susie) | 0.066 (0.013) | 2.85e-07 | 0.261 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | Brain_Cortex | eQTL | 0.296 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNAP91 | scz | Brain_Cortex | sQTL | 0.958 (abf) | 0.125 (0.031) | 5.58e-05 | 0.709 (18) | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | Brain_Hypothalamus | eQTL | 0.071 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNAP91 | scz | Brain_Hypothalamus | sQTL | 0.961 (abf) | 0.123 (0.029) | 2.47e-05 | 0.084 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.223 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNAP91 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.964 (abf) | 0.136 (0.030) | 5.87e-06 | 0.177 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNAP91 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.015 (susie) | 0.178 (0.056) | 0.001 | 0.689 (16) | `neither_tested_null` |
| SNAP91 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.962 (abf) | 0.017 (0.015) | 0.284 | 0.149 (5) | `coloc_smr_not_significant` |
| SYT5 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.519 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SYT5 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.904 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TAOK2 | scz | Brain_Cerebellum | eQTL | 0.105 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TAOK2 | scz | Brain_Cerebellum | sQTL | 0.829 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TEAD4 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.031 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TEAD4 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.908 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| THAP3 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.037 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| THAP3 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.964 (susie) | 0.084 (0.021) | 7e-05 | 0.017 (20) | `coloc_and_smr_heidi_not_rejected` |
| THAP3 | scz | Brain_Cerebellum | eQTL | 0.186 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| THAP3 | scz | Brain_Cerebellum | sQTL | 0.860 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Amygdala | eQTL | 0.833 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Amygdala | sQTL | 0.936 (abf) | -0.059 (0.019) | 0.002 | 0.253 (12) | `coloc_smr_not_significant` |
| TMED4 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.137 (abf) | -0.110 (0.035) | 0.001 | 0.032 (14) | `neither_tested_null` |
| TMED4 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.850 (abf) | -0.072 (0.019) | 0.000188 | 0.481 (11) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.240 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.939 (abf) | -0.077 (0.019) | 3.99e-05 | 0.019 (15) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.173 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.937 (abf) | -0.072 (0.018) | 6.89e-05 | 0.775 (13) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Cerebellum | eQTL | 0.142 (abf) | -0.143 (0.046) | 0.002 | 0.197 (20) | `neither_tested_null` |
| TMED4 | scz | Brain_Cerebellum | sQTL | 0.943 (abf) | 0.072 (0.017) | 2.43e-05 | 0.802 (15) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Cortex | eQTL | 0.274 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Cortex | sQTL | 0.938 (abf) | -0.084 (0.023) | 0.000223 | 0.689 (11) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.569 (abf) | -0.130 (0.039) | 0.000807 | 0.549 (16) | `neither_tested_null` |
| TMED4 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.937 (abf) | -0.090 (0.024) | 0.000144 | 0.208 (15) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Hippocampus | eQTL | 0.889 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Hippocampus | sQTL | 0.939 (abf) | 0.091 (0.025) | 0.000291 | 0.875 (10) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Hypothalamus | eQTL | 0.440 (abf) | -0.162 (0.052) | 0.002 | 0.367 (18) | `neither_tested_null` |
| TMED4 | scz | Brain_Hypothalamus | sQTL | 0.956 (abf) | -0.069 (0.017) | 3.41e-05 | 0.271 (14) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.603 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.941 (abf) | -0.070 (0.017) | 4.88e-05 | 0.840 (16) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.896 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.940 (abf) | -0.064 (0.015) | 1.84e-05 | 0.559 (19) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.951 (abf) | -0.175 (0.047) | 0.000191 | 0.692 (11) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.922 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Substantia_nigra | eQTL | 0.923 (abf) | -0.034 (0.031) | 0.264 | 0.653 (5) | `coloc_smr_not_significant` |
| TMED4 | scz | Brain_Substantia_nigra | sQTL | 0.914 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TSPAN31 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.131 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TSPAN31 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.813 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TUBGCP4 | scz | Brain_Cerebellum | eQTL | 0.026 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TUBGCP4 | scz | Brain_Cerebellum | eQTL | 0.026 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TUBGCP4 | scz | Brain_Cerebellum | sQTL | 0.882 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TUBGCP4 | scz | Brain_Cerebellum | sQTL | 0.882 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YPEL1 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.399 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.899 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YPEL1 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.020 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.984 (abf) | 0.071 (0.018) | 0.000107 | 0.458 (20) | `coloc_and_smr_heidi_not_rejected` |
| YPEL1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.238 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.985 (abf) | 0.074 (0.020) | 0.000248 | 0.832 (20) | `coloc_and_smr_heidi_not_rejected` |
| YPEL1 | scz | Brain_Cerebellum | eQTL | 0.040 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Cerebellum | sQTL | 0.985 (abf) | 0.065 (0.017) | 0.000152 | 0.450 (20) | `coloc_and_smr_heidi_not_rejected` |
| YPEL1 | scz | Brain_Cortex | eQTL | 0.410 (abf) | -0.141 (0.048) | 0.003 | 0.839 (20) | `neither_tested_null` |
| YPEL1 | scz | Brain_Cortex | sQTL | 0.866 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YPEL1 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.323 (abf) | -0.147 (0.055) | 0.007 | 0.152 (20) | `neither_tested_null` |
| YPEL1 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.982 (abf) | 0.093 (0.027) | 0.000472 | 0.782 (12) | `coloc_smr_not_significant` |
| YPEL1 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.904 (abf) | 0.055 (0.017) | 0.000974 | 0.435 (20) | `coloc_smr_not_significant` |
| YPEL1 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.982 (abf) | 0.081 (0.022) | 0.000218 | 0.316 (20) | `coloc_and_smr_heidi_not_rejected` |
| YPEL3 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.000496 (susie) | -0.188 (0.047) | 5.48e-05 | 0.003 (20) | `smr_without_coloc` |
| YPEL3 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.845 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Caudate_basal_ganglia | eQTL | 0.851 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Caudate_basal_ganglia | sQTL | 0.832 (abf) | 0.069 (0.019) | 0.000234 | 0.038 (11) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.821 (abf) | 0.282 (0.079) | 0.000376 | 0.125 (10) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.860 (abf) | 0.043 (0.013) | 0.001 | 0.063 (20) | `coloc_smr_not_significant` |
| YWHAB | scz | Brain_Cerebellum | eQTL | 0.104 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Cerebellum | sQTL | 0.834 (abf) | 0.054 (0.014) | 0.000107 | 0.346 (20) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Cortex | eQTL | 0.813 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Cortex | sQTL | 0.856 (abf) | 0.075 (0.021) | 0.000451 | 0.059 (8) | `coloc_smr_not_significant` |
| YWHAB | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.553 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.891 (abf) | 0.072 (0.018) | 9.61e-05 | 0.059 (11) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Hippocampus | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Hippocampus | sQTL | 0.803 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Hypothalamus | eQTL | 0.846 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Hypothalamus | sQTL | 0.802 (abf) | 0.079 (0.022) | 0.00034 | 0.154 (11) | `coloc_smr_not_significant` |
| YWHAB | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.219 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.881 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.229 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.861 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| ZDHHC12 | scz | Brain_Amygdala | eQTL | 0.053 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Amygdala | sQTL | 0.928 (abf) | -0.052 (0.014) | 0.000118 | 0.077 (19) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.054 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.932 (abf) | 0.055 (0.014) | 8.13e-05 | 0.270 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.029 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.931 (abf) | -0.041 (0.010) | 3.67e-05 | 0.207 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.062 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.931 (abf) | -0.054 (0.014) | 0.000102 | 0.192 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Cerebellum | eQTL | 0.466 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Cerebellum | sQTL | 0.915 (abf) | -0.081 (0.022) | 0.000274 | 0.309 (17) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Cortex | eQTL | 0.030 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Cortex | sQTL | 0.927 (abf) | 0.050 (0.012) | 5.75e-05 | 0.094 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.049 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.923 (abf) | 0.048 (0.012) | 6.3e-05 | 0.204 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Hippocampus | eQTL | 0.041 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Hippocampus | sQTL | 0.930 (abf) | -0.041 (0.010) | 3.84e-05 | 0.147 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Hypothalamus | eQTL | 0.046 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Hypothalamus | sQTL | 0.884 (abf) | -0.055 (0.015) | 0.000256 | 0.105 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.731 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.931 (abf) | 0.042 (0.010) | 2.87e-05 | 0.093 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.056 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.923 (abf) | 0.041 (0.010) | 4.39e-05 | 0.160 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.172 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.998 (abf) | -0.041 (0.010) | 4.05e-05 | 0.394 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZDHHC12 | scz | Brain_Substantia_nigra | eQTL | 0.061 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZDHHC12 | scz | Brain_Substantia_nigra | sQTL | 0.957 (abf) | -0.046 (0.012) | 7.42e-05 | 0.491 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZFYVE21 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZFYVE21 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZFYVE21 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.959 (abf) | -0.220 (0.046) | 1.52e-06 | 0.658 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZFYVE21 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.955 (susie) | -0.220 (0.046) | 1.52e-06 | 0.658 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZFYVE21 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.113 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZFYVE21 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.955 (susie) | -0.172 (0.037) | 3.76e-06 | 0.910 (20) | `coloc_and_smr_heidi_not_rejected` |
| ZFYVE21 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.040 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ZFYVE21 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.924 (susie) | 0.169 (0.039) | 1.9e-05 | 0.482 (20) | `coloc_and_smr_heidi_not_rejected` |

## The multi-SNP arm

`--smr-multi` combines the cis SNPs surviving LD pruning at r2 0.1 instead of testing the top SNP alone. It is a robustness arm for the SMR estimate, not a second discovery pass: it answers whether a signal rests on one lead SNP or on the cis signal more broadly.

`p_SMR` and `smr_status` in this table are still the single-SNP quantities, so rows match the primary arm one for one. The multi-SNP verdict is `p_SMR_multi` / `smr_multi_status`, Bonferroni-corrected over the probes the multi-SNP test actually ran on (`n_multi_family`) -- a smaller denominator than the instrumented count, because SMR skips the test where too few cis SNPs survive pruning. Those probes are `multi_unavailable`, which is **not** a null result.

| smr_multi_status | probes |
|---|---|
| `no_instrument` | 4,661 |
| `multi_tested_null` | 321 |
| `multi_signal_heidi_unavailable` | 4 |
| `multi_signal_heidi_rejects` | 51 |
| `multi_heidi_supported` | 590 |

Of the 642 probes significant on the single-SNP test, 631 are also significant under the multi-SNP test and 0 could not be tested. A probe that does not survive is a signal carried by its lead SNP alone; that is a caveat on the SMR estimate, not a refutation of the colocalization.

Non-primary sQTL probes (the gene's other introns) are in `smr_results.parquet` with `primary_probe = False`; they are reported so the SMR evidence is not a maximum over introns, and they are not summarized here.

