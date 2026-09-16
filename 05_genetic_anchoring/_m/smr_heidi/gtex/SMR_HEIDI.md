# SMR + HEIDI on the signal-level colocalization nominations

QTL source: `gtex`. Instrument threshold `--peqtl-smr 5e-08`. HEIDI SNPs at p < 0.0015654, 3-20 SNPs, 2000 kb cis window. Significance: Bonferroni at 0.05 over the instrumented probes of the row's OWN family (primary confirmatory / secondary event-localization), never pooled across the two. HEIDI rejects a single shared variant at p < 0.01. LD reference: 1000G EUR.

## How these numbers may be read

- `b_SMR` is a signed ratio estimate relating genetically predicted phenotype to disease under a single-causal-variant, no-pleiotropy model. It does **not** establish causal direction and cannot distinguish causality from horizontal pleiotropy.
- An sQTL probe is a LeafCutter intron-excision ratio, which is compositional within its cluster: introns sharing a splice site trade usage, so sibling probes carry opposite `b_SMR` signs by construction. A sign is read relative to its cluster, never alone.
- Failing to reject HEIDI is **not** evidence of a shared variant; HEIDI is underpowered at GTEx brain sample sizes. It reads "not rejected".
- A HEIDI rejection does **not** overrule a strong signal-level colocalization, and SMR significance does not promote a locus coloc did not support. Disagreements stay disagreements.

## What was testable, by family

Two families, corrected apart. The **primary confirmatory family** holds one pre-designated probe per gene; a gene's other introns form a **secondary event-localization family** with its own Bonferroni correction, so the primary threshold is not inflated by introns and the introns do not escape correction when they are discussed. `F` is the instrument strength `(b_eQTL/se_eQTL)^2` of the top cis-QTL SNP, reported on every row and never used to exclude one. It is summarized by its 5th percentile rather than a count below the conventional F < 10: an instrument that clears p < 5e-8 has |z| ≳ 5.4 and so F ≳ 30 (≈ 24 at the relaxed 1e-6 arm), so a weak-instrument count is zero by construction and carries no information; the per-row `weak_instrument` flag stays in `smr_results.parquet` for any run at a looser threshold.

| analysis | modality | family | probes | instrumented | threshold | F median [min-max] | F p5 | `no_instrument` | `instrumented_tested_null` | `smr_signal_heidi_unavailable` | `smr_signal_heidi_rejects` | `smr_heidi_supported` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | eQTL | primary | 17 | 3 | 0.017 | 52 [50-59] | 50 | 14 | 0 | 0 | 1 | 2 |
| aging__ad | sQTL | primary | 17 | 16 | 0.003 | 395 [38-531] | 41 | 1 | 0 | 0 | 0 | 16 |
| aging__ad | sQTL | secondary | 157 | 29 | 0.002 | 197 [32-524] | 54 | 128 | 3 | 0 | 0 | 26 |
| aging__als | eQTL | primary | 16 | 2 | 0.025 | 137 [135-139] | 135 | 14 | 0 | 0 | 0 | 2 |
| aging__als | sQTL | primary | 16 | 15 | 0.003 | 44 [32-135] | 34 | 1 | 0 | 0 | 0 | 15 |
| aging__als | sQTL | secondary | 229 | 5 | 0.010 | 35 [30-84] | 31 | 224 | 0 | 0 | 0 | 5 |
| aging__pd | eQTL | primary | 15 | 7 | 0.007 | 47 [32-55] | 33 | 8 | 0 | 0 | 0 | 7 |
| aging__pd | sQTL | primary | 15 | 13 | 0.004 | 104 [33-309] | 38 | 2 | 0 | 0 | 0 | 13 |
| aging__pd | sQTL | secondary | 127 | 7 | 0.007 | 82 [30-113] | 33 | 120 | 1 | 0 | 0 | 6 |
| aging__scz | eQTL | primary | 101 | 34 | 0.001 | 56 [30-344] | 31 | 67 | 7 | 0 | 4 | 23 |
| aging__scz | sQTL | primary | 101 | 77 | 0.000649 | 55 [30-261] | 32 | 24 | 6 | 0 | 4 | 67 |
| aging__scz | sQTL | secondary | 1274 | 78 | 0.000641 | 49 [30-188] | 30 | 1196 | 30 | 0 | 1 | 47 |
| brainseq-sczd__scz | eQTL | primary | 15 | 7 | 0.007 | 50 [31-344] | 32 | 8 | 2 | 0 | 0 | 5 |
| brainseq-sczd__scz | sQTL | primary | 15 | 13 | 0.004 | 57 [30-225] | 33 | 2 | 2 | 0 | 1 | 10 |
| brainseq-sczd__scz | sQTL | secondary | 170 | 22 | 0.002 | 48 [30-177] | 30 | 148 | 10 | 0 | 0 | 12 |

**`no_instrument` is not a negative result** — the probe was never tested, because no cis-QTL reached the instrument threshold. It is the largest cell in every arm here and must never be read as evidence against a locus.

## Agreement with coloc, primary confirmatory family

| analysis | modality | probes | instrumented | threshold | `coloc_and_smr_heidi_not_rejected` | `coloc_and_smr_heidi_rejected` | `coloc_and_smr_heidi_untestable` | `coloc_smr_not_significant` | `coloc_no_instrument` | `smr_without_coloc` | `smr_coloc_not_scored` | `neither_tested_null` | `neither_no_instrument` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | eQTL | 17 | 3 | 0.017 | 1 | 0 | 0 | 0 | 6 | 2 | 0 | 0 | 8 |
| aging__ad | sQTL | 17 | 16 | 0.003 | 16 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| aging__als | eQTL | 16 | 2 | 0.025 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 14 |
| aging__als | sQTL | 16 | 15 | 0.003 | 15 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| aging__pd | eQTL | 15 | 7 | 0.007 | 2 | 0 | 0 | 0 | 1 | 5 | 0 | 0 | 7 |
| aging__pd | sQTL | 15 | 13 | 0.004 | 13 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 |
| aging__scz | eQTL | 101 | 34 | 0.001 | 19 | 0 | 0 | 0 | 12 | 8 | 0 | 7 | 55 |
| aging__scz | sQTL | 101 | 77 | 0.000649 | 67 | 4 | 0 | 6 | 24 | 0 | 0 | 0 | 0 |
| brainseq-sczd__scz | eQTL | 15 | 7 | 0.007 | 5 | 0 | 0 | 0 | 5 | 0 | 0 | 2 | 3 |
| brainseq-sczd__scz | sQTL | 15 | 13 | 0.004 | 10 | 1 | 0 | 2 | 2 | 0 | 0 | 0 | 0 |

### Secondary family (event localization)

A gene's non-primary introns, corrected within their own family. These localize an event; they are not additional confirmatory evidence for a locus.

| analysis | modality | probes | instrumented | threshold | `smr_heidi_supported` | `smr_signal_heidi_rejects` |
|---|---|---|---|---|---|---|
| aging__ad | sQTL | 157 | 29 | 0.002 | 26 | 0 |
| aging__als | sQTL | 229 | 5 | 0.010 | 5 | 0 |
| aging__pd | sQTL | 127 | 7 | 0.007 | 6 | 0 |
| aging__scz | sQTL | 1274 | 78 | 0.000641 | 47 | 1 |
| brainseq-sczd__scz | sQTL | 170 | 22 | 0.002 | 12 | 0 |

## SNP attrition into HEIDI

HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so what survives is worth stating. Per SMR run; the per-probe end of the trace is `nsnp_HEIDI` in the results table, itself capped at 20.

**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends past its locus loses the remainder. That is not a QC failure and says nothing about the QTL data.

| step | median across runs |
|---|---|
| BESD SNPs in the probe's cis window | 5,769 |
| ... found in the LD panel by rsID | 100.0% |
| ... with matching alleles | 100.0% |
| ... also carried by the GWAS | 4,300 |
| dropped by SMR beyond that (frequency check, duplicates) | 0 |

Genuine harmonization loss — SNPs present in all three and still dropped — peaks at **6** SNPs in any run. The identifier and allele steps are the QC ones, and both sit at 100%.

**Do not compare BESD-SNP retention across QTL sources.** Its denominator is imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far larger set against the same GWAS loci. The ratio measures panel density, not quality.

## Primary probes

| gene | trait | tissue | modality | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | agreement |
|---|---|---|---|---|---|---|---|---|
| NDUFS3 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.000195 (susie) | -0.162 (0.061) | 0.008 | 0.003 (20) | `smr_without_coloc` |
| NDUFS3 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.864 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
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
| TPCN1 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.819 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| TPCN1 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.861 (susie) | 0.051 (0.013) | 9.84e-05 | 0.343 (20) | `coloc_and_smr_heidi_not_rejected` |
| TPCN1 | ad | Brain_Cerebellum | eQTL | 0.243 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TPCN1 | ad | Brain_Cerebellum | sQTL | 0.835 (susie) | -0.062 (0.017) | 0.000207 | 0.745 (10) | `coloc_and_smr_heidi_not_rejected` |
| ZNF232 | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.009 (abf) | -0.039 (0.016) | 0.012 | 0.711 (9) | `smr_without_coloc` |
| ZNF232 | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.899 (abf) | -0.073 (0.017) | 1.17e-05 | 0.020 (14) | `coloc_and_smr_heidi_not_rejected` |
| GGNBP2 | als | Brain_Cerebellum | eQTL | 0.383 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GGNBP2 | als | Brain_Cerebellum | sQTL | 0.912 (susie) | 0.128 (0.035) | 0.000231 | 0.117 (13) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.826 (abf) | -0.074 (0.021) | 0.000423 | 0.814 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Cortex | eQTL | 0.010 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Cortex | sQTL | 0.964 (abf) | -0.059 (0.015) | 9.19e-05 | 0.882 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Hippocampus | eQTL | 0.034 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Hippocampus | sQTL | 0.901 (abf) | -0.097 (0.029) | 0.000879 | 0.616 (18) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.027 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.945 (abf) | -0.059 (0.015) | 8.36e-05 | 0.887 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Putamen_basal_ganglia | eQTL | 0.039 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Putamen_basal_ganglia | sQTL | 0.857 (abf) | -0.124 (0.035) | 0.000415 | 0.914 (20) | `coloc_and_smr_heidi_not_rejected` |
| PRDM2 | als | Brain_Substantia_nigra | eQTL | 0.154 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PRDM2 | als | Brain_Substantia_nigra | sQTL | 0.947 (abf) | -0.097 (0.028) | 0.000463 | 0.891 (19) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cerebellar_Hemisphere | eQTL | 0.031 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cerebellar_Hemisphere | sQTL | 0.890 (susie) | -0.102 (0.026) | 0.000113 | 0.725 (19) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cerebellum | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cerebellum | sQTL | 0.890 (susie) | -0.099 (0.026) | 0.000113 | 0.962 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cortex | eQTL | 0.585 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cortex | sQTL | 0.889 (susie) | -0.120 (0.031) | 0.000118 | 0.663 (12) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Frontal_Cortex_BA9 | eQTL | 0.074 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Frontal_Cortex_BA9 | sQTL | 0.832 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PTPRN | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.964 (abf) | -0.172 (0.043) | 5.97e-05 | 0.444 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.918 (abf) | 0.193 (0.050) | 0.000125 | 0.618 (16) | `coloc_and_smr_heidi_not_rejected` |
| TXNDC15 | als | Brain_Caudate_basal_ganglia | eQTL | 0.347 (abf) | -0.109 (0.034) | 0.001 | 0.066 (20) | `smr_without_coloc` |
| TXNDC15 | als | Brain_Caudate_basal_ganglia | sQTL | 0.969 (abf) | -0.113 (0.029) | 0.000119 | 0.346 (20) | `coloc_and_smr_heidi_not_rejected` |
| UNC13A | als | Brain_Cerebellar_Hemisphere | eQTL | 0.377 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| UNC13A | als | Brain_Cerebellar_Hemisphere | sQTL | 0.936 (susie) | -0.184 (0.028) | 9.5e-11 | 0.205 (6) | `coloc_and_smr_heidi_not_rejected` |
| UNC13A | als | Brain_Cerebellum | eQTL | 0.090 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| UNC13A | als | Brain_Cerebellum | sQTL | 0.959 (susie) | -0.165 (0.026) | 9.03e-11 | 0.166 (6) | `coloc_and_smr_heidi_not_rejected` |
| WIPI2 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.020 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| WIPI2 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.893 (abf) | -0.097 (0.031) | 0.002 | 0.014 (14) | `coloc_and_smr_heidi_not_rejected` |
| CTSB | pd | Brain_Amygdala | eQTL | 0.709 (susie) | -0.247 (0.080) | 0.002 | 0.071 (20) | `smr_without_coloc` |
| CTSB | pd | Brain_Amygdala | sQTL | 0.868 (susie) | 0.078 (0.021) | 0.000254 | 0.193 (20) | `coloc_and_smr_heidi_not_rejected` |
| DDRGK1 | pd | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.902 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DDRGK1 | pd | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.884 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NCOR1 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.032 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NCOR1 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.895 (susie) | -0.169 (0.043) | 8.14e-05 | 0.184 (20) | `coloc_and_smr_heidi_not_rejected` |
| NCOR1 | pd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.775 (susie) | 0.368 (0.096) | 0.000128 | 0.154 (16) | `smr_without_coloc` |
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
| TTC19 | pd | Brain_Caudate_basal_ganglia | eQTL | 0.718 (susie) | -0.418 (0.099) | 2.55e-05 | 0.048 (20) | `smr_without_coloc` |
| TTC19 | pd | Brain_Caudate_basal_ganglia | sQTL | 0.805 (susie) | -0.247 (0.065) | 0.000135 | 0.046 (20) | `coloc_and_smr_heidi_not_rejected` |
| TTC19 | pd | Brain_Putamen_basal_ganglia | eQTL | 0.732 (susie) | -0.346 (0.085) | 4.74e-05 | 0.457 (20) | `smr_without_coloc` |
| TTC19 | pd | Brain_Putamen_basal_ganglia | sQTL | 0.898 (susie) | -0.185 (0.045) | 4.13e-05 | 0.768 (20) | `coloc_and_smr_heidi_not_rejected` |
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
| CDIP1 | scz | Brain_Amygdala | eQTL | 0.624 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Amygdala | sQTL | 0.929 (susie) | 0.066 (0.017) | 0.000166 | 0.726 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.715 (susie) | 0.192 (0.055) | 0.000455 | 0.088 (20) | `smr_without_coloc` |
| CDIP1 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.951 (susie) | -0.099 (0.025) | 7.18e-05 | 0.563 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Cortex | eQTL | 0.803 (susie) | 0.108 (0.027) | 6.92e-05 | 0.548 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Cortex | sQTL | 0.831 (susie) | 0.065 (0.016) | 4.86e-05 | 0.448 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Hippocampus | eQTL | 0.021 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Hippocampus | sQTL | 0.802 (susie) | 0.085 (0.026) | 0.000961 | 0.378 (20) | `coloc_smr_not_significant` |
| CDIP1 | scz | Brain_Hypothalamus | eQTL | 0.037 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Hypothalamus | sQTL | 0.936 (susie) | -0.089 (0.023) | 0.000132 | 0.893 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Substantia_nigra | eQTL | 0.562 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Substantia_nigra | sQTL | 0.928 (abf) | 0.083 (0.022) | 0.000138 | 0.683 (20) | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | Brain_Cerebellum | eQTL | 0.822 (abf) | -0.095 (0.025) | 0.000127 | 0.587 (20) | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | Brain_Cerebellum | sQTL | 0.830 (abf) | -0.114 (0.035) | 0.001 | 0.237 (20) | `coloc_smr_not_significant` |
| DGKZ | scz | Brain_Amygdala | eQTL | 0.072 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Amygdala | sQTL | 0.878 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.882 (susie) | -0.139 (0.027) | 3.63e-07 | 0.034 (20) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Caudate_basal_ganglia | eQTL | 0.495 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Caudate_basal_ganglia | sQTL | 0.935 (susie) | -0.194 (0.043) | 6e-06 | 0.268 (15) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.858 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.924 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Cortex | eQTL | 0.038 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Cortex | sQTL | 0.951 (susie) | -0.145 (0.028) | 1.75e-07 | 0.356 (20) | `coloc_and_smr_heidi_not_rejected` |
| DGKZ | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.953 (susie) | -0.151 (0.028) | 7.34e-08 | 0.004 (20) | `coloc_and_smr_heidi_rejected` |
| DGKZ | scz | Brain_Hippocampus | eQTL | 0.035 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Hippocampus | sQTL | 0.838 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.132 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.928 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DGKZ | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.035 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DGKZ | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.892 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| FGFR1 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FGFR1 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.908 (susie) | 0.092 (0.023) | 6.81e-05 | 0.998 (14) | `coloc_and_smr_heidi_not_rejected` |
| FGFR1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.178 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FGFR1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.957 (susie) | 0.080 (0.018) | 1.19e-05 | 0.882 (11) | `coloc_and_smr_heidi_not_rejected` |
| FGFR1 | scz | Brain_Cerebellum | eQTL | 0.653 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| FGFR1 | scz | Brain_Cerebellum | sQTL | 0.937 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| GABBR2 | scz | Brain_Cerebellum | eQTL | 0.014 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GABBR2 | scz | Brain_Cerebellum | sQTL | 0.976 (susie) | -0.100 (0.026) | 9.97e-05 | 0.604 (5) | `coloc_and_smr_heidi_not_rejected` |
| GLYCTK | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.499 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GLYCTK | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.812 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| IRF3 | scz | Brain_Amygdala | eQTL | 0.274 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Amygdala | sQTL | 0.971 (abf) | 0.077 (0.017) | 6.72e-06 | 0.096 (15) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.144 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.992 (abf) | 0.094 (0.020) | 1.57e-06 | 0.102 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.403 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| IRF3 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.988 (abf) | 0.075 (0.013) | 2.5e-08 | 0.241 (20) | `coloc_and_smr_heidi_not_rejected` |
| IRF3 | scz | Brain_Cortex | eQTL | 0.244 (abf) | -0.196 (0.060) | 0.001 | 0.002 (7) | `smr_without_coloc` |
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
| KLC1 | scz | Brain_Cortex | eQTL | 0.000335 (susie) | 0.150 (0.097) | 0.123 | 0.000978 (20) | `neither_tested_null` |
| KLC1 | scz | Brain_Cortex | sQTL | 0.890 (susie) | -0.293 (0.064) | 4.16e-06 | 0.582 (20) | `coloc_and_smr_heidi_not_rejected` |
| LPCAT4 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.007 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| LPCAT4 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.884 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
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
| NDUFAF7 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.909 (susie) | -0.055 (0.012) | 8.4e-06 | 0.171 (20) | `coloc_and_smr_heidi_not_rejected` |
| NEK4 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.341 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NEK4 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.820 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NUP50 | scz | Brain_Amygdala | eQTL | 0.895 (abf) | -0.091 (0.024) | 0.000166 | 0.678 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Amygdala | sQTL | 0.835 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NUP50 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.925 (abf) | -0.100 (0.024) | 4.51e-05 | 0.978 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.918 (abf) | 0.064 (0.017) | 0.000159 | 0.576 (17) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.972 (abf) | -0.098 (0.023) | 1.64e-05 | 0.746 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.960 (abf) | -0.094 (0.026) | 0.000307 | 0.634 (19) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.948 (abf) | -0.163 (0.044) | 0.000237 | 0.802 (18) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.945 (abf) | -0.092 (0.025) | 0.000191 | 0.831 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cerebellum | eQTL | 0.760 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NUP50 | scz | Brain_Cerebellum | sQTL | 0.953 (abf) | -0.080 (0.021) | 0.000135 | 0.888 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cortex | eQTL | 0.906 (abf) | -0.106 (0.027) | 8.16e-05 | 0.409 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Cortex | sQTL | 0.935 (abf) | 0.066 (0.018) | 0.000345 | 0.813 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.943 (abf) | -0.127 (0.032) | 6.09e-05 | 0.385 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.927 (abf) | 0.058 (0.015) | 0.000167 | 0.951 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Hippocampus | eQTL | 0.932 (abf) | -0.126 (0.031) | 4.96e-05 | 0.774 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Hippocampus | sQTL | 0.910 (abf) | 0.069 (0.020) | 0.000583 | 0.817 (17) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Hypothalamus | eQTL | 0.791 (abf) | -0.115 (0.030) | 0.000156 | 0.806 (20) | `smr_without_coloc` |
| NUP50 | scz | Brain_Hypothalamus | sQTL | 0.918 (abf) | 0.065 (0.018) | 0.000389 | 0.457 (15) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.839 (abf) | -0.088 (0.022) | 7.71e-05 | 0.660 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.956 (abf) | -0.069 (0.018) | 0.000101 | 0.848 (20) | `coloc_and_smr_heidi_not_rejected` |
| NUP50 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.733 (abf) | -0.073 (0.020) | 0.000229 | 0.657 (20) | `smr_without_coloc` |
| NUP50 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.805 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| POLG | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.038 (abf) | -0.157 (0.054) | 0.004 | 0.020 (20) | `neither_tested_null` |
| POLG | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.982 (abf) | -0.048 (0.010) | 4.79e-07 | 0.770 (20) | `coloc_and_smr_heidi_not_rejected` |
| POLG | scz | Brain_Cerebellum | eQTL | 0.047 (abf) | -0.110 (0.034) | 0.001 | 0.009 (20) | `smr_without_coloc` |
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
| PPIL2 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.855 (abf) | -0.040 (0.012) | 0.00052 | 0.032 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIL2 | scz | Brain_Substantia_nigra | eQTL | 0.867 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PPIL2 | scz | Brain_Substantia_nigra | sQTL | 0.936 (abf) | -0.053 (0.014) | 0.000186 | 0.031 (20) | `coloc_and_smr_heidi_not_rejected` |
| PSMD6 | scz | Brain_Hippocampus | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PSMD6 | scz | Brain_Hippocampus | sQTL | 0.976 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RAI1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.717 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RAI1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.945 (susie) | -0.074 (0.017) | 1.07e-05 | 0.439 (20) | `coloc_and_smr_heidi_not_rejected` |
| RAI1 | scz | Brain_Cerebellum | eQTL | 0.308 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RAI1 | scz | Brain_Cerebellum | sQTL | 0.971 (abf) | 0.092 (0.022) | 2.03e-05 | 0.414 (20) | `coloc_and_smr_heidi_not_rejected` |
| RASA1 | scz | Brain_Cerebellum | eQTL | 0.204 (abf) | 0.122 (0.051) | 0.017 | 0.273 (20) | `neither_tested_null` |
| RASA1 | scz | Brain_Cerebellum | sQTL | 0.946 (abf) | -0.047 (0.015) | 0.002 | 0.260 (20) | `coloc_smr_not_significant` |
| RERE | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.119 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| RERE | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.818 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RERE | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.965 (susie) | 0.172 (0.040) | 2.03e-05 | 0.455 (20) | `coloc_and_smr_heidi_not_rejected` |
| RERE | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.856 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
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
| SNAP91 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.016 (susie) | 0.178 (0.056) | 0.001 | 0.689 (16) | `neither_tested_null` |
| SNAP91 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.962 (abf) | 0.017 (0.015) | 0.284 | 0.149 (5) | `coloc_smr_not_significant` |
| SYT5 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.519 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SYT5 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.904 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TAOK2 | scz | Brain_Cerebellum | eQTL | 0.105 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TAOK2 | scz | Brain_Cerebellum | sQTL | 0.829 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
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
| YPEL1 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.982 (abf) | 0.093 (0.027) | 0.000472 | 0.782 (12) | `coloc_and_smr_heidi_not_rejected` |
| YPEL1 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.904 (abf) | 0.055 (0.017) | 0.000974 | 0.435 (20) | `coloc_smr_not_significant` |
| YPEL1 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YPEL1 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.982 (abf) | 0.081 (0.022) | 0.000218 | 0.316 (20) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Caudate_basal_ganglia | eQTL | 0.851 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Caudate_basal_ganglia | sQTL | 0.832 (abf) | 0.069 (0.019) | 0.000234 | 0.038 (11) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.821 (abf) | 0.282 (0.079) | 0.000376 | 0.125 (10) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.860 (abf) | 0.043 (0.013) | 0.001 | 0.063 (20) | `coloc_smr_not_significant` |
| YWHAB | scz | Brain_Cerebellum | eQTL | 0.104 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Cerebellum | sQTL | 0.834 (abf) | 0.054 (0.014) | 0.000107 | 0.346 (20) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Cortex | eQTL | 0.813 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Cortex | sQTL | 0.856 (abf) | 0.075 (0.021) | 0.000451 | 0.059 (8) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.553 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.891 (abf) | 0.072 (0.018) | 9.61e-05 | 0.059 (11) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Hippocampus | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Hippocampus | sQTL | 0.803 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Hypothalamus | eQTL | 0.846 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Hypothalamus | sQTL | 0.802 (abf) | 0.079 (0.022) | 0.00034 | 0.154 (11) | `coloc_and_smr_heidi_not_rejected` |
| YWHAB | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.219 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.881 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| YWHAB | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.229 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| YWHAB | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.861 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |

Non-primary sQTL probes (the gene's other introns) are in `smr_results.parquet` with `primary_probe = False`; they are reported so the SMR evidence is not a maximum over introns, and they are not summarized here.

