# SMR + HEIDI on the signal-level colocalization nominations

QTL source: `gtex`. Instrument threshold `--peqtl-smr 5e-08`. Multi-SNP SMR (`--smr-multi`, LD pruned at r2 0.1) -- **SENSITIVITY ARM**. HEIDI SNPs at p < 0.0015654, 3-20 SNPs, 2000 kb cis window. Significance: Bonferroni at 0.05 over the instrumented probes of the row's OWN family (primary confirmatory / secondary event-localization), never pooled across the two. HEIDI rejects a single shared variant at p < 0.01. LD reference: 1000G EUR.

## How these numbers may be read

- `b_SMR` is a signed ratio estimate relating genetically predicted phenotype to disease under a single-causal-variant, no-pleiotropy model. It does **not** establish causal direction and cannot distinguish causality from horizontal pleiotropy.
- An sQTL probe is a LeafCutter intron-excision ratio, which is compositional within its cluster: introns sharing a splice site trade usage, so sibling probes carry opposite `b_SMR` signs by construction. A sign is read relative to its cluster, never alone.
- Failing to reject HEIDI is **not** evidence of a shared variant; HEIDI is underpowered at GTEx brain sample sizes. It reads "not rejected".
- A HEIDI rejection does **not** overrule a strong signal-level colocalization, and SMR significance does not promote a locus coloc did not support. Disagreements stay disagreements.

## What was testable, by family

Two families, corrected apart. The **primary confirmatory family** holds one pre-designated probe per gene; a gene's other introns form a **secondary event-localization family** with its own Bonferroni correction, so the primary threshold is not inflated by introns and the introns do not escape correction when they are discussed. `F` is the instrument strength `(b_eQTL/se_eQTL)^2` of the top cis-QTL SNP, reported on every row and never used to exclude one.

| analysis | modality | family | probes | instrumented | threshold | F median [min-max] | weak F | `no_instrument` | `instrumented_tested_null` | `smr_signal_heidi_unavailable` | `smr_signal_heidi_rejects` | `smr_heidi_supported` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | eQTL | primary | 34 | 11 | 0.005 | 42 [31-59] | 0 | 23 | 6 | 0 | 0 | 5 |
| aging__ad | sQTL | primary | 34 | 31 | 0.002 | 114 [34-531] | 0 | 3 | 0 | 0 | 0 | 31 |
| aging__ad | sQTL | secondary | 417 | 46 | 0.001 | 107 [30-524] | 0 | 371 | 6 | 0 | 0 | 40 |
| aging__als | eQTL | primary | 26 | 14 | 0.004 | 54 [32-139] | 0 | 12 | 10 | 0 | 1 | 3 |
| aging__als | sQTL | primary | 26 | 24 | 0.002 | 63 [30-192] | 0 | 2 | 0 | 0 | 0 | 24 |
| aging__als | sQTL | secondary | 376 | 9 | 0.006 | 44 [31-143] | 0 | 367 | 1 | 0 | 0 | 8 |
| aging__lbd | eQTL | primary | 11 | 0 | — | — | 0 | 11 | 0 | 0 | 0 | 0 |
| aging__lbd | sQTL | primary | 11 | 11 | 0.005 | 105 [39-221] | 0 | 0 | 0 | 0 | 0 | 11 |
| aging__lbd | sQTL | secondary | 81 | 21 | 0.002 | 66 [32-203] | 0 | 60 | 9 | 0 | 11 | 1 |
| aging__pd | eQTL | primary | 23 | 7 | 0.007 | 53 [33-72] | 0 | 16 | 1 | 0 | 0 | 6 |
| aging__pd | sQTL | primary | 23 | 20 | 0.003 | 104 [33-309] | 0 | 3 | 2 | 0 | 6 | 12 |
| aging__pd | sQTL | secondary | 200 | 19 | 0.003 | 62 [30-135] | 0 | 181 | 1 | 0 | 9 | 9 |
| aging__scz | eQTL | primary | 67 | 25 | 0.002 | 55 [30-124] | 0 | 42 | 10 | 0 | 1 | 14 |
| aging__scz | sQTL | primary | 67 | 50 | 0.001 | 44 [30-291] | 0 | 17 | 3 | 0 | 7 | 40 |
| aging__scz | sQTL | secondary | 745 | 50 | 0.001 | 53 [30-440] | 0 | 695 | 18 | 0 | 7 | 25 |

**`no_instrument` is not a negative result** — the probe was never tested, because no cis-QTL reached the instrument threshold. It is the largest cell in every arm here and must never be read as evidence against a locus.

## Agreement with coloc, primary confirmatory family

| analysis | modality | probes | instrumented | threshold | `coloc_and_smr_heidi_not_rejected` | `coloc_and_smr_heidi_rejected` | `coloc_and_smr_heidi_untestable` | `coloc_smr_not_significant` | `coloc_no_instrument` | `smr_without_coloc` | `smr_coloc_not_scored` | `neither_tested_null` | `neither_no_instrument` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | eQTL | 34 | 11 | 0.005 | 5 | 0 | 0 | 2 | 7 | 0 | 0 | 4 | 16 |
| aging__ad | sQTL | 34 | 31 | 0.002 | 31 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 |
| aging__als | eQTL | 26 | 14 | 0.004 | 1 | 0 | 0 | 0 | 1 | 3 | 0 | 10 | 11 |
| aging__als | sQTL | 26 | 24 | 0.002 | 24 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 |
| aging__lbd | eQTL | 11 | 0 | — | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 11 |
| aging__lbd | sQTL | 11 | 11 | 0.005 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| aging__pd | eQTL | 23 | 7 | 0.007 | 2 | 0 | 0 | 0 | 0 | 4 | 0 | 1 | 16 |
| aging__pd | sQTL | 23 | 20 | 0.003 | 12 | 6 | 0 | 2 | 3 | 0 | 0 | 0 | 0 |
| aging__scz | eQTL | 67 | 25 | 0.002 | 9 | 1 | 0 | 1 | 5 | 5 | 0 | 9 | 37 |
| aging__scz | sQTL | 67 | 50 | 0.001 | 40 | 7 | 0 | 3 | 17 | 0 | 0 | 0 | 0 |

### Secondary family (event localization)

A gene's non-primary introns, corrected within their own family. These localize an event; they are not additional confirmatory evidence for a locus.

| analysis | modality | probes | instrumented | threshold | `smr_heidi_supported` | `smr_signal_heidi_rejects` |
|---|---|---|---|---|---|---|
| aging__ad | sQTL | 417 | 46 | 0.001 | 40 | 0 |
| aging__als | sQTL | 376 | 9 | 0.006 | 8 | 0 |
| aging__lbd | sQTL | 81 | 21 | 0.002 | 1 | 11 |
| aging__pd | sQTL | 200 | 19 | 0.003 | 9 | 9 |
| aging__scz | sQTL | 745 | 50 | 0.001 | 25 | 7 |

## SNP attrition into HEIDI

HEIDI is computed on the SNPs shared by the BESD, the LD reference and the GWAS, so what survives is worth stating. Per SMR run; the per-probe end of the trace is `nsnp_HEIDI` in the results table, itself capped at 20.

**The dominant filter is GWAS locus coverage, by construction** — the `.ma` files carry only the SNPs of the selected GWAS loci, so a gene whose cis window extends past its locus loses the remainder. That is not a QC failure and says nothing about the QTL data.

| step | median across runs |
|---|---|
| BESD SNPs in the probe's cis window | 5,769 |
| ... found in the LD panel by rsID | 100.0% |
| ... with matching alleles | 100.0% |
| ... also carried by the GWAS | 4,408 |
| dropped by SMR beyond that (frequency check, duplicates) | 0 |

Genuine harmonization loss — SNPs present in all three and still dropped — peaks at **4** SNPs in any run. The identifier and allele steps are the QC ones, and both sit at 100%.

**Do not compare BESD-SNP retention across QTL sources.** Its denominator is imputation density: BrainSEQ's TOPMed cis windows carry several times the SNPs of GTEx's pre-filtered all-pairs, so BrainSEQ retains a far smaller *fraction* of a far larger set against the same GWAS loci. The ratio measures panel density, not quality.

## Primary probes

| gene | trait | tissue | modality | coloc PP4 | b_SMR (se) | p_SMR | p_HEIDI (nsnp) | agreement |
|---|---|---|---|---|---|---|---|---|
| CTSH | ad | Brain_Hippocampus | eQTL | 0.986 (abf) | 0.117 (0.032) | 0.000206 | 0.058 (10) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Hippocampus | sQTL | 0.989 (abf) | 0.068 (0.018) | 0.000195 | 0.751 (8) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Hypothalamus | eQTL | 0.959 (abf) | -0.005 (0.029) | 0.855 | 0.631 (17) | `coloc_smr_not_significant` |
| CTSH | ad | Brain_Hypothalamus | sQTL | 0.985 (abf) | -0.078 (0.020) | 0.000101 | 0.723 (14) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.948 (abf) | 0.033 (0.025) | 0.184 | — (—) | `coloc_smr_not_significant` |
| CTSH | ad | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.988 (abf) | -0.083 (0.021) | 4.8e-05 | 0.858 (16) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Putamen_basal_ganglia | eQTL | 0.989 (abf) | 0.101 (0.027) | 0.000134 | 0.404 (10) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Putamen_basal_ganglia | sQTL | 0.977 (abf) | -0.084 (0.023) | 0.000223 | 0.288 (9) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.989 (abf) | 0.125 (0.033) | 0.000135 | 0.117 (7) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.899 (abf) | -0.077 (0.022) | 0.000436 | 0.295 (18) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Substantia_nigra | eQTL | 0.993 (abf) | 0.106 (0.028) | 0.000196 | 0.901 (5) | `coloc_and_smr_heidi_not_rejected` |
| CTSH | ad | Brain_Substantia_nigra | sQTL | 0.985 (abf) | 0.071 (0.019) | 0.000206 | 0.958 (8) | `coloc_and_smr_heidi_not_rejected` |
| DOC2A | ad | Brain_Amygdala | eQTL | 0.030 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | ad | Brain_Amygdala | sQTL | 0.881 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NDUFS3 | ad | Brain_Cerebellar_Hemisphere | eQTL | 0.000215 (susie) | -0.162 (0.061) | 0.008 | 0.003 (20) | `neither_tested_null` |
| NDUFS3 | ad | Brain_Cerebellar_Hemisphere | sQTL | 0.864 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PICALM | ad | Brain_Cortex | eQTL | 0.117 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PICALM | ad | Brain_Cortex | sQTL | 0.818 (abf) | -0.309 (0.059) | 1.49e-07 | 0.170 (20) | `coloc_and_smr_heidi_not_rejected` |
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
| TPCN1 | ad | Brain_Cerebellum | sQTL | 0.835 (susie) | -0.062 (0.017) | 0.000207 | 0.745 (10) | `coloc_and_smr_heidi_not_rejected` |
| ZNF232 | ad | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.009 (abf) | -0.039 (0.016) | 0.012 | 0.711 (9) | `neither_tested_null` |
| ZNF232 | ad | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.899 (abf) | -0.073 (0.017) | 1.17e-05 | 0.020 (14) | `coloc_and_smr_heidi_not_rejected` |
| G2E3 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.012 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| G2E3 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.839 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GGNBP2 | als | Brain_Cerebellum | eQTL | 0.384 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GGNBP2 | als | Brain_Cerebellum | sQTL | 0.911 (susie) | 0.128 (0.035) | 0.000231 | 0.117 (13) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Amygdala | eQTL | 0.033 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PGS1 | als | Brain_Amygdala | sQTL | 0.973 (abf) | -0.122 (0.034) | 0.000342 | 0.047 (13) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.015 (abf) | 0.043 (0.036) | 0.224 | 0.728 (13) | `neither_tested_null` |
| PGS1 | als | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.906 (abf) | -0.061 (0.017) | 0.000468 | 0.031 (19) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Caudate_basal_ganglia | eQTL | 0.058 (abf) | -0.115 (0.039) | 0.003 | 0.000589 (14) | `smr_without_coloc` |
| PGS1 | als | Brain_Caudate_basal_ganglia | sQTL | 0.967 (abf) | -0.093 (0.023) | 5.38e-05 | 0.019 (20) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Cerebellar_Hemisphere | eQTL | 0.020 (abf) | -0.091 (0.043) | 0.037 | 0.002 (14) | `neither_tested_null` |
| PGS1 | als | Brain_Cerebellar_Hemisphere | sQTL | 0.972 (abf) | -0.070 (0.017) | 3.4e-05 | 0.027 (20) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Cerebellum | eQTL | 0.005 (abf) | 0.045 (0.035) | 0.199 | 0.164 (17) | `neither_tested_null` |
| PGS1 | als | Brain_Cerebellum | sQTL | 0.969 (abf) | -0.080 (0.020) | 4.93e-05 | 0.064 (20) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Cortex | eQTL | 0.032 (abf) | -0.096 (0.040) | 0.018 | 0.002 (16) | `neither_tested_null` |
| PGS1 | als | Brain_Cortex | sQTL | 0.976 (abf) | -0.111 (0.029) | 0.000126 | 0.042 (18) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Frontal_Cortex_BA9 | eQTL | 0.165 (abf) | -0.111 (0.039) | 0.004 | 0.002 (16) | `neither_tested_null` |
| PGS1 | als | Brain_Frontal_Cortex_BA9 | sQTL | 0.942 (abf) | -0.065 (0.019) | 0.000518 | 0.015 (18) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Hippocampus | eQTL | 0.006 (abf) | 0.049 (0.045) | 0.270 | 0.180 (16) | `neither_tested_null` |
| PGS1 | als | Brain_Hippocampus | sQTL | 0.974 (abf) | -0.090 (0.022) | 5.49e-05 | 0.071 (19) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Hypothalamus | eQTL | 0.006 (abf) | 0.044 (0.046) | 0.338 | 0.149 (13) | `neither_tested_null` |
| PGS1 | als | Brain_Hypothalamus | sQTL | 0.954 (abf) | -0.080 (0.024) | 0.000773 | 0.029 (17) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.007 (abf) | 0.048 (0.031) | 0.123 | 0.713 (15) | `neither_tested_null` |
| PGS1 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.945 (abf) | -0.106 (0.027) | 8.68e-05 | 0.028 (20) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Putamen_basal_ganglia | eQTL | 0.006 (abf) | 0.042 (0.032) | 0.179 | 0.587 (13) | `neither_tested_null` |
| PGS1 | als | Brain_Putamen_basal_ganglia | sQTL | 0.970 (abf) | -0.129 (0.034) | 0.000168 | 0.020 (20) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.006 (abf) | 0.040 (0.032) | 0.223 | 0.734 (15) | `neither_tested_null` |
| PGS1 | als | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.968 (abf) | -0.087 (0.022) | 6.95e-05 | 0.044 (20) | `coloc_and_smr_heidi_not_rejected` |
| PGS1 | als | Brain_Substantia_nigra | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PGS1 | als | Brain_Substantia_nigra | sQTL | 0.970 (abf) | -0.093 (0.024) | 8.62e-05 | 0.044 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cerebellar_Hemisphere | eQTL | 0.031 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cerebellar_Hemisphere | sQTL | 0.876 (susie) | 0.121 (0.034) | 0.000323 | 0.564 (8) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cerebellum | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cerebellum | sQTL | 0.875 (susie) | -0.099 (0.026) | 0.000113 | 0.962 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Cortex | eQTL | 0.585 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Cortex | sQTL | 0.874 (susie) | -0.120 (0.031) | 0.000118 | 0.663 (12) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Frontal_Cortex_BA9 | eQTL | 0.074 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PTPRN | als | Brain_Frontal_Cortex_BA9 | sQTL | 0.832 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PTPRN | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.964 (abf) | -0.172 (0.043) | 5.97e-05 | 0.444 (20) | `coloc_and_smr_heidi_not_rejected` |
| PTPRN | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.918 (abf) | 0.193 (0.050) | 0.000125 | 0.618 (16) | `coloc_and_smr_heidi_not_rejected` |
| SCFD1 | als | Brain_Hypothalamus | eQTL | 0.929 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| SCFD1 | als | Brain_Hypothalamus | sQTL | 0.856 (abf) | -0.184 (0.042) | 9.15e-06 | 0.401 (17) | `coloc_and_smr_heidi_not_rejected` |
| TPP1 | als | Brain_Cerebellum | eQTL | 0.774 (abf) | 0.213 (0.064) | 0.000948 | 0.538 (7) | `smr_without_coloc` |
| TPP1 | als | Brain_Cerebellum | sQTL | 0.996 (abf) | -0.124 (0.030) | 2.99e-05 | 0.950 (11) | `coloc_and_smr_heidi_not_rejected` |
| TXNDC15 | als | Brain_Caudate_basal_ganglia | eQTL | 0.347 (abf) | -0.109 (0.034) | 0.001 | 0.066 (20) | `smr_without_coloc` |
| TXNDC15 | als | Brain_Caudate_basal_ganglia | sQTL | 0.969 (abf) | -0.113 (0.029) | 0.000119 | 0.346 (20) | `coloc_and_smr_heidi_not_rejected` |
| UNC13A | als | Brain_Cerebellar_Hemisphere | eQTL | 0.377 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| UNC13A | als | Brain_Cerebellar_Hemisphere | sQTL | 0.939 (susie) | -0.184 (0.028) | 9.5e-11 | 0.205 (6) | `coloc_and_smr_heidi_not_rejected` |
| UNC13A | als | Brain_Cerebellum | eQTL | 0.090 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| UNC13A | als | Brain_Cerebellum | sQTL | 0.961 (susie) | -0.165 (0.026) | 9.03e-11 | 0.166 (6) | `coloc_and_smr_heidi_not_rejected` |
| WIPI2 | als | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.020 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| WIPI2 | als | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.893 (abf) | -0.097 (0.031) | 0.002 | 0.014 (14) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Amygdala | eQTL | 0.113 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Amygdala | sQTL | 0.955 (susie) | 0.245 (0.044) | 2.93e-08 | 0.866 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.103 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.952 (susie) | 0.218 (0.037) | 3.43e-09 | 0.177 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Cerebellar_Hemisphere | eQTL | 0.009 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Cerebellar_Hemisphere | sQTL | 0.952 (susie) | 0.267 (0.047) | 1.28e-08 | 0.479 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Cerebellum | eQTL | 0.067 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Cerebellum | sQTL | 0.941 (susie) | 0.288 (0.055) | 1.35e-07 | 0.627 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Cortex | eQTL | 0.009 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Cortex | sQTL | 0.974 (susie) | 0.324 (0.059) | 4.88e-08 | 0.294 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Frontal_Cortex_BA9 | eQTL | 0.040 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Frontal_Cortex_BA9 | sQTL | 0.956 (susie) | 0.275 (0.049) | 2.22e-08 | 0.101 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Hippocampus | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Hippocampus | sQTL | 0.954 (susie) | 0.307 (0.056) | 5.63e-08 | 0.206 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Hypothalamus | eQTL | 0.069 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Hypothalamus | sQTL | 0.971 (susie) | 0.259 (0.046) | 2.46e-08 | 0.500 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.045 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.956 (susie) | 0.743 (0.166) | 7.45e-06 | 0.094 (17) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.044 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.962 (susie) | 0.279 (0.052) | 7.17e-08 | 0.379 (20) | `coloc_and_smr_heidi_not_rejected` |
| SNCA | lbd | Brain_Substantia_nigra | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | lbd | Brain_Substantia_nigra | sQTL | 0.971 (susie) | 0.303 (0.059) | 2.78e-07 | 0.396 (20) | `coloc_and_smr_heidi_not_rejected` |
| AZI2 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.015 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| AZI2 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.803 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NCOR1 | pd | Brain_Cerebellar_Hemisphere | eQTL | 0.032 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NCOR1 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.895 (susie) | -0.169 (0.043) | 8.14e-05 | 0.184 (20) | `coloc_and_smr_heidi_not_rejected` |
| NCOR1 | pd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.774 (susie) | 0.368 (0.096) | 0.000128 | 0.154 (16) | `smr_without_coloc` |
| NCOR1 | pd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.805 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PITPNM2 | pd | Brain_Cerebellar_Hemisphere | eQTL | 2.84e-08 (susie) | 0.090 (0.054) | 0.096 | 0.104 (20) | `neither_tested_null` |
| PITPNM2 | pd | Brain_Cerebellar_Hemisphere | sQTL | 0.989 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
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
| SNCA | pd | Brain_Cerebellum | eQTL | 0.058 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Cerebellum | sQTL | 0.974 (susie) | 0.072 (0.022) | 0.001 | 8.35e-09 (20) | `coloc_and_smr_heidi_rejected` |
| SNCA | pd | Brain_Cortex | eQTL | 0.057 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Cortex | sQTL | 0.959 (susie) | 0.074 (0.024) | 0.002 | 3.93e-12 (20) | `coloc_and_smr_heidi_rejected` |
| SNCA | pd | Brain_Frontal_Cortex_BA9 | eQTL | 0.007 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Frontal_Cortex_BA9 | sQTL | 0.968 (susie) | 0.069 (0.021) | 0.001 | 1.21e-13 (20) | `coloc_and_smr_heidi_rejected` |
| SNCA | pd | Brain_Hippocampus | eQTL | 0.016 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Hippocampus | sQTL | 0.894 (susie) | 0.070 (0.023) | 0.002 | 2.94e-11 (20) | `coloc_and_smr_heidi_rejected` |
| SNCA | pd | Brain_Hypothalamus | eQTL | 0.026 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Hypothalamus | sQTL | 0.903 (susie) | 0.060 (0.019) | 0.002 | 1.22e-12 (20) | `coloc_and_smr_heidi_rejected` |
| SNCA | pd | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.005 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.950 (susie) | 0.171 (0.060) | 0.004 | 1.47e-07 (17) | `coloc_smr_not_significant` |
| SNCA | pd | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.382 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.894 (susie) | 0.064 (0.021) | 0.002 | 1.09e-10 (20) | `coloc_and_smr_heidi_rejected` |
| SNCA | pd | Brain_Substantia_nigra | eQTL | 0.011 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SNCA | pd | Brain_Substantia_nigra | sQTL | 0.970 (susie) | 0.070 (0.023) | 0.003 | 4.15e-07 (20) | `coloc_smr_not_significant` |
| TTC19 | pd | Brain_Caudate_basal_ganglia | eQTL | 0.717 (susie) | -0.418 (0.099) | 2.55e-05 | 0.048 (20) | `smr_without_coloc` |
| TTC19 | pd | Brain_Caudate_basal_ganglia | sQTL | 0.805 (susie) | -0.247 (0.065) | 0.000135 | 0.046 (20) | `coloc_and_smr_heidi_not_rejected` |
| TTC19 | pd | Brain_Putamen_basal_ganglia | eQTL | 0.732 (susie) | -0.346 (0.085) | 4.74e-05 | 0.457 (20) | `smr_without_coloc` |
| TTC19 | pd | Brain_Putamen_basal_ganglia | sQTL | 0.898 (susie) | -0.185 (0.045) | 4.13e-05 | 0.768 (20) | `coloc_and_smr_heidi_not_rejected` |
| ASB3 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.150 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| ASB3 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.859 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| CDIP1 | scz | Brain_Amygdala | eQTL | 0.624 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Amygdala | sQTL | 0.928 (susie) | 0.066 (0.017) | 0.000166 | 0.726 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.723 (susie) | 0.192 (0.055) | 0.000455 | 0.088 (20) | `smr_without_coloc` |
| CDIP1 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.949 (susie) | -0.099 (0.025) | 7.18e-05 | 0.563 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Cortex | eQTL | 0.808 (susie) | 0.108 (0.027) | 6.92e-05 | 0.548 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Cortex | sQTL | 0.842 (susie) | 0.065 (0.016) | 4.86e-05 | 0.448 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Hippocampus | eQTL | 0.021 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Hippocampus | sQTL | 0.804 (susie) | 0.085 (0.026) | 0.000961 | 0.378 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Hypothalamus | eQTL | 0.037 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Hypothalamus | sQTL | 0.934 (susie) | -0.089 (0.023) | 0.000132 | 0.893 (20) | `coloc_and_smr_heidi_not_rejected` |
| CDIP1 | scz | Brain_Substantia_nigra | eQTL | 0.562 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| CDIP1 | scz | Brain_Substantia_nigra | sQTL | 0.928 (abf) | 0.083 (0.022) | 0.000138 | 0.683 (20) | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | Brain_Cerebellum | eQTL | 0.822 (abf) | -0.095 (0.025) | 0.000127 | 0.587 (20) | `coloc_and_smr_heidi_not_rejected` |
| COPA | scz | Brain_Cerebellum | sQTL | 0.830 (abf) | -0.114 (0.035) | 0.001 | 0.237 (20) | `coloc_smr_not_significant` |
| DOC2A | scz | Brain_Amygdala | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Amygdala | sQTL | 0.947 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DOC2A | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.009 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.959 (abf) | -0.140 (0.028) | 9.29e-07 | 0.061 (20) | `coloc_and_smr_heidi_not_rejected` |
| DOC2A | scz | Brain_Caudate_basal_ganglia | eQTL | 0.671 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| DOC2A | scz | Brain_Caudate_basal_ganglia | sQTL | 0.871 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| DOC2A | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.035 (susie) | — (—) | — | — (—) | `neither_no_instrument` |
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
| FAM221A | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.001 (susie) | 0.036 (0.026) | 0.165 | 0.008 (20) | `neither_tested_null` |
| FAM221A | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.817 (susie) | -0.041 (0.009) | 3.97e-06 | 0.003 (20) | `coloc_and_smr_heidi_rejected` |
| FAM221A | scz | Brain_Cortex | eQTL | 0.008 (susie) | 0.078 (0.028) | 0.005 | 0.288 (20) | `neither_tested_null` |
| FAM221A | scz | Brain_Cortex | sQTL | 0.803 (susie) | 0.040 (0.009) | 6.86e-06 | 0.003 (20) | `coloc_and_smr_heidi_rejected` |
| FAM221A | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.009 (susie) | 0.059 (0.019) | 0.003 | 0.033 (20) | `neither_tested_null` |
| FAM221A | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.812 (susie) | -0.056 (0.018) | 0.002 | 0.145 (20) | `coloc_smr_not_significant` |
| GABBR2 | scz | Brain_Cerebellum | eQTL | 0.014 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GABBR2 | scz | Brain_Cerebellum | sQTL | 0.976 (susie) | -0.100 (0.026) | 9.97e-05 | 0.604 (5) | `coloc_and_smr_heidi_not_rejected` |
| GPM6A | scz | Brain_Amygdala | eQTL | 0.060 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Amygdala | sQTL | 0.841 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| GPM6A | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.061 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.990 (susie) | 0.188 (0.043) | 1.18e-05 | 0.004 (13) | `coloc_and_smr_heidi_rejected` |
| GPM6A | scz | Brain_Caudate_basal_ganglia | eQTL | 0.044 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Caudate_basal_ganglia | sQTL | 0.991 (susie) | 0.167 (0.037) | 4.97e-06 | 0.036 (8) | `coloc_and_smr_heidi_not_rejected` |
| GPM6A | scz | Brain_Putamen_basal_ganglia | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| GPM6A | scz | Brain_Putamen_basal_ganglia | sQTL | 0.908 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| HMOX2 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.022 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| HMOX2 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.831 (susie) | -0.089 (0.025) | 0.000371 | 0.138 (20) | `coloc_and_smr_heidi_not_rejected` |
| KLC1 | scz | Brain_Cortex | eQTL | 0.000317 (susie) | 0.150 (0.097) | 0.123 | 0.000978 (20) | `neither_tested_null` |
| KLC1 | scz | Brain_Cortex | sQTL | 0.882 (susie) | -0.293 (0.064) | 4.16e-06 | 0.582 (20) | `coloc_and_smr_heidi_not_rejected` |
| NEK4 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.341 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NEK4 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.820 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NT5C2 | scz | Brain_Amygdala | eQTL | 0.225 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Amygdala | sQTL | 0.823 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NT5C2 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.243 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.901 (susie) | 0.155 (0.031) | 6.91e-07 | 0.243 (12) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Cerebellar_Hemisphere | eQTL | 1.84e-08 (susie) | 0.107 (0.051) | 0.035 | 0.001 (20) | `neither_tested_null` |
| NT5C2 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.943 (susie) | 0.113 (0.019) | 5.87e-09 | 0.645 (20) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Cerebellum | eQTL | 1.79e-08 (susie) | 0.079 (0.035) | 0.024 | 0.005 (20) | `neither_tested_null` |
| NT5C2 | scz | Brain_Cerebellum | sQTL | 0.958 (susie) | 0.088 (0.014) | 8.8e-10 | 0.774 (20) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Cortex | eQTL | 0.963 (susie) | 0.029 (0.030) | 0.328 | 0.097 (20) | `coloc_smr_not_significant` |
| NT5C2 | scz | Brain_Cortex | sQTL | 0.906 (susie) | 0.134 (0.029) | 3.03e-06 | 0.552 (17) | `coloc_and_smr_heidi_not_rejected` |
| NT5C2 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.004 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.822 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| NT5C2 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.406 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| NT5C2 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.907 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| PLCB2 | scz | Brain_Hypothalamus | eQTL | 0.148 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PLCB2 | scz | Brain_Hypothalamus | sQTL | 0.975 (abf) | -0.097 (0.023) | 1.79e-05 | 0.013 (11) | `coloc_and_smr_heidi_not_rejected` |
| PLCB2 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.019 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PLCB2 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.826 (susie) | -0.098 (0.022) | 1.14e-05 | 0.061 (10) | `coloc_and_smr_heidi_not_rejected` |
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
| PPIP5K1 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.018 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIP5K1 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.926 (abf) | 0.114 (0.026) | 1.34e-05 | 0.422 (20) | `coloc_and_smr_heidi_not_rejected` |
| PPIP5K1 | scz | Brain_Cerebellum | eQTL | 0.003 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| PPIP5K1 | scz | Brain_Cerebellum | sQTL | 0.869 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RNASEH2C | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.975 (susie) | -0.151 (0.037) | 3.74e-05 | 0.119 (20) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.981 (susie) | 0.091 (0.024) | 0.000171 | 0.253 (17) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Caudate_basal_ganglia | eQTL | 0.946 (susie) | -0.170 (0.041) | 3.32e-05 | 0.022 (20) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Caudate_basal_ganglia | sQTL | 0.967 (susie) | 0.106 (0.029) | 0.000195 | 0.012 (17) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.957 (susie) | -0.204 (0.050) | 4.54e-05 | 0.026 (20) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.895 (susie) | -0.081 (0.023) | 0.000329 | 0.071 (20) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Cerebellum | eQTL | 0.701 (susie) | -0.103 (0.035) | 0.003 | 0.014 (20) | `neither_tested_null` |
| RNASEH2C | scz | Brain_Cerebellum | sQTL | 0.946 (susie) | -0.090 (0.022) | 4.18e-05 | 0.416 (20) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.945 (susie) | -0.146 (0.033) | 1.37e-05 | 0.040 (20) | `coloc_and_smr_heidi_not_rejected` |
| RNASEH2C | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.866 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| RNASEH2C | scz | Brain_Putamen_basal_ganglia | eQTL | 0.836 (abf) | -0.158 (0.040) | 8e-05 | 0.006 (20) | `coloc_and_smr_heidi_rejected` |
| RNASEH2C | scz | Brain_Putamen_basal_ganglia | sQTL | 0.961 (susie) | 0.086 (0.021) | 5.53e-05 | 0.521 (20) | `coloc_and_smr_heidi_not_rejected` |
| SYT5 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.519 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| SYT5 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.904 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Amygdala | eQTL | 0.833 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Amygdala | sQTL | 0.954 (susie) | -0.059 (0.019) | 0.002 | 0.253 (12) | `coloc_smr_not_significant` |
| TMED4 | scz | Brain_Anterior_cingulate_cortex_BA24 | eQTL | 0.167 (susie) | -0.110 (0.035) | 0.001 | 0.032 (14) | `smr_without_coloc` |
| TMED4 | scz | Brain_Anterior_cingulate_cortex_BA24 | sQTL | 0.914 (susie) | 0.083 (0.024) | 0.00055 | 0.328 (15) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Caudate_basal_ganglia | eQTL | 0.240 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Caudate_basal_ganglia | sQTL | 0.955 (susie) | 0.101 (0.027) | 0.000175 | 0.631 (13) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Cerebellar_Hemisphere | eQTL | 0.173 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Cerebellar_Hemisphere | sQTL | 0.963 (susie) | 0.105 (0.027) | 0.000116 | 0.727 (14) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Cerebellum | eQTL | 0.179 (susie) | -0.143 (0.046) | 0.002 | 0.197 (20) | `smr_without_coloc` |
| TMED4 | scz | Brain_Cerebellum | sQTL | 0.961 (susie) | -0.070 (0.018) | 6.41e-05 | 0.597 (10) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Cortex | eQTL | 0.308 (susie) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Cortex | sQTL | 0.952 (susie) | -0.084 (0.023) | 0.000223 | 0.689 (11) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Frontal_Cortex_BA9 | eQTL | 0.714 (susie) | -0.130 (0.039) | 0.000807 | 0.549 (16) | `smr_without_coloc` |
| TMED4 | scz | Brain_Frontal_Cortex_BA9 | sQTL | 0.954 (susie) | -0.090 (0.024) | 0.000144 | 0.208 (15) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Hippocampus | eQTL | 0.913 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Hippocampus | sQTL | 0.954 (susie) | 0.091 (0.025) | 0.000291 | 0.875 (10) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Hypothalamus | eQTL | 0.610 (susie) | -0.162 (0.052) | 0.002 | 0.367 (18) | `smr_without_coloc` |
| TMED4 | scz | Brain_Hypothalamus | sQTL | 0.965 (susie) | -0.069 (0.017) | 3.41e-05 | 0.271 (14) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Nucleus_accumbens_basal_ganglia | eQTL | 0.603 (abf) | — (—) | — | — (—) | `neither_no_instrument` |
| TMED4 | scz | Brain_Nucleus_accumbens_basal_ganglia | sQTL | 0.955 (susie) | -0.070 (0.017) | 4.88e-05 | 0.840 (16) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Putamen_basal_ganglia | eQTL | 0.896 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Putamen_basal_ganglia | sQTL | 0.956 (susie) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Spinal_cord_cervical_c-1 | eQTL | 0.954 (susie) | -0.175 (0.047) | 0.000191 | 0.692 (11) | `coloc_and_smr_heidi_not_rejected` |
| TMED4 | scz | Brain_Spinal_cord_cervical_c-1 | sQTL | 0.922 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |
| TMED4 | scz | Brain_Substantia_nigra | eQTL | 0.002 (susie) | -0.034 (0.031) | 0.264 | 0.653 (5) | `neither_tested_null` |
| TMED4 | scz | Brain_Substantia_nigra | sQTL | 0.914 (abf) | — (—) | — | — (—) | `coloc_no_instrument` |

## The multi-SNP arm

`--smr-multi` combines the cis SNPs surviving LD pruning at r2 0.1 instead of testing the top SNP alone. It is a robustness arm for the SMR estimate, not a second discovery pass: it answers whether a signal rests on one lead SNP or on the cis signal more broadly.

`p_SMR` and `smr_status` in this table are still the single-SNP quantities, so rows match the primary arm one for one. The multi-SNP verdict is `p_SMR_multi` / `smr_multi_status`, Bonferroni-corrected over the probes the multi-SNP test actually ran on (`n_multi_family`) -- a smaller denominator than the instrumented count, because SMR skips the test where too few cis SNPs survive pruning. Those probes are `multi_unavailable`, which is **not** a null result.

| smr_multi_status | probes |
|---|---|
| `no_instrument` | 1,803 |
| `multi_tested_null` | 63 |
| `multi_signal_heidi_unavailable` | 1 |
| `multi_signal_heidi_rejects` | 43 |
| `multi_heidi_supported` | 231 |

Of the 271 probes significant on the single-SNP test, 269 are also significant under the multi-SNP test and 0 could not be tested. A probe that does not survive is a signal carried by its lead SNP alone; that is a caveat on the SMR estimate, not a refutation of the colocalization.

Non-primary sQTL probes (the gene's other introns) are in `smr_results.parquet` with `primary_probe = False`; they are reported so the SMR evidence is not a maximum over introns, and they are not summarized here.

