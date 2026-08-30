# Stratified LD-score regression: genetic anchoring of the switch layer

Partitioned SNP-heritability of the switch layer's GTEx brain sQTL / eQTL annotations on top of baselineLD v2.2 — the size-robust step beyond MAGMA. `coef_p` is the one-sided p that the annotation's per-SNP heritability coefficient (tau) exceeds 0, conditional on baselineLD (single-annotation models) or on baselineLD + the other QTL layer (joint model).

## Disease case (annotation: brainseq-sczd)

| trait | model | annot | prop_snps | enrichment | enrichment_p | coef_z | coef_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scz | sqtl_only | sqtl_switch | 0.00323 | 0.707 | 0.615 | -0.0771 | 0.531 |
| scz | eqtl_only | eqtl_switch | 0.0058 | 2.04 | 0.0364 | 2.15 | 0.0158 |
| scz | cis_only | cis_switch | 0.00745 | 1.6 | 0.164 | 1.44 | 0.0749 |
| scz | joint | sqtl_switch | 0.00323 | 0.537 | 0.398 | -1.48 | 0.931 |
| scz | joint | eqtl_switch | 0.0058 | 2.11 | 0.0302 | 2.49 | 0.0064 |

## Aging case (annotation: aging)

| trait | model | annot | prop_snps | enrichment | enrichment_p | coef_z | coef_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scz | sqtl_only | sqtl_switch | 0.0238 | 1.74 | 0.00513 | 1.57 | 0.0581 |
| scz | eqtl_only | eqtl_switch | 0.054 | 1.85 | 9.63e-07 | 3.38 | 0.000357 |
| scz | cis_only | cis_switch | 0.0639 | 1.75 | 1.46e-06 | 3.13 | 0.000863 |
| scz | joint | sqtl_switch | 0.0238 | 1.65 | 0.0136 | 0.138 | 0.445 |
| scz | joint | eqtl_switch | 0.054 | 1.84 | 1.34e-06 | 2.98 | 0.00146 |
| ad | sqtl_only | sqtl_switch | 0.0238 | 3.37 | 0.00211 | 1.81 | 0.0355 |
| ad | eqtl_only | eqtl_switch | 0.054 | 2.42 | 0.000401 | 1.49 | 0.0683 |
| ad | cis_only | cis_switch | 0.0639 | 2.57 | 0.000116 | 1.84 | 0.0325 |
| ad | joint | sqtl_switch | 0.0238 | 3.35 | 0.00291 | 1.47 | 0.0702 |
| ad | joint | eqtl_switch | 0.054 | 2.31 | 0.000841 | 0.509 | 0.306 |
| pd | sqtl_only | sqtl_switch | 0.0238 | 3.8 | 0.0489 | 1.47 | 0.0708 |
| pd | eqtl_only | eqtl_switch | 0.054 | 2.73 | 0.000502 | 2.61 | 0.00458 |
| pd | cis_only | cis_switch | 0.0639 | 2.79 | 0.00331 | 2 | 0.0226 |
| pd | joint | sqtl_switch | 0.0238 | 3.7 | 0.0666 | 1.07 | 0.142 |
| pd | joint | eqtl_switch | 0.054 | 2.51 | 0.000475 | 1.05 | 0.148 |
| lbd | sqtl_only | sqtl_switch | 0.0238 | 7.17 | 0.0684 | 1.56 | 0.0599 |
| lbd | eqtl_only | eqtl_switch | 0.054 | 5.35 | 0.00747 | 2.2 | 0.014 |
| lbd | cis_only | cis_switch | 0.0639 | 5.73 | 0.0104 | 2.16 | 0.0153 |
| lbd | joint | sqtl_switch | 0.0238 | 6.89 | 0.104 | 0.99 | 0.161 |
| lbd | joint | eqtl_switch | 0.054 | 4.69 | 0.0215 | 0.956 | 0.17 |
| als | sqtl_only | sqtl_switch | 0.0238 | 2.96 | 0.00222 | 1.2 | 0.114 |
| als | eqtl_only | eqtl_switch | 0.054 | 2.62 | 0.00253 | 1.42 | 0.0781 |
| als | cis_only | cis_switch | 0.0639 | 2.48 | 0.00206 | 1.22 | 0.112 |
| als | joint | sqtl_switch | 0.0238 | 2.85 | 0.00565 | 0.423 | 0.336 |
| als | joint | eqtl_switch | 0.054 | 2.58 | 0.00559 | 0.962 | 0.168 |

## Reading

- **SCZ** (disease): switch-layer sQTL enrichment 0.71× (tau p=0.531); eQTL 2.04× (p=0.0158); total cis 1.60× (p=0.0749). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.931, eQTL=0.0064.
- **SCZ** (aging): switch-layer sQTL enrichment 1.74× (tau p=0.0581); eQTL 1.85× (p=0.000357); total cis 1.75× (p=0.000863). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.445, eQTL=0.00146.
- **AD** (aging): switch-layer sQTL enrichment 3.37× (tau p=0.0355); eQTL 2.42× (p=0.0683); total cis 2.57× (p=0.0325). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.0702, eQTL=0.306.
- **PD** (aging): switch-layer sQTL enrichment 3.80× (tau p=0.0708); eQTL 2.73× (p=0.00458); total cis 2.79× (p=0.0226). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.142, eQTL=0.148.
- **LBD** (aging): switch-layer sQTL enrichment 7.17× (tau p=0.0599); eQTL 5.35× (p=0.014); total cis 5.73× (p=0.0153). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.161, eQTL=0.17.
- **ALS** (aging): switch-layer sQTL enrichment 2.96× (tau p=0.114); eQTL 2.62× (p=0.0781); total cis 2.48× (p=0.112). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.336, eQTL=0.168.

S-LDSC is robust to the gene-size confound that inflates MAGMA on giant modules, so a positive coefficient is heritability genuinely concentrated in the switch layer's cis-regulatory variants, not an artifact of annotation size. Because a gene's sQTL and eQTL SNPs overlap, the joint model splits signal between them and understates each; the single-annotation coefficients are the primary enrichment test and the joint model is only the head-to-head contrast. GTEx brain QTLs are bulk-tissue, so cell-type-specific splicing (e.g. microglial for AD, dopaminergic for PD) is under-sampled — enrichment here is a floor, not a ceiling.
