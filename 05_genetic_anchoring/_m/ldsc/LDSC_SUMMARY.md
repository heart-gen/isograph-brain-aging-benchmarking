# Stratified LD-score regression: genetic anchoring of the switch layer

Partitioned SNP-heritability of the switch layer's GTEx brain sQTL / eQTL annotations on top of baselineLD v2.2 — the size-robust step beyond MAGMA. `coef_p` is the one-sided p that the annotation's per-SNP heritability coefficient (tau) exceeds 0, conditional on baselineLD (single-annotation models) or on baselineLD + the other QTL layer (joint model).

## Disease case (annotation: brainseq-sczd)

| trait | model | annot | prop_snps | enrichment | enrichment_p | coef_z | coef_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scz | sqtl_only | sqtl_switch | 0.00473 | 1.34 | 0.388 | -0.405 | 0.657 |
| scz | eqtl_only | eqtl_switch | 0.0105 | 1.53 | 0.0722 | 0.405 | 0.343 |
| scz | cis_only | cis_switch | 0.0127 | 1.47 | 0.0716 | 0.171 | 0.432 |
| scz | joint | sqtl_switch | 0.00473 | 1.31 | 0.427 | -0.751 | 0.774 |
| scz | joint | eqtl_switch | 0.0105 | 1.55 | 0.0602 | 0.69 | 0.245 |

## Aging case (annotation: aging)

| trait | model | annot | prop_snps | enrichment | enrichment_p | coef_z | coef_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scz | sqtl_only | sqtl_switch | 0.0247 | 1.71 | 0.000188 | 1.81 | 0.0353 |
| scz | eqtl_only | eqtl_switch | 0.0561 | 1.78 | 1.67e-07 | 3.01 | 0.00129 |
| scz | cis_only | cis_switch | 0.0661 | 1.73 | 6.08e-08 | 3.03 | 0.00121 |
| scz | joint | sqtl_switch | 0.0247 | 1.65 | 0.000784 | 0.295 | 0.384 |
| scz | joint | eqtl_switch | 0.0561 | 1.78 | 3.64e-07 | 2.52 | 0.00593 |
| ad | sqtl_only | sqtl_switch | 0.0247 | 2.26 | 0.00977 | 0.00313 | 0.499 |
| ad | eqtl_only | eqtl_switch | 0.0561 | 1.75 | 0.00208 | -1.05 | 0.853 |
| ad | cis_only | cis_switch | 0.0661 | 1.84 | 0.00123 | -0.721 | 0.764 |
| ad | joint | sqtl_switch | 0.0247 | 2.32 | 0.00714 | 0.407 | 0.342 |
| ad | joint | eqtl_switch | 0.0561 | 1.73 | 0.00184 | -1.44 | 0.926 |
| pd | sqtl_only | sqtl_switch | 0.0247 | 4.17 | 0.00132 | 2.58 | 0.00495 |
| pd | eqtl_only | eqtl_switch | 0.0561 | 2.93 | 0.000628 | 2.51 | 0.00611 |
| pd | cis_only | cis_switch | 0.0661 | 2.65 | 0.00101 | 2.12 | 0.0171 |
| pd | joint | sqtl_switch | 0.0247 | 4.04 | 0.00138 | 2.23 | 0.0129 |
| pd | joint | eqtl_switch | 0.0561 | 2.71 | 0.000633 | 1.5 | 0.0672 |
| lbd | sqtl_only | sqtl_switch | 0.0247 | 4.16 | 0.164 | 0.791 | 0.214 |
| lbd | eqtl_only | eqtl_switch | 0.0561 | 3.84 | 0.0589 | 1.13 | 0.129 |
| lbd | cis_only | cis_switch | 0.0661 | 3.36 | 0.0966 | 0.855 | 0.196 |
| lbd | joint | sqtl_switch | 0.0247 | 3.89 | 0.215 | 0.334 | 0.369 |
| lbd | joint | eqtl_switch | 0.0561 | 3.71 | 0.0691 | 0.857 | 0.196 |
| als | sqtl_only | sqtl_switch | 0.0247 | 2.75 | 0.00354 | 0.732 | 0.232 |
| als | eqtl_only | eqtl_switch | 0.0561 | 2.01 | 0.0186 | -0.0868 | 0.535 |
| als | cis_only | cis_switch | 0.0661 | 1.98 | 0.0114 | -0.258 | 0.602 |
| als | joint | sqtl_switch | 0.0247 | 2.79 | 0.00298 | 0.894 | 0.186 |
| als | joint | eqtl_switch | 0.0561 | 1.95 | 0.0264 | -0.452 | 0.674 |

## Reading

- **SCZ** (disease): switch-layer sQTL enrichment 1.34× (tau p=0.657); eQTL 1.53× (p=0.343); total cis 1.47× (p=0.432). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.774, eQTL=0.245.
- **SCZ** (aging): switch-layer sQTL enrichment 1.71× (tau p=0.0353); eQTL 1.78× (p=0.00129); total cis 1.73× (p=0.00121). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.384, eQTL=0.00593.
- **AD** (aging): switch-layer sQTL enrichment 2.26× (tau p=0.499); eQTL 1.75× (p=0.853); total cis 1.84× (p=0.764). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.342, eQTL=0.926.
- **PD** (aging): switch-layer sQTL enrichment 4.17× (tau p=0.00495); eQTL 2.93× (p=0.00611); total cis 2.65× (p=0.0171). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.0129, eQTL=0.0672.
- **LBD** (aging): switch-layer sQTL enrichment 4.16× (tau p=0.214); eQTL 3.84× (p=0.129); total cis 3.36× (p=0.196). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.369, eQTL=0.196.
- **ALS** (aging): switch-layer sQTL enrichment 2.75× (tau p=0.232); eQTL 2.01× (p=0.535); total cis 1.98× (p=0.602). Joint sQTL-vs-eQTL conditional tau p: sQTL=0.186, eQTL=0.674.

S-LDSC is robust to the gene-size confound that inflates MAGMA on giant modules, so a positive coefficient is heritability genuinely concentrated in the switch layer's cis-regulatory variants, not an artifact of annotation size. Because a gene's sQTL and eQTL SNPs overlap, the joint model splits signal between them and understates each; the single-annotation coefficients are the primary enrichment test and the joint model is only the head-to-head contrast. GTEx brain QTLs are bulk-tissue, so cell-type-specific splicing (e.g. microglial for AD, dopaminergic for PD) is under-sampled — enrichment here is a floor, not a ceiling.
