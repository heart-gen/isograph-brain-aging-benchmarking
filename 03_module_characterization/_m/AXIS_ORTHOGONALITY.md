# Abundance-vs-switch axis orthogonality, all analyses

Per gene, the Pearson correlation across samples between IsoGraph's abundance channel and its switch channel (`feature_scores.parquet`), computed by `abundance_structure_separation.py --orthogonality-only` for every store and pooled by `--rollup`. Genes with a constant channel are dropped. The manuscript's separability statement quotes the range of the per-analysis medians and of the fraction of genes with |r| < 0.1; it does not quote one analysis.

- analyses: 17 of 17
- genes per analysis: 12,042-13,222
- median |r|: 0.111-0.243 (median of medians 0.189)
- fraction |r| < 0.1: 0.23-0.46
- fraction |r| < 0.3: 0.59-0.91
- fraction |r| > 0.5: 0.011-0.168

| analysis | region | n genes | median \|r\| | IQR | \|r\| < 0.1 | \|r\| < 0.3 | \|r\| > 0.5 | median r |
|---|---|---|---|---|---|---|---|---|
| brainseq-sczd | caudate | 13,222 | 0.124 | 0.057-0.236 | 0.42 | 0.83 | 0.044 | +0.005 |
| brainseq-aging | caudate | 13,177 | 0.130 | 0.060-0.240 | 0.40 | 0.83 | 0.046 | +0.009 |
| brainseq-aging | hippocampus | 13,130 | 0.111 | 0.051-0.197 | 0.46 | 0.91 | 0.011 | +0.011 |
| brainseq-aging | dlpfc | 12,931 | 0.118 | 0.054-0.214 | 0.43 | 0.86 | 0.031 | +0.015 |
| gtex-aging | amygdala | 12,042 | 0.192 | 0.090-0.334 | 0.28 | 0.70 | 0.087 | +0.009 |
| gtex-aging | anterior_cingulate_cortex_ba24 | 12,190 | 0.203 | 0.094-0.365 | 0.27 | 0.67 | 0.120 | +0.018 |
| gtex-aging | caudate_basal_ganglia | 12,424 | 0.214 | 0.097-0.371 | 0.26 | 0.65 | 0.122 | +0.008 |
| gtex-aging | cerebellar_hemisphere | 12,416 | 0.204 | 0.091-0.360 | 0.27 | 0.67 | 0.116 | -0.001 |
| gtex-aging | cerebellum | 12,480 | 0.148 | 0.069-0.279 | 0.35 | 0.78 | 0.065 | +0.007 |
| gtex-aging | cortex | 12,426 | 0.161 | 0.072-0.306 | 0.34 | 0.74 | 0.092 | +0.010 |
| gtex-aging | frontal_cortex_ba9 | 12,322 | 0.190 | 0.085-0.348 | 0.29 | 0.68 | 0.111 | +0.008 |
| gtex-aging | hippocampus | 12,313 | 0.189 | 0.090-0.323 | 0.28 | 0.72 | 0.067 | +0.009 |
| gtex-aging | hypothalamus | 12,648 | 0.212 | 0.100-0.368 | 0.25 | 0.66 | 0.118 | +0.010 |
| gtex-aging | nucleus_accumbens_basal_ganglia | 12,492 | 0.243 | 0.112-0.420 | 0.23 | 0.59 | 0.168 | +0.007 |
| gtex-aging | putamen_basal_ganglia | 12,100 | 0.210 | 0.096-0.372 | 0.26 | 0.65 | 0.120 | +0.009 |
| gtex-aging | spinal_cord_cervical_c_1 | 12,305 | 0.178 | 0.083-0.310 | 0.30 | 0.74 | 0.059 | +0.003 |
| gtex-aging | substantia_nigra | 12,265 | 0.173 | 0.081-0.304 | 0.30 | 0.74 | 0.056 | +0.007 |

Outputs: `axis_orthogonality_all.parquet` (per gene per analysis), `axis_orthogonality_summary.{parquet,csv}` (this table). Figure: `manuscript/_h/abundance_structure_figure.R` panel A (S-real-4).
