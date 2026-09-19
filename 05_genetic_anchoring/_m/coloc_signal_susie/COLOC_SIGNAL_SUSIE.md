# Signal-level colocalization (coloc.susie)

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `representative`.

Each gene contributes GTEx's single grouped-permutation representative intron, so a locus whose disease-relevant event is not that intron cannot colocalize here however real it is. The `all` arm is what tests a named event.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 79,245
- cells surviving `gtex_matched`: 9,415 (1,210 genes)
- cells surviving `all_signals`: 11,246 (1,294 genes)
- estimator hierarchy: **9,415 coloc.susie**, 153,489 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_locus_over_max_snps`: 70,191
    - `gwas_no_credible_set`: 49,625
    - `no_qtl_credible_set`: 31,842
    - `qtl_cs_not_matching_gtex`: 1,831
- prior robustness of the 1,097 cells calling at the primary prior: abf `intermediate` 386, abf `primary_prior` 320, abf `robust` 148, susie `intermediate` 141, susie `primary_prior` 56, susie `robust` 46

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 1287 | 19 | 45 | 9 | 35 | 0.000106 |
| aging__als | als | 591 | 17 | 19 | 10 | 12 | 0.832 |
| aging__lbd | lbd | 218 | 1 | 2 | 1 | 2 | 1.000 |
| aging__pd | pd | 536 | 8 | 16 | 4 | 12 | 0.077 |
| aging__scz | scz | 3133 | 65 | 114 | 34 | 83 | 6.78e-06 |
| brainseq-sczd__scz | scz | 579 | 11 | 24 | 4 | 17 | 0.007 |
| POOLED | ALL | 6344 | 121 | 220 | 62 | 161 | 2.49e-11 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| ZDHHC12 | scz | 0.998 | 0.731 | 13/13 | Brain_Spinal_cord_cervical_c-1 |
| FOXN2 | scz | 0.997 | 0.997 | 2/13 | Brain_Cerebellar_Hemisphere |
| ACTR1B | scz | 0.997 | 0.996 | 11/13 | Brain_Caudate_basal_ganglia |
| ACTR1B | scz | 0.997 | 0.996 | 11/13 | Brain_Caudate_basal_ganglia |
| TPP1 | als | 0.996 | 0.774 | 1/9 | Brain_Cerebellum |
| BCKDK | ad | 0.995 | 0.258 | 4/13 | Brain_Spinal_cord_cervical_c-1 |
| CDHR3 | pd | 0.995 | 0.019 | 2/13 | Brain_Cerebellar_Hemisphere |
| MRPS33 | scz | 0.994 | 0.970 | 3/13 | Brain_Hypothalamus |
| IRF3 | scz | 0.994 | 0.977 | 11/13 | Brain_Cortex |
| GPM6A | scz | 0.993 | 0.992 | 4/13 | Brain_Caudate_basal_ganglia |
| MAD1L1 | scz | 0.992 | 0.870 | 3/13 | Brain_Anterior_cingulate_cortex_BA24 |
| RAD51C | ad | 0.992 | 0.031 | 1/13 | Brain_Cortex |
| C9orf72 | als | 0.991 | 0.293 | 8/13 | Brain_Cerebellum |
| SLC39A13 | ad | 0.987 | 0.922 | 2/13 | Brain_Spinal_cord_cervical_c-1 |
| YPEL1 | scz | 0.985 | 0.864 | 8/13 | Brain_Cerebellar_Hemisphere |
| PBRM1 | scz | 0.985 | 0.523 | 1/13 | Brain_Cerebellum |
| POLG | scz | 0.982 | 0.274 | 2/13 | Brain_Cerebellar_Hemisphere |
| PSMD6 | scz | 0.976 | 0.228 | 1/13 | Brain_Hippocampus |
| SNAP91 | scz | 0.974 | 0.930 | 6/13 | Brain_Cerebellum |
| SNCA | lbd | 0.974 | 0.113 | 8/13 | Brain_Cortex |
| RAI1 | scz | 0.971 | 0.717 | 2/13 | Brain_Cerebellum |
| TXNDC15 | als | 0.969 | 0.379 | 1/13 | Brain_Caudate_basal_ganglia |
| DOC2A | scz | 0.966 | 0.755 | 9/13 | Brain_Cortex |
| SIRPA | ad | 0.965 | 0.947 | 13/13 | Brain_Cerebellar_Hemisphere |
| EFHB | scz | 0.965 | 0.000984 | 1/11 | Brain_Frontal_Cortex_BA9 |

