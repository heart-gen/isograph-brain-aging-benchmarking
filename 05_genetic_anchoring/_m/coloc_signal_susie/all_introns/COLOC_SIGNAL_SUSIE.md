# Signal-level colocalization (coloc.susie)

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `all`.

Every intron phenotype of each gene is tested, not only GTEx's grouped-permutation representative, so a named literature event is testable. The `coloc.abf` fallback rows still come from the representative-intron layer.

A cell's PP4 here is a **maximum over that gene's introns**, and `n_phenotypes_tested` records how many it was taken over (1 in the representative arm, which is what makes the arms comparable). Read a PP4 against that count, not on its own; `phenotype_id` names the intron that won, which is what the event audit checks against the curated coordinates.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 123,200
- cells surviving `gtex_matched`: 9,598 (1,214 genes)
- cells surviving `all_signals`: 11,445 (1,301 genes)
- estimator hierarchy: **9,598 coloc.susie**, 153,306 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_locus_over_max_snps`: 70,191
    - `gwas_no_credible_set`: 49,625
    - `no_qtl_credible_set`: 31,643
    - `qtl_cs_not_matching_gtex`: 1,847
- prior robustness of the 1,136 cells calling at the primary prior: abf `intermediate` 383, abf `primary_prior` 316, abf `robust` 148, susie `intermediate` 171, susie `primary_prior` 66, susie `robust` 52

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 1287 | 21 | 45 | 11 | 35 | 0.000536 |
| aging__als | als | 590 | 18 | 19 | 11 | 12 | 1.000 |
| aging__lbd | lbd | 218 | 1 | 2 | 1 | 2 | 1.000 |
| aging__pd | pd | 536 | 9 | 16 | 4 | 11 | 0.118 |
| aging__scz | scz | 3132 | 67 | 114 | 36 | 83 | 1.96e-05 |
| brainseq-sczd__scz | scz | 579 | 11 | 24 | 4 | 17 | 0.007 |
| POOLED | ALL | 6342 | 127 | 220 | 67 | 160 | 5.91e-10 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| TMEM175 | als | 0.999 | 0.139 | 3/13 | Brain_Cerebellar_Hemisphere |
| PILRB | ad | 0.999 | 2.8e-14 | 1/13 | Brain_Putamen_basal_ganglia |
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
| ITGB1BP1 | ad | 0.979 | 0.215 | 5/13 | Brain_Cerebellar_Hemisphere |
| PSMD6 | scz | 0.976 | 0.228 | 1/13 | Brain_Hippocampus |
| GABBR2 | scz | 0.975 | 0.085 | 1/11 | Brain_Cerebellum |
| SNAP91 | scz | 0.974 | 0.930 | 7/13 | Brain_Cerebellum |
| SNCA | lbd | 0.974 | 0.113 | 11/13 | Brain_Cortex |
| RAI1 | scz | 0.971 | 0.717 | 2/13 | Brain_Cerebellum |

