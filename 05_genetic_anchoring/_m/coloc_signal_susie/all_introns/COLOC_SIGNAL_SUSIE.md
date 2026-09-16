# Signal-level colocalization (coloc.susie)

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `all`.

Every intron phenotype of each gene is tested, not only GTEx's grouped-permutation representative, so a named literature event is testable. The `coloc.abf` fallback rows still come from the representative-intron layer.

A cell's PP4 here is a **maximum over that gene's introns**, and `n_phenotypes_tested` records how many it was taken over (1 in the representative arm, which is what makes the arms comparable). Read a PP4 against that count, not on its own; `phenotype_id` names the intron that won, which is what the event audit checks against the curated coordinates.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 52,775
- cells surviving `gtex_matched`: 3,830 (541 genes)
- cells surviving `all_signals`: 4,714 (596 genes)
- estimator hierarchy: **3,830 coloc.susie**, 48,054 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_no_credible_set`: 19,956
    - `no_qtl_credible_set`: 14,506
    - `gwas_locus_over_max_snps`: 12,708
    - `qtl_cs_not_matching_gtex`: 884
- prior robustness of the 407 cells calling at the primary prior: abf `intermediate` 141, abf `primary_prior` 98, abf `robust` 64, susie `intermediate` 59, susie `primary_prior` 19, susie `robust` 26

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 403 | 4 | 13 | 2 | 11 | 0.022 |
| aging__als | als | 188 | 6 | 4 | 4 | 2 | 0.688 |
| aging__lbd | lbd | 78 | 0 | 1 | 0 | 1 | 1.000 |
| aging__pd | pd | 167 | 5 | 11 | 1 | 7 | 0.070 |
| aging__scz | scz | 1040 | 24 | 41 | 8 | 25 | 0.005 |
| brainseq-sczd__scz | scz | 152 | 3 | 6 | 1 | 4 | 0.375 |
| POOLED | ALL | 2028 | 42 | 76 | 16 | 50 | 3.33e-05 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| ACTR1B | scz | 0.997 | 0.996 | 11/13 | Brain_Caudate_basal_ganglia |
| ACTR1B | scz | 0.997 | 0.996 | 11/13 | Brain_Caudate_basal_ganglia |
| MRPS33 | scz | 0.994 | 0.970 | 3/13 | Brain_Hypothalamus |
| IRF3 | scz | 0.994 | 0.977 | 11/13 | Brain_Cortex |
| YPEL1 | scz | 0.985 | 0.864 | 8/13 | Brain_Cerebellar_Hemisphere |
| POLG | scz | 0.982 | 0.274 | 2/13 | Brain_Cerebellar_Hemisphere |
| PSMD6 | scz | 0.976 | 0.228 | 1/13 | Brain_Hippocampus |
| GABBR2 | scz | 0.976 | 0.085 | 1/11 | Brain_Cerebellum |
| SNAP91 | scz | 0.974 | 0.930 | 7/13 | Brain_Cerebellum |
| RAI1 | scz | 0.971 | 0.717 | 2/13 | Brain_Cerebellum |
| TXNDC15 | als | 0.969 | 0.379 | 1/13 | Brain_Caudate_basal_ganglia |
| SIRPA | ad | 0.965 | 0.947 | 13/13 | Brain_Cerebellar_Hemisphere |
| PRDM2 | als | 0.964 | 0.495 | 6/13 | Brain_Cortex |
| NUP50 | scz | 0.960 | 0.972 | 11/13 | Brain_Caudate_basal_ganglia |
| UNC13A | als | 0.959 | 0.377 | 2/13 | Brain_Cerebellum |
| FGFR1 | scz | 0.957 | 0.653 | 3/13 | Brain_Cerebellar_Hemisphere |
| DGKZ | scz | 0.953 | 0.858 | 9/13 | Brain_Frontal_Cortex_BA9 |
| CDIP1 | scz | 0.951 | 0.848 | 6/13 | Brain_Caudate_basal_ganglia |
| RASA1 | scz | 0.946 | 0.856 | 1/13 | Brain_Cerebellum |
| PPIL2 | scz | 0.943 | 0.938 | 8/13 | Brain_Amygdala |
| PTPRN | als | 0.918 | 0.964 | 5/13 | Brain_Nucleus_accumbens_basal_ganglia |
| GGNBP2 | als | 0.912 | 0.988 | 1/13 | Brain_Cerebellum |
| NDUFAF7 | scz | 0.909 | 0.980 | 2/13 | Brain_Putamen_basal_ganglia |
| SYT5 | scz | 0.904 | 0.911 | 1/13 | Brain_Frontal_Cortex_BA9 |
| ZNF232 | ad | 0.899 | 0.750 | 1/13 | Brain_Spinal_cord_cervical_c-1 |

