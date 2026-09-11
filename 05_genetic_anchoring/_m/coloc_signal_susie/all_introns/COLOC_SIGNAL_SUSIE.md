# Signal-level colocalization (coloc.susie)

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `all`.

Every intron phenotype of each gene is tested, not only GTEx's grouped-permutation representative, so a named literature event is testable. The `coloc.abf` fallback rows still come from the representative-intron layer.

A cell's PP4 here is a **maximum over that gene's introns**, and `n_phenotypes_tested` records how many it was taken over (1 in the representative arm, which is what makes the arms comparable). Read a PP4 against that count, not on its own; `phenotype_id` names the intron that won, which is what the event audit checks against the curated coordinates.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 36,955
- cells surviving `gtex_matched`: 2,933 (435 genes)
- cells surviving `all_signals`: 3,608 (468 genes)
- estimator hierarchy: **2,933 coloc.susie**, 39,197 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_no_credible_set`: 16,043
    - `no_qtl_credible_set`: 11,705
    - `gwas_locus_over_max_snps`: 10,774
    - `qtl_cs_not_matching_gtex`: 675
- prior robustness of the 349 cells calling at the primary prior: abf `intermediate` 121, abf `primary_prior` 70, abf `robust` 24, susie `intermediate` 86, susie `primary_prior` 27, susie `robust` 21

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 348 | 9 | 11 | 4 | 6 | 0.754 |
| aging__als | als | 159 | 9 | 9 | 5 | 5 | 1.000 |
| aging__lbd | lbd | 58 | 1 | 1 | 1 | 1 | 1.000 |
| aging__pd | pd | 155 | 6 | 7 | 4 | 5 | 1.000 |
| aging__scz | scz | 852 | 17 | 33 | 8 | 24 | 0.007 |
| brainseq-sczd__scz | scz | 66 | 0 | 4 | 0 | 4 | 0.125 |
| POOLED | ALL | 1638 | 42 | 65 | 22 | 45 | 0.007 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| TPP1 | als | 0.996 | 0.774 | 1/9 | Brain_Cerebellum |
| GPM6A | scz | 0.991 | 0.990 | 4/13 | Brain_Caudate_basal_ganglia |
| CTSH | ad | 0.989 | 0.993 | 6/13 | Brain_Hippocampus |
| PITPNM2 | pd | 0.989 | 0.614 | 1/13 | Brain_Cerebellar_Hemisphere |
| RNASEH2C | scz | 0.981 | 0.975 | 6/13 | Brain_Anterior_cingulate_cortex_BA24 |
| GABBR2 | scz | 0.976 | 0.085 | 1/11 | Brain_Cerebellum |
| PGS1 | als | 0.976 | 0.165 | 13/13 | Brain_Cortex |
| PLCB2 | scz | 0.975 | 0.157 | 2/13 | Brain_Hypothalamus |
| SNCA | lbd | 0.974 | 0.113 | 11/13 | Brain_Cortex |
| SNCA | pd | 0.974 | 0.382 | 8/13 | Brain_Cerebellum |
| DOC2A | scz | 0.970 | 0.755 | 9/13 | Brain_Frontal_Cortex_BA9 |
| TXNDC15 | als | 0.969 | 0.379 | 1/13 | Brain_Caudate_basal_ganglia |
| SIRPA | ad | 0.965 | 0.947 | 13/13 | Brain_Cerebellar_Hemisphere |
| TMED4 | scz | 0.965 | 0.954 | 13/13 | Brain_Hypothalamus |
| UNC13A | als | 0.961 | 0.377 | 2/13 | Brain_Cerebellum |
| NT5C2 | scz | 0.958 | 0.963 | 7/13 | Brain_Cerebellum |
| CDIP1 | scz | 0.949 | 0.854 | 6/13 | Brain_Caudate_basal_ganglia |
| PPIL2 | scz | 0.943 | 0.938 | 8/13 | Brain_Amygdala |
| PPIP5K1 | scz | 0.926 | 0.761 | 2/13 | Brain_Cerebellar_Hemisphere |
| SPI1 | ad | 0.926 | 0.881 | 1/13 | Brain_Cerebellar_Hemisphere |
| PTPRN | als | 0.918 | 0.964 | 5/13 | Brain_Nucleus_accumbens_basal_ganglia |
| GGNBP2 | als | 0.911 | 0.989 | 1/13 | Brain_Cerebellum |
| SYT5 | scz | 0.904 | 0.911 | 1/13 | Brain_Frontal_Cortex_BA9 |
| ZNF232 | ad | 0.899 | 0.750 | 1/13 | Brain_Spinal_cord_cervical_c-1 |
| TTC19 | pd | 0.898 | 0.794 | 2/13 | Brain_Putamen_basal_ganglia |

