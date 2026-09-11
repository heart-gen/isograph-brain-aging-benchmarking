# Signal-level colocalization (coloc.susie)

**SENSITIVITY ARM** (`/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/05_genetic_anchoring/_m/coloc_signal_susie/sensitivity/max_snps_30000`): the GWAS SuSiE SNP guard was raised above the primary 12,000 for the analyses in this directory only. Loci recovered this way were empirically enriched for GWAS-reference-LD inconsistency, and none may enter the biological narrative without a locus-specific LD audit (`26.locus_ld_robustness.sh`). This is not the primary grid.

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `representative`.

Each gene contributes GTEx's single grouped-permutation representative intron, so a locus whose disease-relevant event is not that intron cannot colocalize here however real it is. The `all` arm is what tests a named event.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 10,310
- cells surviving `gtex_matched`: 889 (160 genes)
- cells surviving `all_signals`: 1,136 (175 genes)
- estimator hierarchy: **889 coloc.susie**, 8,004 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_no_credible_set`: 4,153
    - `no_qtl_credible_set`: 3,604
    - `qtl_cs_not_matching_gtex`: 247
- prior robustness of the 73 cells calling at the primary prior: abf `intermediate` 31, abf `primary_prior` 9, abf `robust` 15, susie `intermediate` 4, susie `primary_prior` 5, susie `robust` 9

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 346 | 7 | 11 | 2 | 6 | 0.289 |
| POOLED | ALL | 346 | 7 | 11 | 2 | 6 | 0.289 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| CTSH | ad | 0.989 | 0.993 | 6/13 | Brain_Hippocampus |
| SIRPA | ad | 0.965 | 0.947 | 13/13 | Brain_Cerebellar_Hemisphere |
| SPI1 | ad | 0.926 | 0.881 | 1/13 | Brain_Cerebellar_Hemisphere |
| DOC2A | ad | 0.881 | 0.986 | 1/13 | Brain_Amygdala |
| NDUFS3 | ad | 0.864 | 0.099 | 1/13 | Brain_Cerebellar_Hemisphere |
| TPCN1 | ad | 0.829 | 0.819 | 1/13 | Brain_Cerebellum |
| PICALM | ad | 0.813 | 0.491 | 1/13 | Brain_Cortex |
| VSTM2A | ad | 0.792 | 0.044 | 0/11 | Brain_Substantia_nigra |
| ZNF232 | ad | 0.771 | 0.788 | 0/13 | Brain_Spinal_cord_cervical_c-1 |
| ADAM17 | ad | 0.728 | 0.095 | 0/13 | Brain_Cerebellum |
| NDUFA2 | ad | 0.673 | 0.033 | 0/13 | Brain_Substantia_nigra |
| BIN1 | ad | 0.620 | 0.261 | 0/13 | Brain_Amygdala |
| STX4 | ad | 0.610 | 0.444 | 0/13 | Brain_Nucleus_accumbens_basal_ganglia |
| CLU | ad | 0.604 | 0.989 | 0/13 | Brain_Putamen_basal_ganglia |
| CARF | ad | 0.589 | 0.823 | 0/13 | Brain_Nucleus_accumbens_basal_ganglia |
| FERMT3 | ad | 0.582 | 0.061 | 0/11 | Brain_Putamen_basal_ganglia |
| TMEM106B | ad | 0.531 | 0.923 | 0/13 | Brain_Cerebellar_Hemisphere |
| BCAR1 | ad | 0.507 | 0.031 | 0/13 | Brain_Cortex |
| FCER1G | ad | 0.493 | 0.983 | 0/13 | Brain_Hippocampus |
| SLC43A2 | ad | 0.421 | 0.067 | 0/13 | Brain_Putamen_basal_ganglia |
| PLEKHM1 | ad | 0.413 | 0.384 | 0/13 | Brain_Cerebellar_Hemisphere |
| MKNK2 | ad | 0.402 | 0.084 | 0/13 | Brain_Substantia_nigra |
| FMNL1 | ad | 0.398 | 0.358 | 0/13 | Brain_Anterior_cingulate_cortex_BA24 |
| UBTF | ad | 0.381 | 0.096 | 0/13 | Brain_Cortex |
| MAPT | ad | 0.381 | 0.358 | 0/13 | Brain_Cerebellum |

