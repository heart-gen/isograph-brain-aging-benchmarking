# Signal-level colocalization (coloc.susie)

**SENSITIVITY ARM** (`/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/05_genetic_anchoring/_m/coloc_signal_susie/sensitivity/max_snps_30000`): the GWAS SuSiE SNP guard was raised above the primary 12,000 for the analyses in this directory only. Loci recovered this way were empirically enriched for GWAS-reference-LD inconsistency, and none may enter the biological narrative without a locus-specific LD audit (`08d.locus_ld_robustness.sh`). This is not the primary grid.

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `representative`.

Each gene contributes GTEx's single grouped-permutation representative intron, so a locus whose disease-relevant event is not that intron cannot colocalize here however real it is. The `all` arm is what tests a named event.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 36,925
- cells surviving `gtex_matched`: 3,282 (525 genes)
- cells surviving `all_signals`: 4,021 (561 genes)
- estimator hierarchy: **3,282 coloc.susie**, 29,626 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_no_credible_set`: 15,523
    - `no_qtl_credible_set`: 11,634
    - `gwas_locus_over_max_snps`: 1,730
    - `qtl_cs_not_matching_gtex`: 739
- prior robustness of the 205 cells calling at the primary prior: abf `intermediate` 55, abf `primary_prior` 49, abf `robust` 31, susie `intermediate` 11, susie `primary_prior` 22, susie `robust` 37

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 1278 | 18 | 45 | 8 | 35 | 4.19e-05 |
| POOLED | ALL | 1278 | 18 | 45 | 8 | 35 | 4.19e-05 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| BCKDK | ad | 0.995 | 0.258 | 4/13 | Brain_Spinal_cord_cervical_c-1 |
| RAD51C | ad | 0.992 | 0.031 | 1/13 | Brain_Cortex |
| SLC39A13 | ad | 0.986 | 0.974 | 2/13 | Brain_Spinal_cord_cervical_c-1 |
| SIRPA | ad | 0.965 | 0.947 | 13/13 | Brain_Cerebellar_Hemisphere |
| YPEL3 | ad | 0.935 | 0.932 | 4/13 | Brain_Substantia_nigra |
| INTS8 | ad | 0.931 | 0.086 | 3/13 | Brain_Cerebellum |
| SPI1 | ad | 0.926 | 0.881 | 1/13 | Brain_Cerebellar_Hemisphere |
| SERPINB1 | ad | 0.917 | 0.967 | 1/13 | Brain_Spinal_cord_cervical_c-1 |
| ZNF232 | ad | 0.899 | 0.750 | 1/13 | Brain_Spinal_cord_cervical_c-1 |
| DOC2A | ad | 0.881 | 0.986 | 1/13 | Brain_Amygdala |
| IFNAR2 | ad | 0.870 | 0.829 | 7/13 | Brain_Nucleus_accumbens_basal_ganglia |
| NDUFS3 | ad | 0.864 | 0.099 | 1/13 | Brain_Cerebellar_Hemisphere |
| VWA5B2 | ad | 0.864 | 0.642 | 1/13 | Brain_Amygdala |
| AKT1 | ad | 0.840 | 0.157 | 1/13 | Brain_Caudate_basal_ganglia |
| COG7 | ad | 0.837 | 0.939 | 1/13 | Brain_Amygdala |
| INO80E | ad | 0.832 | 0.860 | 1/13 | Brain_Nucleus_accumbens_basal_ganglia |
| TPCN1 | ad | 0.830 | 0.819 | 1/13 | Brain_Cerebellum |
| PICALM | ad | 0.813 | 0.491 | 1/13 | Brain_Cortex |
| VSTM2A | ad | 0.792 | 0.044 | 0/11 | Brain_Substantia_nigra |
| PLEKHA1 | ad | 0.778 | 0.943 | 0/13 | Brain_Spinal_cord_cervical_c-1 |
| KANSL1 | ad | 0.771 | 0.434 | 0/13 | Brain_Putamen_basal_ganglia |
| TMEM219 | ad | 0.726 | 0.847 | 0/13 | Brain_Cortex |
| REEP6 | ad | 0.675 | 0.012 | 0/13 | Brain_Cortex |
| NDUFA2 | ad | 0.673 | 0.033 | 0/13 | Brain_Substantia_nigra |
| VEGFB | ad | 0.673 | 0.017 | 0/5 | Brain_Hippocampus |

