# Signal-level colocalization (coloc.susie)

**SENSITIVITY ARM** (`/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/05_genetic_anchoring/_m/coloc_signal_susie/sensitivity/max_snps_30000`): the GWAS SuSiE SNP guard was raised above the primary 12,000 for the analyses in this directory only. Loci recovered this way were empirically enriched for GWAS-reference-LD inconsistency, and none may enter the biological narrative without a locus-specific LD audit (`08d.locus_ld_robustness.sh`). This is not the primary grid.

Gene pool arm: `switch`. Primary signal filter: `gtex_matched`. sQTL phenotypes: `representative`.

Each gene contributes GTEx's single grouped-permutation representative intron, so a locus whose disease-relevant event is not that intron cannot colocalize here however real it is. The `all` arm is what tests a named event.

`coloc.susie` fine-maps both traits and colocalizes credible set against credible set, so a locus carrying more than one causal signal is not forced into the single-causal-variant assumption `coloc.abf` makes. The QTL side is re-fit with `susie_rss` because GTEx v11 ships credible-set summaries, not SuSiE objects.

## What was fit

- signal-pair posteriors: 7,790
- cells surviving `gtex_matched`: 770 (150 genes)
- cells surviving `all_signals`: 1,006 (164 genes)
- estimator hierarchy: **770 coloc.susie**, 9,540 coloc.abf fallback
  - why coloc.abf, first place the cell left the signal-level pipeline:
    - `gwas_no_credible_set`: 5,118
    - `no_qtl_credible_set`: 3,562
    - `gwas_locus_over_max_snps`: 624
    - `qtl_cs_not_matching_gtex`: 236
- prior robustness of the 62 cells calling at the primary prior: abf `intermediate` 25, abf `primary_prior` 9, abf `robust` 6, susie `intermediate` 1, susie `primary_prior` 3, susie `robust` 18

## Reference-LD caveat

The QTL side is fine-mapped under the 1000G EUR Phase 3 panel, because that is the only LD available for GTEx. GTEx brain donors are not a pure EUR sample, so a re-fit credible set can be a reference-LD artefact. The primary arm therefore keeps only QTL signals that share a variant with **GTEx's own** shipped credible set (`cs_matches_gtex`); the `all_signals` arm drops that requirement and is reported as a sensitivity, never as the headline. `kriging_rss` diagnostics are in `signal_pairs.parquet`.

## Paired modality contrast

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| aging__ad | ad | 403 | 4 | 13 | 2 | 11 | 0.022 |
| POOLED | ALL | 403 | 4 | 13 | 2 | 11 | 0.022 |

Conditioning on `PP4_sQTL >= 0.8` and then reading `PP4_eQTL` is a **selection**, so a splicing-preferential count is a set of locus nominations, not an unbiased splicing-specificity estimate. The unbiased test is the paired McNemar / Wilcoxon above.

## Per-tissue consistency

`genes.parquet` carries `n_tissue_sQTL_coloc`, `frac_tissue_sQTL_coloc`, `max_tissue_sQTL` and the full `tissue_pp4_sQTL` vector. The headline PP4 is a maximum over 13 correlated tissues; quote it with the consistency count beside it, never alone.

| gene | trait | PP4 sQTL | PP4 eQTL | tissues coloc | max tissue |
|---|---|---|---|---|---|
| SIRPA | ad | 0.965 | 0.947 | 13/13 | Brain_Cerebellar_Hemisphere |
| ZNF232 | ad | 0.899 | 0.750 | 1/13 | Brain_Spinal_cord_cervical_c-1 |
| NDUFS3 | ad | 0.864 | 0.099 | 1/13 | Brain_Cerebellar_Hemisphere |
| TPCN1 | ad | 0.829 | 0.819 | 1/13 | Brain_Cerebellum |
| ZMAT2 | ad | 0.627 | 0.017 | 0/13 | Brain_Substantia_nigra |
| MTMR4 | ad | 0.622 | 0.205 | 0/13 | Brain_Frontal_Cortex_BA9 |
| CLU | ad | 0.606 | 0.989 | 0/13 | Brain_Putamen_basal_ganglia |
| STX4 | ad | 0.602 | 0.394 | 0/13 | Brain_Nucleus_accumbens_basal_ganglia |
| TMEM163 | ad | 0.600 | 0.855 | 0/13 | Brain_Cortex |
| BCAR1 | ad | 0.507 | 0.031 | 0/13 | Brain_Cortex |
| ALKBH5 | ad | 0.421 | 0.512 | 0/13 | Brain_Spinal_cord_cervical_c-1 |
| PLEKHM1 | ad | 0.415 | 0.385 | 0/13 | Brain_Cerebellar_Hemisphere |
| FMNL1 | ad | 0.398 | 0.359 | 0/13 | Brain_Anterior_cingulate_cortex_BA24 |
| MAPT | ad | 0.381 | 0.358 | 0/13 | Brain_Cerebellum |
| UBTF | ad | 0.381 | 0.096 | 0/13 | Brain_Cortex |
| NDUFA11 | ad | 0.373 | 0.040 | 0/13 | Brain_Anterior_cingulate_cortex_BA24 |
| YWHAZ | ad | 0.355 | 0.003 | 0/13 | Brain_Amygdala |
| SCN3B | ad | 0.349 | 0.014 | 0/9 | Brain_Anterior_cingulate_cortex_BA24 |
| EFCAB7 | ad | 0.333 | 0.017 | 0/13 | Brain_Putamen_basal_ganglia |
| PRRT2 | ad | 0.314 | 0.133 | 0/13 | Brain_Nucleus_accumbens_basal_ganglia |
| NMT1 | ad | 0.297 | 0.00062 | 0/13 | Brain_Hippocampus |
| PACSIN1 | ad | 0.293 | 0.068 | 0/13 | Brain_Hippocampus |
| TOP3A | ad | 0.280 | 0.421 | 0/13 | Brain_Substantia_nigra |
| MNT | ad | 0.277 | 0.063 | 0/13 | Brain_Cerebellar_Hemisphere |
| SLC22A23 | ad | 0.270 | 0.080 | 0/13 | Brain_Frontal_Cortex_BA9 |

