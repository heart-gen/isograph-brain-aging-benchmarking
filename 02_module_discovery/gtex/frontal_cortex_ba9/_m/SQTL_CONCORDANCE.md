# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Frontal_Cortex_BA9 sQTLs (gtex-aging/frontal_cortex_ba9)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region frontal_cortex_ba9`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 365 | 0.3475 | 0.3409 | 0.2521 | -0.0157 | 0.273 |
| pheno_sig_modules | 277 | 0.355 | 0.3509 | 0.2635 | -0.0168 | 0.369 |
| go_invisible_modules | 49 | 0.3353 | 0.2932 | 0.1633 | -0.0241 | 0.0695 |
| go_visible_modules | 228 | 0.3592 | 0.3636 | 0.2851 | -0.0152 | 0.615 |

## Reading

- mean |rho| above the permutation null would indicate the sQTL and the isoform switch move the *same* transcripts in a coupled direction. Note this within-gene test is underpowered: a single lead sQTL variant often tags introns with near-constant per-transcript direction, so per-cohort |rho| hugs the null; the pooled meta and the diagnosis live in `05_genetic_anchoring/_m/sqtl_concordance_meta/`.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable. The disease-anchored directional test is colocalization (`sqtl_coloc`), not this one.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts.
