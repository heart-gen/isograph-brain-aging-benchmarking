# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Frontal_Cortex_BA9 sQTLs (gtex-aging/frontal_cortex_ba9)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region frontal_cortex_ba9`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 304 | 0.352 | 0.3587 | 0.2664 | -0.0236 | 0.705 |
| pheno_sig_modules | 235 | 0.3577 | 0.3586 | 0.2723 | -0.0221 | 0.517 |
| go_invisible_modules | 115 | 0.3491 | 0.3524 | 0.2261 | -0.0292 | 0.561 |
| go_visible_modules | 120 | 0.3661 | 0.3636 | 0.3167 | -0.0153 | 0.458 |

## Reading

- mean |rho| above the permutation null would indicate the sQTL and the isoform switch move the *same* transcripts in a coupled direction. Note this within-gene test is underpowered: a single lead sQTL variant often tags introns with near-constant per-transcript direction, so per-cohort |rho| hugs the null; the pooled meta and the diagnosis live in `05_genetic_anchoring/_m/sqtl_concordance_meta/`.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable. The disease-anchored directional test is colocalization (`sqtl_coloc`), not this one.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts.
