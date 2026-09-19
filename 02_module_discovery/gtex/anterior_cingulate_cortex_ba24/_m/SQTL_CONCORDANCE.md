# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Anterior_cingulate_cortex_BA24 sQTLs (gtex-aging/anterior_cingulate_cortex_ba24)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 154 | 0.3258 | 0.3404 | 0.2273 | 0.0106 | 0.795 |
| pheno_sig_modules | 93 | 0.2969 | 0.3397 | 0.1828 | -0.005 | 0.974 |
| go_visible_modules | 93 | 0.2969 | 0.3403 | 0.1828 | -0.005 | 0.975 |

## Reading

- mean |rho| above the permutation null would indicate the sQTL and the isoform switch move the *same* transcripts in a coupled direction. Note this within-gene test is underpowered: a single lead sQTL variant often tags introns with near-constant per-transcript direction, so per-cohort |rho| hugs the null; the pooled meta and the diagnosis live in `05_genetic_anchoring/_m/sqtl_concordance_meta/`.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable. The disease-anchored directional test is colocalization (`sqtl_coloc`), not this one.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts.
