# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Anterior_cingulate_cortex_BA24 sQTLs (gtex-aging/anterior_cingulate_cortex_ba24)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 166 | 0.3031 | 0.3395 | 0.1807 | 0.0154 | 0.988 |
| pheno_sig_modules | 108 | 0.2938 | 0.352 | 0.1759 | 0.0217 | 0.998 |
| go_invisible_modules | 24 | 0.2923 | 0.3968 | 0.1667 | 0.0824 | 0.989 |
| go_visible_modules | 84 | 0.2942 | 0.3387 | 0.1786 | 0.0043 | 0.974 |

## Reading

- mean |rho| above the permutation null would indicate the sQTL and the isoform switch move the *same* transcripts in a coupled direction. Note this within-gene test is underpowered: a single lead sQTL variant often tags introns with near-constant per-transcript direction, so per-cohort |rho| hugs the null; the pooled meta and the diagnosis live in `05_genetic_anchoring/_m/sqtl_concordance_meta/`.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable. The disease-anchored directional test is colocalization (`sqtl_coloc`), not this one.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts.
