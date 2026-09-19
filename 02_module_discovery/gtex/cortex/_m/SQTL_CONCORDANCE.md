# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Cortex sQTLs (gtex-aging/cortex)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region cortex`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 297 | 0.3342 | 0.349 | 0.2424 | 0.0125 | 0.888 |
| pheno_sig_modules | 174 | 0.3228 | 0.337 | 0.2241 | 0.0266 | 0.809 |
| go_invisible_modules | 63 | 0.3137 | 0.321 | 0.2063 | 0.0794 | 0.615 |
| go_visible_modules | 111 | 0.328 | 0.3459 | 0.2342 | -0.0033 | 0.797 |

## Reading

- mean |rho| above the permutation null would indicate the sQTL and the isoform switch move the *same* transcripts in a coupled direction. Note this within-gene test is underpowered: a single lead sQTL variant often tags introns with near-constant per-transcript direction, so per-cohort |rho| hugs the null; the pooled meta and the diagnosis live in `05_genetic_anchoring/_m/sqtl_concordance_meta/`.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable. The disease-anchored directional test is colocalization (`sqtl_coloc`), not this one.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts.
