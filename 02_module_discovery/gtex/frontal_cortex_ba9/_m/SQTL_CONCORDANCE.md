# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Frontal_Cortex_BA9 sQTLs (gtex-aging/frontal_cortex_ba9)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region frontal_cortex_ba9`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 355 | 0.3198 | 0.3421 | 0.2366 | -0.0272 | 0.973 |
| pheno_sig_modules | 322 | 0.3247 | 0.3454 | 0.2453 | -0.0217 | 0.958 |
| go_invisible_modules | 172 | 0.2882 | 0.2987 | 0.186 | 0.0241 | 0.765 |
| go_visible_modules | 150 | 0.3666 | 0.3971 | 0.3133 | -0.0742 | 0.942 |

## Reading

- On-thesis: phenotype-associated / GO-invisible module sets show **mean |rho| above the permutation null (small p)** => the sQTL and the isoform switch move the *same* transcripts in a coupled direction, not merely in the same genes.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable, and its magnitude is what the null calibrates.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts; it is supporting evidence for genetic anchoring of the switch direction, not proof that the *co-switching* across genes is genetically coordinated.
