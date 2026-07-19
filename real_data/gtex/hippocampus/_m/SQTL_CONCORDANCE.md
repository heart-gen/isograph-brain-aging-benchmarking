# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Hippocampus sQTLs (gtex-aging/hippocampus)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region hippocampus`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 149 | 0.3026 | 0.2866 | 0.2081 | -0.0307 | 0.165 |
| pheno_sig_modules | 106 | 0.2962 | 0.2809 | 0.1792 | -0.0334 | 0.212 |
| go_invisible_modules | 53 | 0.2783 | 0.2881 | 0.1509 | 0.0147 | 0.638 |
| go_visible_modules | 53 | 0.3141 | 0.2736 | 0.2075 | -0.0814 | 0.074 |

## Reading

- On-thesis: phenotype-associated / GO-invisible module sets show **mean |rho| above the permutation null (small p)** => the sQTL and the isoform switch move the *same* transcripts in a coupled direction, not merely in the same genes.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable, and its magnitude is what the null calibrates.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts; it is supporting evidence for genetic anchoring of the switch direction, not proof that the *co-switching* across genes is genetically coordinated.
