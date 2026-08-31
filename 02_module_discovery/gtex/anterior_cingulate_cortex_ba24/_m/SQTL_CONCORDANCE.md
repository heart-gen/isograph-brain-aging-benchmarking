# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Anterior_cingulate_cortex_BA24 sQTLs (gtex-aging/anterior_cingulate_cortex_ba24)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 178 | 0.2976 | 0.3077 | 0.1629 | -0.0507 | 0.75 |
| pheno_sig_modules | 139 | 0.3143 | 0.3219 | 0.1871 | -0.0723 | 0.664 |
| go_invisible_modules | 75 | 0.2697 | 0.3103 | 0.1467 | -0.065 | 0.963 |
| go_visible_modules | 64 | 0.3666 | 0.3373 | 0.2344 | -0.0808 | 0.142 |

## Reading

- On-thesis: phenotype-associated / GO-invisible module sets show **mean |rho| above the permutation null (small p)** => the sQTL and the isoform switch move the *same* transcripts in a coupled direction, not merely in the same genes.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable, and its magnitude is what the null calibrates.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts; it is supporting evidence for genetic anchoring of the switch direction, not proof that the *co-switching* across genes is genetically coordinated.
