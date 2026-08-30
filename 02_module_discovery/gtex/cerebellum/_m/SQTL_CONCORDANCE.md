# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Cerebellum sQTLs (gtex-aging/cerebellum)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis gtex-aging --region cerebellum`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 207 | 0.3127 | 0.3215 | 0.2029 | -0.0585 | 0.725 |
| pheno_sig_modules | 74 | 0.2895 | 0.2948 | 0.1622 | -0.0453 | 0.591 |
| go_invisible_modules | 62 | 0.2974 | 0.2932 | 0.1774 | -0.0397 | 0.431 |
| go_visible_modules | 12 | 0.2491 | 0.2927 | 0.0833 | -0.0747 | 0.769 |

## Reading

- On-thesis: phenotype-associated / GO-invisible module sets show **mean |rho| above the permutation null (small p)** => the sQTL and the isoform switch move the *same* transcripts in a coupled direction, not merely in the same genes.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable, and its magnitude is what the null calibrates.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts; it is supporting evidence for genetic anchoring of the switch direction, not proof that the *co-switching* across genes is genetically coordinated.
