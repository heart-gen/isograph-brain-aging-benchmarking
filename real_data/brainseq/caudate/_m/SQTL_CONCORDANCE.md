# sQTL direction concordance — IsoGraph switches vs GTEx Brain_Caudate_basal_ganglia sQTLs (brainseq-aging/caudate)

Per gene with a lead sQTL and >= 3 shared transcripts, the within-gene Spearman correlation (rho) between IsoGraph's switch-axis usage change (mean_high - mean_low) and the sQTL lead variant's net intron slope per transcript. Aggregated as mean |rho| per module set; `pvalue` is a within-gene rank-permutation null (seed 13, 2000 draws).

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance --analysis brainseq-aging --region caudate`

## Concordance by module set

| module_set | n_genes | mean_abs_rho | perm_mean_abs_rho | frac_strong | mean_signed_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 144 | 0.3114 | 0.3916 | 0.2083 | 0.0477 | 1 |
| pheno_sig_modules | 56 | 0.3182 | 0.3743 | 0.25 | 0.0242 | 0.979 |
| go_invisible_modules | 56 | 0.3182 | 0.3747 | 0.25 | 0.0242 | 0.979 |

## Reading

- On-thesis: phenotype-associated / GO-invisible module sets show **mean |rho| above the permutation null (small p)** => the sQTL and the isoform switch move the *same* transcripts in a coupled direction, not merely in the same genes.
- The statistic is flip-invariant by construction: the sQTL allele reference and the module switch-axis orientation are both arbitrary, so only the *relative* within-gene direction structure is testable, and its magnitude is what the null calibrates.
- Scope: concordance is evaluated on genes whose introns map to the GENCODE cache (~80% of GTEx brain sQTL introns) and that carry >= 3 interpretable switch transcripts; it is supporting evidence for genetic anchoring of the switch direction, not proof that the *co-switching* across genes is genetically coordinated.
