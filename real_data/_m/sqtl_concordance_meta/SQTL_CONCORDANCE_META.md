# sQTL direction concordance — pooled across brain cohorts

Per-gene Spearman rho between IsoGraph's switch-axis usage change and the lead sQTL's per-transcript intron direction, pooled across the 17 brain cohorts. `pvalue` compares the pooled mean |rho| to a within-gene rank-permutation null (seed 13, 2000 draws); |rho| is used because both the sQTL allele reference and the switch-axis orientation are arbitrary.

Reproduce: `python -m isograph_benchmark.real_data.sqtl_concordance_meta` (after the per-cohort array in `real_data/brainseq/_h/14.sqtl_concordance.sh`).

## Pooled concordance by module set

| module_set | n_obs | n_cohorts | mean_abs_rho | perm_mean_abs_rho | frac_strong | pvalue |
| --- | --- | --- | --- | --- | --- | --- |
| all_modules | 2185 | 17 | 0.3088 | 0.332 | 0.2137 | 1 |
| pheno_sig_modules | 859 | 10 | 0.306 | 0.324 | 0.2002 | 0.994 |
| go_invisible_modules | 517 | 10 | 0.2861 | 0.3106 | 0.176 | 0.997 |
| go_visible_modules | 342 | 8 | 0.336 | 0.3441 | 0.2368 | 0.764 |

## Reading

- **Outcome: this within-gene rank-concordance test is underpowered by construction and returns a null** (pooled mean |rho| at or below the permutation null in every module set). The cause is diagnosed, not biological: a single lead sQTL variant tags introns that map to nearly the same net direction across a gene's transcripts (~two thirds of tested genes have a constant-sign per-transcript genetic direction), so the within-gene correlation collapses to tie-breaking noise and falls below a full-variance null. Absence of concordance here is therefore not evidence of absence of genetic anchoring.
- The directional question is instead resolved by colocalization (`sqtl_coloc`), where the GWAS supplies a disease-anchored allele direction and a shared-causal-variant posterior, rather than by this allele-reference-free relative test. The positive genetic-anchoring evidence is the sQTL/eQTL specificity enrichment (`qtl_anchoring`) plus that colocalization.
- Scope: introns mapped to the GENCODE cache (~80% of GTEx brain sQTL introns), genes with enough interpretable switch transcripts, and cis-sQTL anchoring of member-gene splicing rather than the cross-gene co-switching itself. Each cohort contributes independent measurements; the same gene in two tissues is kept as two observations.
