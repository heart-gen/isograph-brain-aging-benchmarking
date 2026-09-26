# Open issues found during manuscript review (2026-09-22)

## 1. ~~Module-recovery axes are labelled "AUC" but the metric is best-match Jaccard~~ DONE 2026-09-23

Both scripts relabelled to "best-match Jaccard" (nine axis strings) and Fig 1 plus the
synthetic supplement regenerated. Kept for the record:

`metrics_module_recovery` is the best-match Jaccard between each planted module and
its closest inferred module (`isograph_benchmark/benchmark/partition_metrics.py`,
`01_synthetic_benchmark/03_metrics/_m/PARTITION_METRICS.md`). It is not an AUC. Every
recovery axis nevertheless reads "Module recovery (AUC; 1 = perfect)":

- `isograph_benchmark/figures/synthetic_benchmark.R`: lines 221, 725, 1138, 1205, 1770,
  1786, and the null panel at line 1236 ("False module recovery (AUC; 0 = ideal)")
- `manuscript/_h/concept_overview_figure.R`: lines 229 (Fig 1D) and 291 (Fig 1F)

Affected outputs: Fig 1, `03_metrics/figures/figS*` (including the specificity-null panel),
and the re-laid supplementary figures in the manuscript repo (confound robustness,
extended benchmark, non-switching background, unequal abundance, degradation fallback).

Fix: relabel as "Module recovery\n(best-match Jaccard; 1 = perfect)" and regenerate the
figures. The manuscript text and the S1 caption already call it best-match Jaccard.

## 2. ~~Gene-level overlap of switch-unique sets before/after composition adjustment is not in the pipeline~~ DONE 2026-09-26

`_contrast_rows` now derives `n_overlap`/`n_new` from the base and adjusted `gene_level.parquet` for both the BrainSEQ and GTEx arms, `retained_frac` is gone, and the rollup, `tableS13_composition_adjustment.csv` and `figCompositionRobustness` were regenerated. All twelve rows of the table below reproduce exactly. Kept for the record:

The manuscript's aging-composition paragraph now reports how many unadjusted
switch-unique genes stay switch-unique after cell-type adjustment. For example,
BrainSEQ DLPFC keeps 8 of 22, and 20 of the 28 adjusted genes are new. It no longer
reports `retained_frac`, which divides counts and is misleading: the adjusted set is not a
subset of the unadjusted one.

These overlaps were computed **ad hoc** (an inline script, not a committed CLI). That
breaks the rule in `ANALYSIS_MAP.md` that everything reaching the paper comes from a
committed CLI. The calculation:

```python
b = pd.read_parquet(f"{store}/isograph_vae/incremental_association/gene_level.parquet")
a = pd.read_parquet(f"{store}/isograph_vae/incremental_association_composition/gene_level.parquet")
B = set(b.gene_id[b.category == "composition_unique"])
A = set(a.gene_id[a.category == "composition_unique"])
overlap, new = len(B & A), len(A - B)
```

Values cited in the manuscript (base / adj / overlap / new):

| store | base | adj | overlap | new |
|---|---|---|---|---|
| brainseq/caudate | 45 | 14 | 10 | 4 |
| brainseq/hippocampus | 1 | 1 | 1 | 0 |
| brainseq/dlpfc | 22 | 28 | 8 | 20 |
| brainseq/caudate_sczd | 60 | 11 | 5 | 6 |
| gtex/caudate_basal_ganglia | 87 | 50 | 19 | 31 |
| gtex/hippocampus | 44 | 29 | 21 | 8 |
| gtex/anterior_cingulate_cortex_ba24 | 928 | 365 | 293 | 72 |
| gtex/frontal_cortex_ba9 | 797 | 2 | 1 | 1 |
| gtex/cortex | 1169 | 2 | 0 | 2 |
| gtex/amygdala | 33 | 5 | 4 | 1 |
| gtex/putamen_basal_ganglia | 2 | 4 | 2 | 2 |
| gtex/nucleus_accumbens_basal_ganglia | 2 | 2 | 1 | 1 |

Fix: add `n_overlap` and `n_new` columns to the meta rollup in
`isograph_benchmark/real_data/celltype_composition.py` (`_contrast_rows`, and the GTEx
equivalent). Drop or rename `retained_frac`, re-run `03_module_characterization/_h/02c`,
propagate the result to `tableS13_composition_adjustment.csv`, and confirm the numbers
above reproduce.

Done by: `_switch_unique_genes()` + the two new columns in `_contrast_rows`;
`composition_robustness_figure.R` panel B replotted as `n_overlap / comp_unique_base`
(a real fraction, unlike a ratio of the two counts) annotated with `n_new`; the narrative
in `COMPOSITION_ADJUSTMENT_SUMMARY.md` and `GTEX_COMPOSITION_SUMMARY.md` re-derived from
the rollup rather than from the superseded hard-coded counts.

## 3. The switch–abundance separation claim is general, but only one analysis is behind it

Results subsection 2 says the switch and abundance channels capture partly distinct
within-gene variation and quotes `median |r| = 0.12`, `42% of genes |r| < 0.1`
(Fig 2b). Those numbers come from a single store:
`02_module_discovery/brainseq/caudate_sczd/_m/isograph_vae/abundance_structure/axis_orthogonality.parquet`
has 13,222 rows, which is that analysis's `n_tested` — the BrainSEQ caudate
**schizophrenia** analysis. No other store has an `abundance_structure/` directory, so the
other 16 analyses contribute nothing to the claim. The manuscript text and the Fig 2b
legend were corrected on 2026-09-26 to name the single analysis, which is honest but
narrow: a reviewer will reasonably ask whether the channels are separable in the aging
analyses that the rest of the section is about.

`isograph_benchmark/real_data/abundance_structure_separation.py` already takes
`analysis` + `--region`, so no new estimator is needed.

Fix: add a stage script that runs `abundance_structure` across all 17 stores (the same
loop shape as `03_module_characterization/_h/01e`/`01f`), roll the per-analysis
`axis_orthogonality.parquet` up into one table with each analysis's median |r| and
fraction below 0.1, and facet `manuscript/_h/abundance_structure_figure.R` panel A by
analysis. Then either restore the general wording in the manuscript or keep the
single-analysis wording and cite the new supplement for the rest.

## 4. ~~The switch-unique definition has no threshold-sensitivity check~~ DONE 2026-09-26

`switch_unique_threshold.py` recounts all 17 analyses at FDR 0.05 / 0.10 / 0.20, before and after composition adjustment, and `figSwitchUniqueThreshold` (Fig S15) reports it. **The answer is mixed, and the manuscript was updated to say so.** The count ranking is essentially threshold-free (Spearman 0.97-0.98 against alpha = 0.10), and the two claims the text leans on hold at every alpha: GTEx cortex and frontal cortex BA9 collapse (persistence <= 0.004 throughout) while ACC BA24 and BrainSEQ DLPFC retain (>= 0.26 throughout). But the persistence *ordering* in the middle of the range is not stable (0.50 at alpha = 0.05 even after dropping analyses with fewer than 10 unadjusted genes), and GTEx hippocampus in particular runs 0.04 -> 0.48 -> 0.55 across the three alphas, so its "about half" is a property of alpha = 0.10. Kept for the record:

A gene is called switch-unique when its switch coordinate reaches FDR < 0.10 conditional
on abundance **and** the abundance test conditional on the switch coordinate does not.
That is an asymmetric use of one threshold: the second half is an acceptance of the null,
so a gene can move in or out of the class because its abundance test sits just either side
of 0.10. The manuscript already says the classification does not test whether the two
conditional effects differ, but nothing shows how much the regional pattern depends on
where the line is drawn — and the regional pattern is the result of the subsection.

Everything needed is already on disk: `incremental_association/gene_level.parquet` and
`incremental_association_composition/gene_level.parquet` carry `fdr_switch_given_abund`
and `fdr_abund_given_switch` per gene for all 17 analyses, so this is a recount, not a
refit.

Fix: add a CLI that recounts switch-unique genes at FDR ∈ {0.05, 0.10, 0.20} for every
analysis, before and after composition adjustment, together with `n_overlap`/`n_new` at
each threshold; report it as a supplementary figure (counts and gene-level persistence
versus threshold, one line per analysis) plus a supplementary table. The claim to support
is that the *ordering* of regions — and the cortical collapse in particular — is not an
artifact of the 0.10 cut. If it turns out to be threshold-sensitive, that belongs in the
text.

Done by: `real_data/switch_unique_threshold.py` + `_h/04c` (local, no scheduler) +
`manuscript/_h/switch_unique_threshold_figure.R`. The rollup reports rank stability both
over all analyses and restricted to those with at least 10 unadjusted switch-unique genes,
because persistence is a ratio and several analyses sit at 1-2 genes; and a per-analysis
call that is only made when it holds at every alpha.
