# Open issues found during manuscript review (2026-09-22)

## 1. Module-recovery axes are labelled "AUC" but the metric is best-match Jaccard

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

## 2. Gene-level overlap of switch-unique sets before/after composition adjustment is not in the pipeline

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
