# Benchmark Metrics

## Module Recovery Score

**Definition.** Let $\mathcal{T} = \{T_1, \dots, T_K\}$ be the set of ground-truth
modules and $\mathcal{P} = \{P_1, \dots, P_M\}$ be the set of predicted modules.
The module recovery score is the mean best-Jaccard similarity:

$$\text{MRS}(\mathcal{P}, \mathcal{T}) = \frac{1}{K} \sum_{i=1}^{K} \max_{j \in [M]} \frac{|T_i \cap P_j|}{|T_i \cup P_j|}$$

Each ground-truth module $T_i$ is matched to the predicted module $P_j$ that
maximises the Jaccard index.  The score ranges from 0 (no overlap between any
truth and predicted module) to 1 (every truth module is exactly recovered by some
predicted module).

**Properties:**

- **Asymmetric in direction**: measures how well predicted modules cover the truth,
  not how well truth covers predictions.  A method that predicts one giant module
  covering all genes will score poorly because per-module Jaccard is diluted by the
  large union.
- **Penalises fragmentation**: if a truth module $T_i$ is split across two predicted
  modules $P_a$ and $P_b$, the best Jaccard is less than 1 even if $T_i \subseteq P_a \cup P_b$.
- **Independent of module count**: number of predicted modules affects the score
  only through partition quality, not directly.

**Implementation.** `isograph.evaluation.metrics.module_recovery_score(predicted, truth)`.
Both arguments are DataFrames with columns `gene_id` and `module_id`.

**Reference.** The best-match Jaccard approach for cluster recovery evaluation is
described in Lancichinetti & Fortunato (2009) *Physical Review E* 80:056117 and is
used by the igraph `compare_communities` benchmark suite under the name
"best-match F1 / Jaccard".

---

## Switch Gene Detection Rate

Fraction of ground-truth switching genes assigned to any predicted module:

$$\text{SGDR} = \frac{|\hat{G} \cap G_\text{switch}|}{|G_\text{switch}|}$$

where $\hat{G}$ is the set of all genes assigned to any predicted module and
$G_\text{switch}$ is the set of genes whose PSI varies with the simulated trait.

---

## Non-Switch Gene Module Rate (False Positive Rate)

Fraction of non-switching background genes that are assigned to a predicted module:

$$\text{NGMR} = \frac{|\hat{G} \cap G_\text{background}|}{|G_\text{background}|}$$

Lower is better.  A method that assigns every gene to a module scores 1.0 regardless
of specificity.

---

## Statistical Testing

Pairwise comparisons between each IsoGraph variant and WGCNA use a **paired
Wilcoxon signed-rank test** (two-sided).  Pairing is by synthetic dataset identity
(`dataset_id` / seed) so each method pair is evaluated on the same data.
Multiple-comparison correction uses **Benjamini-Hochberg FDR** applied jointly
across all (scenario × metric × method) tests.  Significance thresholds: ** FDR < 0.05,
\* FDR < 0.10.

Implementation: `isograph_benchmark/stats/hypothesis_tests.py::paired_tests()`.
