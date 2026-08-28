"""Fragmentation-sensitive partition metrics for the planted-module benchmark.

Best-match Jaccard (``module_recovery_score``) scores each ground-truth module against
its single best predicted module, so a method that shatters the network into many small
modules is never penalised for the shattering itself — a signal-free WGCNA fit reaches a
recovery of ~0.4 by emitting ~282 tiny modules.  These functions add the two corrections:

* whole-partition comparisons (adjusted Rand, adjusted mutual information, and the
  homogeneity / completeness / V-measure triple that *diagnoses* fragmentation vs merging);
* a module-count-preserving permutation null for the existing best-match Jaccard, so the
  score can be reported as an excess / z-score over what the fitted partition's own module
  size distribution yields by chance.

Pure functions with no I/O, so `run_one.compute_metrics` (live runs) and
`backfill_metrics` (the 13k already-completed runs) share one implementation.

The null reproduces ``module_recovery_score``'s conventions exactly, which is subtler than
it looks: the union in |T∩P|/|T∪P| is taken against the **full** predicted module — all of
its genes, including the background genes that were never planted.  Restricting the
computation to planted genes would shrink every denominator and inflate the score, so
predicted module sizes are always carried at full width and only the intersections are
restricted to planted genes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    homogeneity_completeness_v_measure,
)

# Label for genes a method left out of every module (IsoGraph drops components below
# min_module_size; WGCNA drops grey).  Truth tables carry planted genes only, so in the
# whole-universe variant the non-planted genes get their own truth class.
UNASSIGNED = "__unassigned__"
BACKGROUND = "__background__"

DEFAULT_N_PERM = 200
DEFAULT_SEED = 13

METRIC_KEYS = (
    "ari_planted", "ami_planted", "ari_planted_blob",
    "homogeneity_planted", "completeness_planted", "v_measure_planted",
    "ari_assigned", "ami_assigned",
    "ari_universe", "ami_universe",
    "frac_planted_assigned", "n_truth_modules",
    "module_recovery_recomputed",
    "module_recovery_null_mean", "module_recovery_null_sd",
    "module_recovery_excess", "module_recovery_z", "module_recovery_perm_p",
)


def _pairs(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    return frame["gene_id"].astype(str).to_numpy(), frame["module_id"].astype(str).to_numpy()


# --------------------------------------------------------------------------- #
# Whole-partition comparisons (item 1)
# --------------------------------------------------------------------------- #
def partition_labels(
    predicted: pd.DataFrame,
    truth: pd.DataFrame,
    universe: list[str] | None = None,
    singleton_unassigned: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Aligned (truth, predicted) integer label vectors over a fixed gene universe.

    ``universe`` defaults to the planted genes.  Genes the method did not assign become
    their own singleton class by default — the neutral convention, which neither rewards a
    method for dropping genes nor charges it with one giant false merge.  Setting
    ``singleton_unassigned=False`` pools them into one class instead (WGCNA-grey
    semantics); both are reported rather than one being defended.
    """
    t_genes, t_mods = _pairs(truth)
    truth_map = dict(zip(t_genes, t_mods))
    if predicted is None or predicted.empty:
        pred_map: dict[str, str] = {}
    else:
        p_genes, p_mods = _pairs(predicted)
        pred_map = dict(zip(p_genes, p_mods))

    genes = sorted(truth_map) if universe is None else sorted({str(g) for g in universe})
    t_labels, p_labels = [], []
    for i, gene in enumerate(genes):
        t_labels.append(truth_map.get(gene, BACKGROUND))
        if gene in pred_map:
            p_labels.append(pred_map[gene])
        else:
            p_labels.append(f"{UNASSIGNED}{i}" if singleton_unassigned else UNASSIGNED)
    return (
        pd.factorize(pd.Index(t_labels), sort=True)[0],
        pd.factorize(pd.Index(p_labels), sort=True)[0],
    )


# --------------------------------------------------------------------------- #
# Best-match Jaccard and its module-count-preserving null (item 4)
# --------------------------------------------------------------------------- #
def _mrs_from_counts(counts: np.ndarray, truth_sizes: np.ndarray, pred_sizes: np.ndarray) -> np.ndarray:
    """Mean over truth modules of max_j |T∩P_j| / (|T_i| + |P_j| − |T∩P_j|).

    ``counts`` may be 2-D (n_truth × n_pred) or 3-D (n_draw × n_truth × n_pred).
    """
    union = truth_sizes[..., :, None] + pred_sizes[..., None, :] - counts
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = np.where(union > 0, counts / union, 0.0)
    return jaccard.max(axis=-1).mean(axis=-1)


def _mrs_arrays(predicted: pd.DataFrame, truth: pd.DataFrame):
    """Codes and sizes for the best-match Jaccard, in module_recovery_score's convention.

    Returns ``(pred_codes, truth_codes_of_shared, shared_index, truth_sizes, pred_sizes)``
    where ``shared_index`` positions into the predicted gene vector for genes that are also
    planted.  Predicted sizes are full-width (all assigned genes, planted or not); truth
    sizes are full-width too (all planted genes, assigned or not).
    """
    p_genes, p_mods = _pairs(predicted)
    t_genes, t_mods = _pairs(truth)

    pred_codes, _ = pd.factorize(pd.Index(p_mods), sort=True)
    truth_codes_all, _ = pd.factorize(pd.Index(t_mods), sort=True)
    pred_sizes = np.bincount(pred_codes)
    truth_sizes = np.bincount(truth_codes_all)

    # Positions in the predicted vector for genes that are also planted, and the truth
    # module each of those genes belongs to.  Only these contribute intersections.
    truth_of_gene = dict(zip(t_genes, truth_codes_all))
    shared_index, shared_truth = [], []
    for i, gene in enumerate(p_genes):
        code = truth_of_gene.get(gene)
        if code is not None:
            shared_index.append(i)
            shared_truth.append(code)
    return (
        pred_codes,
        np.asarray(shared_truth, dtype=np.int64),
        np.asarray(shared_index, dtype=np.int64),
        truth_sizes.astype(np.int64),
        pred_sizes.astype(np.int64),
    )


def best_match_jaccard(predicted: pd.DataFrame, truth: pd.DataFrame) -> float:
    """Reimplementation of ``isograph.evaluation.metrics.module_recovery_score``.

    Kept here so the permutation null and the published metric are provably the same
    computation; ``backfill_metrics`` asserts the two agree on every run.
    """
    if predicted is None or predicted.empty or truth is None or truth.empty:
        return 0.0
    pred_codes, shared_truth, shared_index, truth_sizes, pred_sizes = _mrs_arrays(predicted, truth)
    n_t, n_p = len(truth_sizes), len(pred_sizes)
    counts = np.zeros((n_t, n_p), dtype=np.int64)
    if len(shared_index):
        flat = shared_truth * n_p + pred_codes[shared_index]
        counts = np.bincount(flat, minlength=n_t * n_p).reshape(n_t, n_p)
    return float(_mrs_from_counts(counts, truth_sizes, pred_sizes))


def size_preserving_null(
    predicted: pd.DataFrame,
    truth: pd.DataFrame,
    n_perm: int = DEFAULT_N_PERM,
    seed: int = DEFAULT_SEED,
    chunk_cells: int = 4_000_000,
) -> np.ndarray:
    """Best-match Jaccard under relabelling that preserves the predicted module sizes.

    Which genes fill each predicted module is randomised; the module-size multiset and the
    set of assigned genes are both held exactly fixed.  This isolates how much recovery a
    partition of this granularity earns for free — the calibration the ~282-module
    signal-free WGCNA fit needs.
    """
    if predicted is None or predicted.empty or truth is None or truth.empty:
        return np.empty(0)
    pred_codes, shared_truth, shared_index, truth_sizes, pred_sizes = _mrs_arrays(predicted, truth)
    n_t, n_p = len(truth_sizes), len(pred_sizes)
    if not len(shared_index) or n_perm <= 0:
        return np.empty(0)

    rng = np.random.default_rng(seed)
    n_assigned = len(pred_codes)
    chunk = max(1, min(n_perm, chunk_cells // max(n_assigned, 1)))
    out = np.empty(n_perm)
    done = 0
    while done < n_perm:
        b = min(chunk, n_perm - done)
        # Size-preserving relabel: one independent permutation of the label vector per draw.
        perm = np.argsort(rng.random((b, n_assigned)), axis=1)
        labels = pred_codes[perm][:, shared_index]              # b × n_shared
        flat = shared_truth[None, :] * n_p + labels
        counts = np.stack([np.bincount(f, minlength=n_t * n_p).reshape(n_t, n_p) for f in flat])
        out[done:done + b] = _mrs_from_counts(counts, truth_sizes[None, :], pred_sizes[None, :])
        done += b
    return out


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def partition_metrics(
    predicted: pd.DataFrame,
    truth: pd.DataFrame,
    universe: list[str] | None = None,
    seed: int = DEFAULT_SEED,
    n_perm: int = DEFAULT_N_PERM,
) -> dict[str, float | None]:
    """All item-1 and item-4 metrics for one fitted network.

    ``predicted``/``truth`` are ``gene_id``/``module_id`` frames; ``universe`` is the
    simulation's full gene list (``truth_switch.parquet`` covers all genes, whereas
    ``truth_modules.parquet`` holds only the planted ones).  Values are ``None`` when not
    estimable, so degenerate runs drop out of the summaries rather than averaging in as
    zeros.
    """
    out: dict[str, float | None] = dict.fromkeys(METRIC_KEYS)
    if truth is None or truth.empty or predicted is None or predicted.empty:
        return out

    out["n_truth_modules"] = float(truth["module_id"].nunique())

    # A partition comparison needs at least two items to compare.  ``negative_control_noise``
    # plants a *single* gene in a *single* module, so every method scores a vacuous
    # ari/homogeneity/completeness of 1.0 there — precisely the kind of free score these
    # metrics exist to remove.  Return None so the scenario drops out of the partition
    # summaries; its diagnosis comes from the size-preserving null below, which is
    # well-defined regardless, plus the raw predicted module count.
    if truth["gene_id"].nunique() < 2:
        obs = best_match_jaccard(predicted, truth)
        out["module_recovery_recomputed"] = obs
        _add_null(out, obs, predicted, truth, n_perm, seed)
        return out

    t_pl, p_pl = partition_labels(predicted, truth, singleton_unassigned=True)
    out["ari_planted"] = float(adjusted_rand_score(t_pl, p_pl))
    out["ami_planted"] = float(adjusted_mutual_info_score(t_pl, p_pl))
    hom, com, vms = homogeneity_completeness_v_measure(t_pl, p_pl)
    out["homogeneity_planted"] = float(hom)
    out["completeness_planted"] = float(com)
    out["v_measure_planted"] = float(vms)

    _, p_blob = partition_labels(predicted, truth, singleton_unassigned=False)
    out["ari_planted_blob"] = float(adjusted_rand_score(t_pl, p_blob))

    planted = set(truth["gene_id"].astype(str))
    assigned = set(predicted["gene_id"].astype(str))
    covered = planted & assigned
    out["frac_planted_assigned"] = float(len(covered) / len(planted)) if planted else None

    # Coverage-conditional variant: the same comparison restricted to the planted genes the
    # method actually placed in a module.  ``ari_planted`` conflates two very different
    # failures — leaving planted genes unassigned (a min_module_size / coverage property)
    # and grouping the genes it did assign incorrectly (a real accuracy failure).  The gap
    # between ``ari_planted`` and ``ari_assigned`` is exactly the first of those, so the two
    # together say which one a low score is made of.
    if len(covered) >= 2:
        t_cov, p_cov = partition_labels(predicted, truth, universe=sorted(covered))
        out["ari_assigned"] = float(adjusted_rand_score(t_cov, p_cov))
        out["ami_assigned"] = float(adjusted_mutual_info_score(t_cov, p_cov))

    if universe is not None and len(universe):
        t_all, p_all = partition_labels(predicted, truth, universe=universe, singleton_unassigned=True)
        out["ari_universe"] = float(adjusted_rand_score(t_all, p_all))
        out["ami_universe"] = float(adjusted_mutual_info_score(t_all, p_all))

    obs = best_match_jaccard(predicted, truth)
    out["module_recovery_recomputed"] = obs
    _add_null(out, obs, predicted, truth, n_perm, seed)
    return out


def _add_null(out: dict, obs: float, predicted: pd.DataFrame, truth: pd.DataFrame,
              n_perm: int, seed: int) -> None:
    """Calibrate the observed best-match Jaccard against the size-preserving null."""
    null = size_preserving_null(predicted, truth, n_perm=n_perm, seed=seed)
    null = null[np.isfinite(null)]
    if not null.size:
        return
    null_mean, null_sd = float(null.mean()), float(null.std(ddof=1))
    out["module_recovery_null_mean"] = null_mean
    out["module_recovery_null_sd"] = null_sd
    out["module_recovery_excess"] = float(obs - null_mean)
    out["module_recovery_z"] = float((obs - null_mean) / null_sd) if null_sd > 0 else None
    out["module_recovery_perm_p"] = float((1 + int((null >= obs).sum())) / (null.size + 1))
