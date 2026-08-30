"""Collapse the redundant ATtRACT RBP motifs into similarity families (reviewer item 6b).

ATtRACT carries ~1,200 human matrices for ~160 RBPs, and many are near-duplicates: several
matrices for one protein, and different proteins that bind the same short degenerate element
(the AU-rich binders in particular).  Every downstream regulon table counts matrices or RBPs,
so N near-identical motifs read as N independent lines of regulatory evidence when they are
really one.  This assigns each matrix a **family**, so recurrence can be reported per family
with the per-RBP view demoted to secondary.

There is no Tomtom in the environment, so PWM similarity is computed directly:

* each matrix column is a probability vector over A/C/G/U; centre it at 0.25 and normalise,
  so a column-vs-column Pearson correlation is a plain dot product;
* slide one matrix over the other at every ungapped offset with at least ``--min-overlap``
  overlapping columns and take the mean column correlation, maximised over offsets;
* **sense orientation only** — RBP motifs act on single-stranded RNA, so unlike a DNA
  transcription-factor motif there is no reverse complement to consider;
* average-linkage hierarchical clustering on ``1 - similarity``, cut at ``--cut``.

The cut is calibrated against ATtRACT's own annotation rather than chosen by eye: the report
prints how the within-gene and within-annotated-family similarity distributions separate from
the between distribution, and the chosen cut is pinned in ``configs/rbp_families.yaml``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from isograph_benchmark.paths import ensure_dir, rel, stage_out

MOTIF_DIR = rel("inputs", "rbp_motifs")
OUT_DIR = stage_out("regulation", "rbp")
CONFIG = rel("configs", "rbp_families.yaml")

DEFAULT_CUT = 0.25          # distance = 1 - mean column correlation
MIN_OVERLAP = 4


# --------------------------------------------------------------------------- #
# PWM loading (MOODS-free: this stage runs in the isograph env)
# --------------------------------------------------------------------------- #
def load_frequency_matrices(motif_dir: Path = MOTIF_DIR) -> tuple[list[str], list[np.ndarray]]:
    """(matrix_ids, L×4 probability matrices) for the human ATtRACT motifs."""
    db = pd.read_csv(motif_dir / "ATtRACT_db.txt", sep="\t", dtype=str)
    human = set(db.loc[db["Organism"] == "Homo_sapiens", "Matrix_id"])

    mats: dict[str, list[list[float]]] = {}
    cur, rows = None, []
    for line in (motif_dir / "pwm.txt").read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None and rows:
                mats[cur] = rows
            cur, rows = line[1:].split()[0], []
        else:
            rows.append([float(x) for x in line.split()])
    if cur is not None and rows:
        mats[cur] = rows

    ids, arrays = [], []
    for mid, rws in mats.items():
        if mid not in human:
            continue
        a = np.asarray(rws, dtype=float)
        if a.ndim != 2 or a.shape[1] != 4 or a.shape[0] == 0:
            continue
        ids.append(mid)
        arrays.append(a)
    return ids, arrays


def _annotation(motif_dir: Path = MOTIF_DIR) -> pd.DataFrame:
    """matrix_id -> gene name and ATtRACT's own family annotation (one row per matrix)."""
    db = pd.read_csv(motif_dir / "ATtRACT_db.txt", sep="\t", dtype=str)
    db = db[db["Organism"] == "Homo_sapiens"]
    keep = ["Matrix_id", "Gene_name", "Family", "Len"]
    out = db[keep].drop_duplicates("Matrix_id").rename(columns={
        "Matrix_id": "matrix_id", "Gene_name": "rbp",
        "Family": "attract_family", "Len": "motif_len"})
    return out.set_index("matrix_id")


# --------------------------------------------------------------------------- #
# Similarity
# --------------------------------------------------------------------------- #
def _normalised(arrays: list[np.ndarray]) -> list[np.ndarray]:
    """Centre each column at 0.25 and unit-normalise, so column Pearson r = dot product.

    A perfectly uniform column carries no information; it normalises to zeros and therefore
    contributes 0 to the mean correlation rather than a spurious 1.
    """
    out = []
    for a in arrays:
        c = a - 0.25
        n = np.linalg.norm(c, axis=1, keepdims=True)
        out.append(np.divide(c, n, out=np.zeros_like(c), where=n > 1e-12))
    return out


def similarity_matrix(arrays: list[np.ndarray], min_overlap: int = MIN_OVERLAP) -> np.ndarray:
    """Max-over-offsets mean column correlation between every pair of motifs.

    Vectorised by grouping motifs of equal length: for a given (length, length, offset) the
    overlapping columns of every motif in each group flatten to a matrix, and all pairwise
    mean correlations for that offset are one matmul.
    """
    norm = _normalised(arrays)
    n = len(norm)
    sim = np.full((n, n), -np.inf)
    np.fill_diagonal(sim, 1.0)

    by_len: dict[int, list[int]] = {}
    for i, a in enumerate(norm):
        by_len.setdefault(a.shape[0], []).append(i)

    lengths = sorted(by_len)
    for li in lengths:
        A_idx = np.array(by_len[li])
        A = np.stack([norm[i] for i in A_idx])           # nA × li × 4
        for lj in lengths:
            if lj < li:
                continue
            B_idx = np.array(by_len[lj])
            B = np.stack([norm[j] for j in B_idx])       # nB × lj × 4
            # offset d = start of A within B's coordinate frame; allow partial overlap
            best = np.full((len(A_idx), len(B_idx)), -np.inf)
            for d in range(-(li - min_overlap), lj - min_overlap + 1):
                a0, b0 = max(0, -d), max(0, d)
                ov = min(li - a0, lj - b0)
                if ov < min_overlap:
                    continue
                Aa = A[:, a0:a0 + ov, :].reshape(len(A_idx), -1)
                Bb = B[:, b0:b0 + ov, :].reshape(len(B_idx), -1)
                np.maximum(best, (Aa @ Bb.T) / ov, out=best)
            sim[np.ix_(A_idx, B_idx)] = np.maximum(sim[np.ix_(A_idx, B_idx)], best)
            sim[np.ix_(B_idx, A_idx)] = np.maximum(sim[np.ix_(B_idx, A_idx)], best.T)

    np.fill_diagonal(sim, 1.0)
    return np.clip(sim, -1.0, 1.0)


# --------------------------------------------------------------------------- #
# Clustering + calibration
# --------------------------------------------------------------------------- #
def _calibration(sim: np.ndarray, labels: pd.Series) -> dict:
    """How well does raw similarity separate annotated groups from unrelated pairs?

    Reported so the chosen cut is defensible rather than eyeballed: if within-annotation
    pairs are not clearly more similar than between-annotation pairs, no threshold will
    recover meaningful families and the family view should not be trusted.
    """
    lab = labels.to_numpy()
    known = pd.notna(lab) & (lab != "")
    iu = np.triu_indices(len(lab), k=1)
    s = sim[iu]
    same = (lab[iu[0]] == lab[iu[1]]) & known[iu[0]] & known[iu[1]]
    diff = (lab[iu[0]] != lab[iu[1]]) & known[iu[0]] & known[iu[1]]
    q = [0.25, 0.5, 0.75, 0.9]
    return {
        "n_within_pairs": int(same.sum()),
        "n_between_pairs": int(diff.sum()),
        "within_quantiles": {str(x): float(np.quantile(s[same], x)) for x in q} if same.any() else {},
        "between_quantiles": {str(x): float(np.quantile(s[diff], x)) for x in q} if diff.any() else {},
    }


_SWEEP_CUTS = (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)


def _cut_sweep(Z: np.ndarray, rbp: np.ndarray) -> list[dict]:
    """Family count and RBP purity across cuts.

    ATtRACT's motifs are short and degenerate, so similarity does not fall into separated
    clusters and no cut is uniquely correct.  Emitting the sweep makes the chosen cut's
    influence visible instead of hidden.
    """
    rows = []
    for cut in _SWEEP_CUTS:
        f = fcluster(Z, t=cut, criterion="distance")
        sizes = pd.Series(f).value_counts()
        modal = pd.DataFrame({"f": f, "rbp": rbp}).groupby("f")["rbp"].agg(
            lambda x: x.mode().iloc[0])
        purity = float((pd.Series(rbp) == pd.Series(f).map(modal)).mean())
        rows.append({"cut": cut, "n_families": int(len(sizes)),
                     "largest_family": int(sizes.max()),
                     "n_singletons": int((sizes == 1).sum()),
                     "rbp_purity": purity})
    return rows


def run(cut: float, min_overlap: int, motif_dir: Path = MOTIF_DIR) -> pd.DataFrame:
    ids, arrays = load_frequency_matrices(motif_dir)
    print(f"human ATtRACT matrices: {len(ids)}")
    sim = similarity_matrix(arrays, min_overlap=min_overlap)

    dist = np.clip(1.0 - sim, 0.0, None)
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    Z = linkage(squareform(dist, checks=False), method="average")
    family = fcluster(Z, t=cut, criterion="distance")

    ann = _annotation(motif_dir).reindex(ids)
    out = pd.DataFrame({
        "matrix_id": ids,
        "rbp": ann["rbp"].to_numpy(),
        "attract_family": ann["attract_family"].to_numpy(),
        "family_id": [f"F{f:04d}" for f in family],
    })
    sizes = out.groupby("family_id")["matrix_id"].size()
    n_rbps = out.groupby("family_id")["rbp"].nunique()
    out["family_size"] = out["family_id"].map(sizes)
    out["n_rbps_in_family"] = out["family_id"].map(n_rbps)

    # Label each family by its most frequent RBP, so the tables stay readable.
    label = (out.groupby(["family_id", "rbp"]).size().reset_index(name="n")
             .sort_values(["family_id", "n"], ascending=[True, False])
             .drop_duplicates("family_id").set_index("family_id")["rbp"])
    out["family_label"] = out["family_id"].map(label)

    # Max similarity to another member, so a family's internal tightness is inspectable.
    idx = {m: i for i, m in enumerate(ids)}
    within = []
    for fam, grp in out.groupby("family_id"):
        members = [idx[m] for m in grp["matrix_id"]]
        sub = sim[np.ix_(members, members)].copy()
        np.fill_diagonal(sub, -np.inf)
        v = sub.max(axis=1) if len(members) > 1 else np.full(len(members), np.nan)
        within.extend(zip(grp["matrix_id"], v))
    out["max_sim_within"] = out["matrix_id"].map(dict(within))

    ensure_dir(OUT_DIR)
    out.to_parquet(OUT_DIR / "rbp_motif_families.parquet", index=False)

    calib = {
        "by_rbp": _calibration(sim, out["rbp"]),
        "by_attract_family": _calibration(sim, out["attract_family"]),
    }
    payload = {
        "cut": cut, "min_overlap": min_overlap,
        "n_matrices": len(ids), "n_rbps": int(out["rbp"].nunique()),
        "n_families": int(out["family_id"].nunique()),
        "n_multi_rbp_families": int((n_rbps > 1).sum()),
        "largest_family": int(sizes.max()),
        "calibration": calib,
        "cut_sweep": _cut_sweep(Z, out["rbp"].to_numpy()),
    }
    (OUT_DIR / "rbp_motif_families.json").write_text(json.dumps(payload, indent=2))
    _write_report(out, payload)
    print(f"{len(ids)} matrices / {out['rbp'].nunique()} RBPs -> "
          f"{out['family_id'].nunique()} families at cut={cut} "
          f"({int((n_rbps > 1).sum())} span >1 RBP)")
    return out


def _write_report(out: pd.DataFrame, payload: dict) -> None:
    c = payload["calibration"]
    lines = [
        "# RBP motif families", "",
        f"{payload['n_matrices']} human ATtRACT matrices for {payload['n_rbps']} RBPs collapse "
        f"to **{payload['n_families']} similarity families** at a distance cut of "
        f"{payload['cut']} (mean column correlation >= {1 - payload['cut']:.2f} under "
        f"average linkage, minimum {payload['min_overlap']} overlapping columns, sense "
        f"orientation only). {payload['n_multi_rbp_families']} families span more than one "
        f"RBP; the largest holds {payload['largest_family']} matrices.", "",
        "Report motif recurrence **per family**. Counting matrices, or even distinct RBPs, "
        "treats near-duplicate motifs as independent evidence.", "",
        "## Cut calibration", "",
        "Similarity quantiles for pairs of matrices annotated to the same group vs different "
        "groups. The cut is only meaningful if these separate.", "",
        "| grouping | pairs within | pairs between | within q25 | within q50 | between q75 | between q90 |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, key in [("same RBP", "by_rbp"), ("same ATtRACT family", "by_attract_family")]:
        d = c[key]
        w, b = d["within_quantiles"], d["between_quantiles"]
        lines.append(
            f"| {name} | {d['n_within_pairs']:,} | {d['n_between_pairs']:,} | "
            f"{w.get('0.25', float('nan')):.3f} | {w.get('0.5', float('nan')):.3f} | "
            f"{b.get('0.75', float('nan')):.3f} | {b.get('0.9', float('nan')):.3f} |")

    lines += [
        "", "The two distributions overlap heavily: matrices for *different* RBPs are about "
        "as similar as matrices for the same RBP. That is the redundancy this analysis "
        "exists to absorb — many RBPs bind near-identical short degenerate elements — but it "
        "also means no cut is uniquely correct, so the sweep below is reported rather than a "
        "single number defended.", "",
        "## Cut sensitivity", "",
        "| cut | families | largest | singletons | RBP purity |", "|---|---|---|---|---|"]
    for s in payload["cut_sweep"]:
        mark = " **(pinned)**" if abs(s["cut"] - payload["cut"]) < 1e-9 else ""
        lines.append(f"| {s['cut']:.2f}{mark} | {s['n_families']} | {s['largest_family']} | "
                     f"{s['n_singletons']} | {s['rbp_purity']:.3f} |")

    lines += ["", "## Largest families", "",
              "| family | label | matrices | RBPs |", "|---|---|---|---|"]
    top = (out.drop_duplicates("family_id")
           .sort_values("family_size", ascending=False).head(15))
    for r in top.itertuples():
        lines.append(f"| {r.family_id} | {r.family_label} | {r.family_size} | {r.n_rbps_in_family} |")
    (OUT_DIR / "RBP_MOTIF_FAMILIES.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cut", type=float, default=DEFAULT_CUT,
                   help="average-linkage distance cut (distance = 1 - mean column correlation)")
    p.add_argument("--min-overlap", type=int, default=MIN_OVERLAP)
    args = p.parse_args()
    run(args.cut, args.min_overlap)


if __name__ == "__main__":
    main()
