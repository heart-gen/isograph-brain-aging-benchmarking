"""Donor and tissue-sample sharing across the 16 regions and 17 analyses.

The manuscript is careful to say that the 16 aging analyses and the schizophrenia analysis
are 17 *applications*, not 17 independent cohorts. Nothing in the repository showed that,
so a reader had to take the caveat on trust, and a reviewer reading "replicates across
regions" had no way to see how much of that is re-measurement of the same donors.

Three identifier levels matter here, and conflating any two of them gives the wrong count:

* **Donor** -- ``BrNum`` (BrainSEQ) or ``SUBJID`` (GTEx). One person, shared across every
  region that person contributed tissue to.
* **Tissue sample** -- ``sample_id``, which is the ``RNum`` for BrainSEQ. Specific to a
  donor *and* a region, so it never repeats across regions.
* **Analysis** -- one model fit. Not the same as a region: BrainSEQ **caudate** is a single
  region used by two analyses. The aging analysis takes the 238 controls; the schizophrenia
  analysis takes those same 238 controls plus 152 patients, re-using the identical RNums.

So the design is **16 regions** (3 BrainSEQ + 13 GTEx) across **17 analyses**, and the
sample-analysis rows cover fewer distinct tissue samples than their sum, because the
caudate control samples are counted by both caudate analyses. Both levels are reported:
region overlap is the sampling design, analysis overlap is what a "replicates across
regions" claim is actually built on.

Outputs (all pure reads over the bundle sample tables -- nothing is refit):

1. ``donor_overlap.parquet``/``.csv`` -- every ordered analysis pair, shared donors and
   Jaccard; the diagonal carries the analysis's own donor count.
2. ``region_overlap.csv`` -- the same at region level, with each region's donor set taken
   as the union over the analyses that use it.
3. ``donor_incidence.csv`` -- long donor x region membership, and
   ``donor_analysis_counts.csv`` -- per donor, how many regions and analyses they appear in.
4. ``analysis_donor_counts.csv`` / ``region_donor_counts.csv`` and ``DONOR_STRUCTURE.md``.

Donor identifiers are only comparable within a cohort: a BrainSEQ ``BrNum`` and a GTEx
``SUBJID`` name different people and are never matched across cohorts, so cross-cohort
pairs are reported as zero overlap by construction rather than by test.

Run (local, no scheduler): ``bash 02_module_discovery/_h/03a.donor_structure.sh``
"""
from __future__ import annotations

import argparse
import itertools
import json

import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out

DONOR_COL = {"BrainSEQ": "BrNum", "GTEx": "SUBJID"}

# (analysis, region, cohort, bundle, phenotype). `bundle` is relative to inputs/bundles.
# Note the two caudate rows: one region, two analyses, nested sample sets.
ANALYSES: list[tuple[str, str, str, str, str]] = [
    ("BrainSEQ caudate (aging)", "BrainSEQ caudate", "BrainSEQ",
     "brainseq_v1/caudate", "aging"),
    ("BrainSEQ caudate (SCZD)", "BrainSEQ caudate", "BrainSEQ",
     "brainseq_sczd/caudate", "diagnosis"),
    ("BrainSEQ hippocampus", "BrainSEQ hippocampus", "BrainSEQ",
     "brainseq_v1/hippocampus", "aging"),
    ("BrainSEQ DLPFC", "BrainSEQ DLPFC", "BrainSEQ",
     "brainseq_v1/dlpfc", "aging"),
    ("GTEx amygdala", "GTEx amygdala", "GTEx",
     "gtex_v11_brain/amygdala", "aging"),
    ("GTEx ACC BA24", "GTEx ACC BA24", "GTEx",
     "gtex_v11_brain/anterior_cingulate_cortex_ba24", "aging"),
    ("GTEx caudate", "GTEx caudate", "GTEx",
     "gtex_v11_brain/caudate_basal_ganglia", "aging"),
    ("GTEx cerebellar hem.", "GTEx cerebellar hem.", "GTEx",
     "gtex_v11_brain/cerebellar_hemisphere", "aging"),
    ("GTEx cerebellum", "GTEx cerebellum", "GTEx",
     "gtex_v11_brain/cerebellum", "aging"),
    ("GTEx cortex", "GTEx cortex", "GTEx",
     "gtex_v11_brain/cortex", "aging"),
    ("GTEx frontal ctx BA9", "GTEx frontal ctx BA9", "GTEx",
     "gtex_v11_brain/frontal_cortex_ba9", "aging"),
    ("GTEx hippocampus", "GTEx hippocampus", "GTEx",
     "gtex_v11_brain/hippocampus", "aging"),
    ("GTEx hypothalamus", "GTEx hypothalamus", "GTEx",
     "gtex_v11_brain/hypothalamus", "aging"),
    ("GTEx n. accumbens", "GTEx n. accumbens", "GTEx",
     "gtex_v11_brain/nucleus_accumbens_basal_ganglia", "aging"),
    ("GTEx putamen", "GTEx putamen", "GTEx",
     "gtex_v11_brain/putamen_basal_ganglia", "aging"),
    ("GTEx spinal cord C1", "GTEx spinal cord C1", "GTEx",
     "gtex_v11_brain/spinal_cord_cervical_c_1", "aging"),
    ("GTEx substantia nigra", "GTEx substantia nigra", "GTEx",
     "gtex_v11_brain/substantia_nigra", "aging"),
]


def _load() -> tuple[dict[str, set[str]], dict[str, set[str]], pd.DataFrame]:
    """Donor ids and tissue-sample ids per analysis, plus the per-analysis count table.

    Bundles are already filtered to the analysed samples (``build_bundles`` drops
    ``dropped != "f"`` rows), so these are exactly the samples the modules were fit on.
    """
    donors, samples, rows = {}, {}, []
    for analysis, region, cohort, bundle, phenotype in ANALYSES:
        p = rel("inputs", "bundles", *bundle.split("/")) / "samples.parquet"
        if not p.exists():
            print(f"[donors] skipping {analysis}: no bundle at {p}", flush=True)
            continue
        s = pd.read_parquet(p)
        col = DONOR_COL[cohort]
        if col not in s.columns:
            raise SystemExit(f"{analysis}: sample table has no {col!r} column")
        donors[analysis] = set(s[col].dropna().astype(str))
        samples[analysis] = set(s["sample_id"].dropna().astype(str))
        rows.append({"analysis": analysis, "region": region, "cohort": cohort,
                     "phenotype": phenotype, "n_samples": int(len(s)),
                     "n_donors": len(donors[analysis])})
    return donors, samples, pd.DataFrame(rows)


def _check_identifier_levels(samples: dict[str, set[str]],
                             counts: pd.DataFrame) -> list[str]:
    """Assert the identifier semantics the rest of the module relies on.

    A tissue-sample id must never be shared by two *different* regions -- if it were, either
    the id is not region-specific or two bundles are the same tissue under two names, and
    every region-level count below would be wrong. Sharing within a region is expected, and
    is exactly the caudate case.
    """
    region_of = dict(zip(counts["analysis"], counts["region"]))
    notes = []
    for a, b in itertools.combinations(samples, 2):
        shared = samples[a] & samples[b]
        if not shared:
            continue
        if region_of[a] != region_of[b]:
            raise SystemExit(
                f"tissue-sample ids shared across regions ({a} / {b}): {len(shared)} ids. "
                "sample_id is expected to be region-specific (RNum for BrainSEQ)."
            )
        notes.append(f"{a} and {b} share {len(shared)} tissue samples "
                     f"within {region_of[a]}")
    # A donor appearing twice within one bundle would break the one-sample-per-donor
    # reading of the counts; say so rather than silently averaging over it.
    for analysis, row in counts.set_index("analysis").iterrows():
        if row["n_samples"] != row["n_donors"]:
            notes.append(f"{analysis}: {row['n_samples']} samples from "
                         f"{row['n_donors']} donors (repeat sampling within the analysis)")
    return notes


def _union_by(keys: pd.Series, values: pd.Series,
              sets: dict[str, set[str]]) -> dict[str, set[str]]:
    """Union the per-analysis sets up to whatever grouping `keys` names (here, region)."""
    out: dict[str, set[str]] = {}
    for key, analysis in zip(keys, values):
        out.setdefault(key, set()).update(sets.get(analysis, set()))
    return out


def _overlap(sets: dict[str, set[str]], cohort_of: dict[str, str],
             name: str) -> pd.DataFrame:
    """Shared-donor count and Jaccard for every ordered pair, diagonal included."""
    rows = []
    for a, b in itertools.product(sets, repeat=2):
        # Donor namespaces are cohort-local; a BrNum is never the same person as a SUBJID.
        same_cohort = cohort_of[a] == cohort_of[b]
        shared = len(sets[a] & sets[b]) if same_cohort else 0
        union = len(sets[a] | sets[b]) if same_cohort else len(sets[a]) + len(sets[b])
        rows.append({
            f"{name}_a": a, f"{name}_b": b,
            "cohort_a": cohort_of[a], "cohort_b": cohort_of[b],
            "n_donors_a": len(sets[a]), "n_donors_b": len(sets[b]),
            "n_shared_donors": shared,
            "frac_of_a_shared": shared / len(sets[a]) if sets[a] else float("nan"),
            "jaccard": shared / union if union else float("nan"),
        })
    return pd.DataFrame(rows)


def _incidence(region_donors: dict[str, set[str]],
               cohort_of: dict[str, str]) -> pd.DataFrame:
    """Long donor x region membership: one row per (donor, region) they contributed to."""
    return pd.DataFrame([{"cohort": cohort_of[region], "donor_id": d, "region": region}
                         for region, ds in region_donors.items() for d in sorted(ds)])


def _per_donor(incidence: pd.DataFrame, donors: dict[str, set[str]],
               counts: pd.DataFrame) -> pd.DataFrame:
    """Per donor: how many regions, and how many analyses, they appear in.

    The two differ exactly for the BrainSEQ caudate donors, who contribute one region but
    are carried by two analyses.
    """
    cohort_a = dict(zip(counts["analysis"], counts["cohort"]))
    n_analyses: dict[tuple[str, str], int] = {}
    for analysis, ds in donors.items():
        for d in ds:
            key = (cohort_a[analysis], d)
            n_analyses[key] = n_analyses.get(key, 0) + 1
    out = (incidence.groupby(["cohort", "donor_id"], as_index=False)
           .agg(n_regions=("region", "size")))
    out["n_analyses"] = [n_analyses[(c, d)] for c, d in
                         zip(out["cohort"], out["donor_id"])]
    return out.sort_values(["cohort", "n_regions", "n_analyses"],
                           ascending=[True, False, False], ignore_index=True)


def _md_table(df: pd.DataFrame) -> str:
    head = "| " + " | ".join(df.columns) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    body = ["| " + " | ".join("" if pd.isna(v) else
                              (f"{v:.3g}" if isinstance(v, float) else str(v))
                              for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def _top_pairs(overlap: pd.DataFrame, name: str, n: int = 10) -> pd.DataFrame:
    """The n most-overlapping unordered pairs; the matrix is symmetric by construction."""
    off = overlap[overlap[f"{name}_a"] != overlap[f"{name}_b"]]
    off = off.assign(_pair=[frozenset((a, b)) for a, b in
                            zip(off[f"{name}_a"], off[f"{name}_b"])])
    return (off.sort_values("n_shared_donors", ascending=False)
            .drop_duplicates(subset="_pair")
            .head(n)[[f"{name}_a", f"{name}_b", "n_shared_donors", "jaccard"]])


def run() -> None:
    donors, samples, counts = _load()
    if counts.empty:
        raise SystemExit("no bundles found under inputs/bundles")
    notes = _check_identifier_levels(samples, counts)

    cohort_r = dict(zip(counts["region"], counts["cohort"]))
    cohort_a = dict(zip(counts["analysis"], counts["cohort"]))
    region_donors = _union_by(counts["region"], counts["analysis"], donors)
    region_samples = _union_by(counts["region"], counts["analysis"], samples)

    region_counts = (pd.DataFrame([
        {"region": r, "cohort": cohort_r[r], "n_donors": len(region_donors[r]),
         "n_tissue_samples": len(region_samples[r]),
         "n_analyses": int((counts["region"] == r).sum())}
        for r in region_donors])
        .sort_values(["cohort", "region"], ignore_index=True))

    analysis_overlap = _overlap(donors, cohort_a, "analysis")
    region_overlap = _overlap(region_donors, cohort_r, "region")
    incidence = _incidence(region_donors, cohort_r)
    per_donor = _per_donor(incidence, donors, counts)

    out = ensure_dir(stage_out("modules", "_m", "donor_structure"))
    analysis_overlap.to_parquet(out / "donor_overlap.parquet", index=False,
                                compression="zstd")
    # CSVs alongside: the figure reads these, so it builds with an `arrow` compiled
    # without zstd support.
    analysis_overlap.to_csv(out / "donor_overlap.csv", index=False)
    region_overlap.to_csv(out / "region_overlap.csv", index=False)
    counts.to_csv(out / "analysis_donor_counts.csv", index=False)
    region_counts.to_csv(out / "region_donor_counts.csv", index=False)
    per_donor.to_csv(out / "donor_analysis_counts.csv", index=False)
    incidence.to_csv(out / "donor_incidence.csv", index=False)

    n_unique_samples = len(set().union(*samples.values()))
    by_cohort = (per_donor.groupby("cohort")
                 .agg(n_donors=("donor_id", "size"),
                      median_regions=("n_regions", "median"),
                      max_regions=("n_regions", "max"),
                      n_in_one_region_only=("n_regions", lambda s: int((s == 1).sum())))
                 .reset_index())
    summary = {
        "n_regions": int(len(region_counts)),
        "n_analyses": int(len(counts)),
        "n_sample_analysis_rows": int(counts["n_samples"].sum()),
        "n_unique_tissue_samples": n_unique_samples,
        "n_unique_donors": int(len(per_donor)),
        "max_pairwise_shared_donors_regions": int(
            region_overlap.loc[region_overlap["region_a"] != region_overlap["region_b"],
                               "n_shared_donors"].max()),
        "by_cohort": by_cohort.to_dict(orient="records"),
        "identifier_notes": notes,
    }
    (out / "donor_structure_summary.json").write_text(json.dumps(summary, indent=2))

    lines = [
        "# Donor and tissue-sample sharing across the 16 regions and 17 analyses", "",
        f"{len(region_counts)} regions are analysed by {len(counts)} analyses. The "
        f"{int(counts['n_samples'].sum()):,} sample-analysis rows cover "
        f"{n_unique_samples:,} distinct tissue samples from {len(per_donor):,} donors: "
        "BrainSEQ caudate is one region used by two analyses, so its control samples are "
        "counted by both. Donor ids are cohort-local (`BrNum` for BrainSEQ, `SUBJID` for "
        "GTEx) and are never matched across cohorts; `sample_id` is the tissue-specific "
        "id (`RNum` for BrainSEQ) and is checked above to be region-specific.", "",
        "## Per region", "", _md_table(region_counts), "",
        "## Per analysis", "", _md_table(counts), "",
        "## Donors per cohort", "", _md_table(by_cohort), "",
        "## Most-overlapping region pairs", "",
        _md_table(_top_pairs(region_overlap, "region")), "",
        "## Identifier notes", "",
        *([f"- {n}" for n in notes] or ["- none"]), "",
        "## Interpretation", "",
        "- **These are 17 applications over 16 regions, not 17 independent cohorts.** The "
        "GTEx regions are largely re-measurements of one donor pool, the three BrainSEQ "
        "regions are largely the same brains, and the two BrainSEQ caudate analyses are "
        "the same region: the aging analysis is the 238 controls, and the schizophrenia "
        "analysis is those same 238 control samples plus 152 patients.",
        "- Any statement that a result 'replicates across regions' within a cohort is "
        "therefore a statement about the same donors measured in different tissue, not "
        "about independent samples. The cross-cohort comparisons (BrainSEQ vs GTEx) are "
        "the ones that carry independent donors, and they are also the ones that cross a "
        "transcript-processing pipeline boundary.",
        "- The split-half and permutation nulls elsewhere in the repository are computed "
        "within an analysis, so they are unaffected by this sharing; it bounds how much "
        "*between-analysis* agreement should be read as independent replication.",
    ]
    (out / "DONOR_STRUCTURE.md").write_text("\n".join(lines))
    print(f"[donors] {len(region_counts)} regions, {len(counts)} analyses, "
          f"{n_unique_samples:,} tissue samples, {len(per_donor):,} donors", flush=True)
    for n in notes:
        print(f"[donors]   note: {n}", flush=True)
    print(f"[donors] wrote {out/'DONOR_STRUCTURE.md'}", flush=True)


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    run()


if __name__ == "__main__":
    main()
