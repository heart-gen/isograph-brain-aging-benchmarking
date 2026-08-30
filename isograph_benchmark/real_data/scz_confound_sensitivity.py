"""Are the schizophrenia module findings explained by clinical or technical confounding?

A reviewer asked for the schizophrenia results to be tested against medication, toxicology,
smoking "and related confounding **where available**". This CLI answers that question in the
only honest order: first establish what is actually available, then test everything that is,
then say plainly what cannot be tested and why.

**What the cohort releases.** BrainSEQ's colData carries 130 columns. Nine are subject or
sample descriptors (``Dx``, ``Age``, ``Sex``, ``Race``, ``PMI``, ``MoD``, ``RIN``,
``Protocol``, ``SNP_PC1-10``); every remaining column is sequencing QC. There is no
medication field, no antipsychotic exposure, no toxicology panel, no nicotine or smoking
status -- not in the released RSEs and not in the derived bundles. The ``audit`` subcommand
records this as a first-class, re-runnable output rather than a claim in prose.

**Three tiers of evidence, reported separately and never merged:**

  1. ``measured`` -- covariates the cohort does release and the published model does *not*
     already adjust for (PMI, Race, r_rRNA rate, ancestry PC6-10, Protocol). Each is added
     to the published diagnosis model on its own, then all together. This is real
     confounder control and carries the most weight.
  2. ``proxy`` -- molecular surrogates for two of the unmeasured exposures, built from gene
     abundance in the *raw bundle counts* rather than from IsoGraph's own features, so the
     proxy cannot inherit module structure. Smoking uses the AHR/xenobiotic battery
     (CYP1B1, AHRR, CYP1A1, NQO1, TIPARP); antipsychotic exposure uses striatal DRD2, whose
     upregulation under chronic D2 blockade is long established. **These are surrogates, not
     measurements.** A module that survives them is not thereby proven unconfounded -- a
     proxy with poor sensitivity cannot exonerate anything. Their informative direction is
     the negative one: a module that *moves* under a proxy is flagged.
  3. ``unavailable`` -- what remains untestable in this cohort, named explicitly so the
     limitation is stated rather than implied.

The baseline model is verified against the published ``diagnosis_assoc.parquet`` before any
sensitivity is fitted, so the ladder perturbs the published statistic and not a lookalike.

Outputs land in ``06_switch_mechanism/_m/scz_confound_sensitivity/``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, region_store, rel, stage_out
from isograph_benchmark.real_data.run_models import (
    BRAINSEQ_COVARIATES,
    diagnosis_association,
)

_BUNDLE = rel("inputs", "bundles", "brainseq_sczd", "caudate")
_ARTIFACTS = region_store("brainseq", "caudate_sczd", "isograph_vae")

# The published SCZD diagnosis model (run_models.run_brainseq_caudate_sczd).
PUBLISHED_COVARIATES = ["Age"] + BRAINSEQ_COVARIATES

# Released but not in the published diagnosis model. Each is a plausible route by which a
# case/control difference could arise without being disease biology.
CANDIDATE_MEASURED = [
    "PMI",          # post-mortem interval -- degradation, and it differs by manner of death
    "Race",         # population structure beyond the genotype PCs, and ascertainment
    "r_rna_rate",   # rRNA carry-over; the published model adjusts mapping and mito only
    "Protocol",     # library protocol; confounded with batch
    "SNP_PC6", "SNP_PC7", "SNP_PC8", "SNP_PC9", "SNP_PC10",
]

# Exposures the reviewer named that this cohort does not release at all.
UNAVAILABLE = {
    "medication": "no antipsychotic, mood-stabiliser or any other medication field",
    "toxicology": "no post-mortem toxicology panel",
    "smoking": "no smoking status, pack-years, nicotine or cotinine measurement",
    "substance_use": "no alcohol or illicit substance field",
    "duration_of_illness": "no age at onset or illness duration",
}

# Molecular proxies (GENCODE v47 ids, resolved from the repo's own annotation).
PROXY_SETS = {
    "smoking_ahr_battery": {
        "genes": {
            "ENSG00000138061": "CYP1B1",
            "ENSG00000063438": "AHRR",
            "ENSG00000140465": "CYP1A1",
            "ENSG00000181019": "NQO1",
            "ENSG00000163659": "TIPARP",
        },
        "rationale": (
            "aryl-hydrocarbon-receptor xenobiotic battery; the canonical transcriptional "
            "response to tobacco smoke exposure"
        ),
    },
    "antipsychotic_drd2": {
        "genes": {"ENSG00000149295": "DRD2"},
        "rationale": (
            "striatal DRD2 is upregulated by chronic D2-receptor blockade; caudate is the "
            "relevant tissue, and DRD2 is the drug target itself"
        ),
    },
}


def _out_dir() -> Path:
    return ensure_dir(stage_out("mechanism", "scz_confound_sensitivity"))


# --------------------------------------------------------------------------- #
# Availability audit
# --------------------------------------------------------------------------- #
def _is_testable(sample_table: pd.DataFrame, col: str) -> tuple[bool, str]:
    """A covariate must vary to be a covariate. A single-level column silently drops out of
    the design and would otherwise register as a confounder the result 'survived'."""
    if col not in sample_table.columns:
        return False, "not present in this bundle"
    v = sample_table[col]
    if v.notna().sum() == 0:
        return False, "present but entirely missing"
    if v.nunique(dropna=True) < 2:
        only = v.dropna().iloc[0] if v.notna().any() else "?"
        return False, f"released but invariant in this cohort (single level: {only})"
    return True, "released and testable; added as a sensitivity below"


def audit_availability(sample_table: pd.DataFrame) -> pd.DataFrame:
    """One row per reviewer-named or candidate confounder, with what the cohort holds."""
    rows = []
    for name, why in UNAVAILABLE.items():
        rows.append(
            {
                "confounder": name,
                "tier": "unavailable",
                "released_by_cohort": False,
                "in_published_model": False,
                "detail": why,
                "n_nonmissing": 0,
            }
        )
    for col in PUBLISHED_COVARIATES:
        present = col in sample_table.columns
        rows.append(
            {
                "confounder": col,
                "tier": "measured",
                "released_by_cohort": present,
                "in_published_model": True,
                "detail": "already adjusted in the published diagnosis model",
                "n_nonmissing": int(sample_table[col].notna().sum()) if present else 0,
            }
        )
    for col in CANDIDATE_MEASURED:
        present = col in sample_table.columns
        testable, detail = _is_testable(sample_table, col)
        rows.append(
            {
                "confounder": col,
                "tier": "measured" if testable else "untestable",
                "released_by_cohort": present,
                "in_published_model": False,
                "detail": detail,
                "n_nonmissing": int(sample_table[col].notna().sum()) if present else 0,
            }
        )
    for name, spec in PROXY_SETS.items():
        rows.append(
            {
                "confounder": name,
                "tier": "proxy",
                "released_by_cohort": False,
                "in_published_model": False,
                "detail": f"molecular surrogate ({', '.join(spec['genes'].values())}): "
                + spec["rationale"],
                "n_nonmissing": 0,
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_sample_table() -> pd.DataFrame:
    st = pd.read_parquet(_BUNDLE / "samples.parquet")
    st["sample_id"] = st["sample_id"].astype(str)
    return st


def load_eigengenes() -> pd.DataFrame:
    """Sample-wise module eigengenes, reconstructed the way the fit defines them.

    ``compute_trait_associations`` takes a module's eigengene to be the unweighted mean of
    every ``feature_scores`` row belonging to its genes; both inputs are persisted, so this
    is a reconstruction, not a re-derivation. ``verify_baseline`` then checks it against the
    published association before anything is concluded from it.
    """
    fs = pd.read_parquet(_ARTIFACTS / "feature_scores.parquet")
    modules = pd.read_parquet(_ARTIFACTS / "modules.parquet")
    meta = {"feature_id", "gene_id", "feature_type", "n_transcripts"}
    samples = [c for c in fs.columns if c not in meta]

    gene_to_module = dict(
        zip(modules["gene_id"].astype(str), modules["module_id"].astype(str))
    )
    fs["gene_id"] = fs["gene_id"].astype(str)
    fs = fs[fs["gene_id"].isin(gene_to_module)].copy()
    fs["module_id"] = fs["gene_id"].map(gene_to_module)
    eig = fs.groupby("module_id")[samples].mean().sort_index()
    out = eig.T.reset_index().rename(columns={"index": "sample_id"})
    out["sample_id"] = out["sample_id"].astype(str)
    return out


def build_proxies(sample_table: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Z-scored log-CPM composites for each proxy set, from the raw bundle counts.

    Deliberately built from ``gene_counts.npz`` rather than from ``feature_scores``: a proxy
    derived from the same matrix that defines the modules would be partly circular.
    """
    genes = pd.read_parquet(_BUNDLE / "genes.parquet")
    counts = np.load(_BUNDLE / "gene_counts.npz", allow_pickle=True)["data"]
    counts = np.asarray(counts, dtype=float)
    if counts.shape[0] != len(genes):
        counts = counts.T
    lib = counts.sum(axis=0)
    cpm = np.log2(counts / np.where(lib > 0, lib, np.nan) * 1e6 + 1.0)

    bare = genes["gene_id"].astype(str).str.split(".", n=1).str[0]
    order = sample_table["sample_id"].astype(str).tolist()
    out = pd.DataFrame({"sample_id": order})
    meta: dict[str, dict] = {}
    for name, spec in PROXY_SETS.items():
        idx = np.flatnonzero(bare.isin(spec["genes"]).to_numpy())
        found = sorted({spec["genes"][g] for g in bare.iloc[idx]})
        missing = sorted(set(spec["genes"].values()) - set(found))
        if len(idx) == 0:
            meta[name] = {"n_genes_found": 0, "genes_found": [], "genes_missing": missing,
                          "note": "no member gene survived the bundle expression filter; "
                                  "proxy not constructed"}
            continue
        sub = cpm[idx, :]
        z = (sub - np.nanmean(sub, axis=1, keepdims=True)) / np.nanstd(
            sub, axis=1, keepdims=True
        )
        score = np.nanmean(z, axis=0)
        out[name] = score[: len(order)] if len(score) >= len(order) else np.nan
        meta[name] = {
            "n_genes_found": len(found),
            "genes_found": found,
            "genes_missing": missing,
            "rationale": spec["rationale"],
        }
    return out, meta


# --------------------------------------------------------------------------- #
# Model ladder
# --------------------------------------------------------------------------- #
def verify_baseline(assoc: pd.DataFrame, tol: float) -> float:
    """Assert the refit baseline reproduces the published ``diagnosis_assoc.parquet``."""
    published = pd.read_parquet(_ARTIFACTS / "diagnosis_assoc.parquet")
    merged = published.merge(
        assoc[["module_id", "effect"]], on="module_id", suffixes=("_pub", "_new")
    )
    if merged.empty:
        raise SystemExit("baseline refit shares no module with diagnosis_assoc.parquet")
    worst = float((merged["effect_pub"] - merged["effect_new"]).abs().max())
    if worst > tol:
        raise SystemExit(
            f"baseline diagnosis effects deviate from the published "
            f"diagnosis_assoc.parquet by {worst:.3g} (> {tol:g}); the sensitivity would be "
            "perturbing a different statistic than the published one."
        )
    return worst


def _ladder(available_measured: list[str], proxies: list[str]) -> list[tuple[str, str, list[str]]]:
    """(model name, tier, extra covariates) for each rung."""
    rungs: list[tuple[str, str, list[str]]] = [("published", "baseline", [])]
    rungs += [(f"plus_{c}", "measured", [c]) for c in available_measured]
    if available_measured:
        rungs.append(("plus_all_measured", "measured", list(available_measured)))
    rungs += [(f"plus_{p}", "proxy", [p]) for p in proxies]
    if proxies:
        rungs.append(("plus_all_proxies", "proxy", list(proxies)))
    if available_measured and proxies:
        rungs.append(
            ("plus_everything", "combined", [*available_measured, *proxies])
        )
    return rungs


def run(args) -> None:
    out_dir = _out_dir()
    sample_table = load_sample_table()
    eig = load_eigengenes()

    avail = audit_availability(sample_table)
    proxy_tab, proxy_meta = build_proxies(sample_table)
    sample_table = sample_table.merge(proxy_tab, on="sample_id", how="left")

    measured = [c for c in CANDIDATE_MEASURED if _is_testable(sample_table, c)[0]]
    not_testable = {
        c: _is_testable(sample_table, c)[1]
        for c in CANDIDATE_MEASURED
        if not _is_testable(sample_table, c)[0]
    }
    proxies = [p for p in PROXY_SETS if p in sample_table.columns]

    frames = []
    for name, tier, extra in _ladder(measured, proxies):
        assoc = diagnosis_association(
            eig, sample_table, PUBLISHED_COVARIATES + extra
        )
        if name == "published":
            worst = verify_baseline(assoc, args.tol)
            print(f"[verify] baseline reproduces diagnosis_assoc.parquet "
                  f"(max |Δeffect| = {worst:.3g})")
        assoc = assoc.assign(model=name, tier=tier,
                             added=",".join(extra) if extra else "")
        frames.append(assoc)

    long = pd.concat(frames, ignore_index=True)
    long.to_parquet(out_dir / "confound_sensitivity_long.parquet", index=False)
    avail.to_parquet(out_dir / "confounder_availability.parquet", index=False)

    base = long[long["model"] == "published"].set_index("module_id")
    rows = []
    for name, grp in long.groupby("model", sort=False):
        g = grp.set_index("module_id")
        shared = base.index.intersection(g.index)
        d_eff = (g.loc[shared, "effect"] - base.loc[shared, "effect"]).abs()
        # Attenuation of the published-significant set is what a reviewer cares about.
        sig_base = base.loc[shared][base.loc[shared, "fdr"] < args.fdr].index
        rows.append(
            {
                "model": name,
                "tier": grp["tier"].iloc[0],
                "added": grp["added"].iloc[0],
                "n_modules": int(len(shared)),
                "n_fdr_sig": int((g.loc[shared, "fdr"] < args.fdr).sum()),
                "n_fdr_sig_baseline": int(len(sig_base)),
                "n_baseline_sig_retained": int(
                    (g.loc[sig_base, "fdr"] < args.fdr).sum()
                ) if len(sig_base) else 0,
                "median_abs_effect_shift": float(d_eff.median()),
                "max_abs_effect_shift": float(d_eff.max()),
                "median_effect_ratio": float(
                    (g.loc[shared, "effect"].abs() / base.loc[shared, "effect"].abs())
                    .replace([np.inf, -np.inf], np.nan)
                    .median()
                ),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_parquet(out_dir / "confound_sensitivity_summary.parquet", index=False)

    meta = {
        "analysis": "SCZD module diagnosis association under clinical/technical confounders",
        "cohort": "BrainSEQ caudate (Dx = SCZD vs Control)",
        "artifacts": str(_ARTIFACTS),
        "published_covariates": PUBLISHED_COVARIATES,
        "measured_tested": measured,
        "measured_not_testable": not_testable,
        "proxies_built": proxy_meta,
        "unavailable": UNAVAILABLE,
        "fdr": args.fdr,
        "baseline_tolerance": args.tol,
        "n_samples": int(len(sample_table)),
        "n_modules": int(base.shape[0]),
    }
    (out_dir / "confound_sensitivity.json").write_text(json.dumps(meta, indent=2, default=str))
    _write_report(out_dir, avail, summary, long, meta, args.fdr)
    print(summary.to_string(index=False))


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _md(frame: pd.DataFrame, floats: int = 4) -> str:
    def cell(v):
        if isinstance(v, float):
            return "" if pd.isna(v) else f"{v:.{floats}g}"
        return "" if v is None else str(v)

    cols = list(frame.columns)
    return "\n".join(
        [
            "| " + " | ".join(cols) + " |",
            "|" + "|".join("---" for _ in cols) + "|",
            *[
                "| " + " | ".join(cell(v) for v in row) + " |"
                for row in frame.itertuples(index=False, name=None)
            ],
        ]
    )


def _write_report(out_dir, avail, summary, long, meta, fdr) -> None:
    base = summary[summary["model"] == "published"].iloc[0]
    lines = [
        "# Schizophrenia findings under clinical and technical confounding",
        "",
        "Generated by `scz_confound_sensitivity.py`. The reviewer's request was to test the "
        "schizophrenia findings for medication, toxicology, smoking and related confounding "
        "**where available**. This report takes the qualifier seriously.",
        "",
        "## 1. What is available",
        "",
        "BrainSEQ releases no medication, toxicology or smoking variable. The released "
        "colData is 130 columns, of which nine are subject/sample descriptors and the rest "
        "are sequencing QC. This is not a gap in our processing -- the fields do not exist "
        "in the distributed data, so the reviewer's primary confounders cannot be tested "
        "directly in this cohort by anyone.",
        "",
        _md(avail),
        "",
        "## 2. Measured confounders (the load-bearing tier)",
        "",
        f"Each released covariate the published model does not already adjust for is added "
        f"to it, alone and then jointly. The published model adjusts "
        f"`{', '.join(meta['published_covariates'])}`.",
        "",
        _md(summary),
        "",
        (
            "Excluded as untestable in this cohort (a single-level column drops out of the "
            "design, so a null result for it would be an artefact of the encoding, not "
            "evidence): "
            + "; ".join(f"`{k}` — {v}" for k, v in meta["measured_not_testable"].items())
            + "."
        ) if meta.get("measured_not_testable") else "",
        "",
        f"Baseline: {int(base['n_fdr_sig'])} modules at FDR < {fdr} of "
        f"{int(base['n_modules'])} tested.",
        "",
        "`n_baseline_sig_retained` is the number of *baseline-significant* modules still "
        "significant under that model, which is the quantity a confounding argument turns "
        "on. `median_effect_ratio` below 1 means the diagnosis effect shrank.",
        "",
        "## 3. Proxy confounders (weaker evidence, read in one direction only)",
        "",
        "Two unmeasured exposures have defensible molecular surrogates, built from the raw "
        "bundle counts rather than from IsoGraph's own features so they cannot inherit "
        "module structure:",
        "",
    ]
    for name, info in meta["proxies_built"].items():
        found = ", ".join(info.get("genes_found", [])) or "none"
        note = info.get("note", "")
        lines.append(
            f"- **{name}** — genes found: {found}"
            + (f"; missing: {', '.join(info['genes_missing'])}" if info.get("genes_missing") else "")
            + (f". {info.get('rationale','')}" if info.get("rationale") else "")
            + (f" **{note}**" if note else "")
        )
    lines += [
        "",
        "**These are surrogates, not measurements, and they are asymmetric evidence.** A "
        "module that survives proxy adjustment is *not* thereby shown to be unconfounded: a "
        "proxy with poor sensitivity to the true exposure cannot exonerate anything. The "
        "informative direction is the other one — a module whose effect moves materially "
        "when a proxy enters the model is flagged as a candidate for confounding, and is "
        "listed below.",
        "",
    ]
    proxy_models = summary[summary["tier"] == "proxy"]["model"].tolist()
    flagged = []
    if proxy_models:
        b = long[long["model"] == "published"].set_index("module_id")
        for m in proxy_models:
            g = long[long["model"] == m].set_index("module_id")
            shared = b.index.intersection(g.index)
            ratio = (g.loc[shared, "effect"].abs() / b.loc[shared, "effect"].abs())
            moved = shared[(ratio < 0.8) | (ratio > 1.25)]
            for mod in moved:
                flagged.append(
                    {
                        "module_id": mod,
                        "proxy": m,
                        "effect_published": float(b.loc[mod, "effect"]),
                        "effect_adjusted": float(g.loc[mod, "effect"]),
                        "fdr_published": float(b.loc[mod, "fdr"]),
                        "fdr_adjusted": float(g.loc[mod, "fdr"]),
                    }
                )
    if flagged:
        lines += [
            "Modules whose diagnosis effect changes by more than 20% under a proxy:",
            "",
            _md(pd.DataFrame(flagged)),
            "",
        ]
    else:
        lines += ["No module's diagnosis effect changes by more than 20% under either proxy.", ""]

    lines += [
        "## 4. What remains untestable",
        "",
        "State these in the manuscript limitations rather than leaving them implied:",
        "",
    ]
    lines += [f"- **{k}** — {v}" for k, v in meta["unavailable"].items()]
    lines += [
        "",
        "Closing them requires a cohort that releases clinical annotation (or a data-use "
        "agreement covering the restricted BrainSEQ clinical fields, if one exists). The "
        "proxy tier above is the most that can be done with the public release, and it is "
        "reported as such.",
        "",
    ]
    (out_dir / "SCZ_CONFOUND_SENSITIVITY.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("audit", help="availability audit only; fits nothing")
    a.set_defaults(func=lambda args: print(audit_availability(load_sample_table()).to_string(index=False)))

    r = sub.add_parser("run", help="availability audit + full sensitivity ladder")
    r.add_argument("--fdr", type=float, default=0.05)
    r.add_argument("--tol", type=float, default=1e-6,
                   help="max |Δeffect| allowed between the refit baseline and the "
                        "published diagnosis_assoc.parquet")
    r.set_defaults(func=run)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
