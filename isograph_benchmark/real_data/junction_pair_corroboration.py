"""Reproducible Figure 6c: donor support for nominated, exact transcript pairs.

Table-only analysis: no IsoGraph model fitting, BAM access, or HPC dependency.
Run from the repository root via 06_switch_mechanism/_h/02g.junction_pair_corroboration.sh.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy import stats

REGIONS = ("dlpfc", "hippocampus", "caudate")
# Generic cortex is deliberately distinct from an exact BA9/DLPFC match.
MAPPINGS = {
    "Brain_Frontal_Cortex_BA9": [("dlpfc", "exact"), ("caudate", "secondary")],
    "Brain_Cortex": [("dlpfc", "cortex_proxy"), ("caudate", "secondary")],
    "Brain_Hippocampus": [("hippocampus", "exact")],
    "Brain_Caudate_basal_ganglia": [("caudate", "exact")],
    "Brain_Anterior_cingulate_cortex_BA24": [("dlpfc", "adjacent")],
    "Brain_Putamen_basal_ganglia": [("caudate", "adjacent")],
    "Brain_Nucleus_accumbens_basal_ganglia": [("caudate", "adjacent")],
}
MATCH_ORDER = {"exact": 0, "cortex_proxy": 1, "adjacent": 2, "secondary": 3}
PRIMARY = {"exact", "cortex_proxy"}
COVARIATES = ["Age", "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
              "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"]
EVENT_PATH = "05_genetic_anchoring/_m/coloc_signal_susie/all_introns/coloc_isoform_events.parquet"
OUT_PATH = "06_switch_mechanism/_m/junction_pair_corroboration"


def canonical_pair(t1, t2):
    """Unordered identity, retaining transcript versions (no silent version rescue)."""
    return "|".join(sorted((str(t1), str(t2))))


def nominated_pairs(event, pairs):
    """Use actual tissue-specific edges, not combinations of flattened members."""
    gene = str(event["gene"]).split(".")[0]
    carriers = {t.strip() for t in str(event["junction_transcripts"]).split(",") if t.strip()}
    keep = pairs.gene_id.str.split(".").str[0].eq(gene)
    keep &= pairs.transcript_id_1.isin(carriers) | pairs.transcript_id_2.isin(carriers)
    return pairs.loc[keep].drop_duplicates(["transcript_id_1", "transcript_id_2"])


def junction_match(junction, pair, junctions):
    """LeafCutter flanking exon bases -> 1-based closed intron, exact strand/ends."""
    m = re.fullmatch(r"(chr[^:]+):(\d+)-(\d+)\(([+-])\)", str(junction))
    if not m:
        return False, False
    chrom, start, end, strand = m.groups()
    if strand != str(pair["strand"]):
        return False, False
    hit = junctions[(junctions.chrom == chrom)
                    & (junctions.start == int(start) + 1)
                    & (junctions.end == int(end) - 1)]
    return not hit.empty, bool(hit.specific.any())


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


class Inputs:
    def __init__(self, root):
        self.root, self.paths = root, set()

    def path(self, relative):
        p = self.root / relative
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p}")
        with p.open("rb") as f:
            if f.read(42).startswith(b"version https://git-lfs.github.com"):
                raise ValueError(f"Git LFS pointer, not data: {p}")
        self.paths.add(p)
        return p

    def frame(self, relative, **kwargs):
        return pd.read_parquet(self.path(relative), **kwargs)


def build_targets(inputs):
    events = inputs.frame(EVENT_PATH)
    events = events[events.junction_in_switch_pair.fillna(False)].copy()
    if events.empty:
        raise ValueError("No signal-level anchored nominations")
    events = events.reset_index(drop=True)
    events["event_id"] = [f"signal_{i:04d}" for i in range(len(events))]
    ledgers = {}
    rows = []
    for event in events.to_dict("records"):
        tissue = event["iso_region"]
        if tissue not in ledgers:
            ledgers[tissue] = inputs.frame(
                f"02_module_discovery/gtex/{tissue}/_m/isograph_vae/module_interpret/structure_switch_pairs.parquet")
        pairs = nominated_pairs(event, ledgers[tissue])
        mappings = MAPPINGS.get(event["tissue"], [("", "none")])
        # Keep even a nomination whose current pair ledger cannot resolve its members.
        records = pairs.to_dict("records") or [{"transcript_id_1": "", "transcript_id_2": ""}]
        for pair in records:
            for region, match in mappings:
                rows.append({
                    "event_id": event["event_id"], "gene": event["gene"],
                    "gene_name": event["gene_name"], "trait": event["trait"],
                    "qtl_tissue": event["tissue"], "iso_region": tissue,
                    "PP4_sQTL": event["PP4_sQTL"], "junction": event["junction"],
                    "t1": pair["transcript_id_1"], "t2": pair["transcript_id_2"],
                    "pair_key": (canonical_pair(pair["transcript_id_1"], pair["transcript_id_2"])
                                 if not pairs.empty else f"{event['gene']}|unresolved"),
                    "region": region, "match": match,
                    "in_nomination_pair_ledger": not pairs.empty,
                })
    return events, pd.DataFrame(rows).drop_duplicates()


def resolve_recount(targets, inputs):
    result = targets.copy()
    result["structural_status"] = "no_matched_region"
    result["count_pair_id"] = ""
    result["anchor_specific"] = False
    result["two_sided"] = False
    for region in REGIONS:
        base = f"06_switch_mechanism/_m/ase_junction_switch/{region}"
        pairs = inputs.frame(f"{base}/pairs.parquet")
        pairs["pair_key"] = [canonical_pair(a, b) for a, b in
                             zip(pairs.transcript_id_1, pairs.transcript_id_2)]
        if pairs.pair_key.duplicated().any():
            raise ValueError(f"Ambiguous recount orientation: {region}")
        by_pair = pairs.set_index("pair_key").to_dict("index")
        wanted = set(result.loc[result.region.eq(region), "pair_key"])
        wanted_ids = set(pairs.loc[pairs.pair_key.isin(wanted), "pair_id"])
        jn = inputs.frame(f"{base}/junctions.parquet", filters=[("pair_id", "in", sorted(wanted_ids))])
        jby = {k: v for k, v in jn.groupby("pair_id", sort=False)}
        for idx, row in result[result.region == region].iterrows():
            status = "pair_not_counted"
            pair = by_pair.get(row.pair_key)
            if not row.in_nomination_pair_ledger:
                status = "nomination_pair_unresolved"
            elif pair is not None:
                result.at[idx, "count_pair_id"] = pair["pair_id"]
                # Retain the recount's orientation, even if source ledger was reversed.
                result.at[idx, "t1"] = pair["transcript_id_1"]
                result.at[idx, "t2"] = pair["transcript_id_2"]
                if str(pair["gene_id"]).split(".")[0] != str(row.gene).split(".")[0]:
                    raise ValueError("Gene mismatch between nomination and recount")
                hit, specific = junction_match(row.junction, pair,
                                              jby.get(pair["pair_id"], jn.iloc[:0]))
                two = bool(pd.notna(pair["two_sided"]) and pair["two_sided"])
                result.at[idx, "anchor_specific"] = specific
                result.at[idx, "two_sided"] = two
                status = ("anchor_not_in_pair_model" if not hit else
                          "not_two_sided" if not two else "eligible")
            result.at[idx, "structural_status"] = status
    return result


def aggregate_counts(path, eligible_samples, pair_ids, batch_size=250_000):
    """Stream all regions independently; pool hap 0/1/2; exclude other samples."""
    keys = ["sample_id", "donor_id", "pair_id", "isoform"]
    parts = []
    pset, sset = pa.array(sorted(pair_ids)), pa.array(sorted(eligible_samples))
    if not pair_ids:
        return pd.DataFrame(columns=keys + ["n_frag"])
    for batch in pq.ParquetFile(path).iter_batches(
            batch_size=batch_size, columns=keys + ["hap", "n_frag"]):
        mask = pc.and_(pc.is_in(batch.column("pair_id"), value_set=pset),
                       pc.is_in(batch.column("sample_id"), value_set=sset))
        chosen = batch.filter(mask).to_pandas()
        if chosen.empty:
            continue
        if not chosen.hap.isin([0, 1, 2]).all():
            raise ValueError("Unexpected haplotype labels")
        if not chosen.isoform.isin([1, 2]).all() or (chosen.n_frag < 0).any():
            raise ValueError("Invalid fragment counts")
        parts.append(chosen.groupby(keys, as_index=False, observed=True).n_frag.sum())
    if not parts:
        return pd.DataFrame(columns=keys + ["n_frag"])
    return pd.concat(parts).groupby(keys, as_index=False, observed=True).n_frag.sum()


def sample_membership(inputs, region):
    sample = inputs.frame(f"inputs/bundles/brainseq_v1/{region}/samples.parquet")
    sample["sample_id"] = sample.sample_id.astype(str)
    sample["donor_id"] = sample.BrNum.astype(str)
    if sample.sample_id.duplicated().any() or sample.donor_id.duplicated().any():
        raise ValueError(f"Repeated samples/donors need explicit aggregation: {region}")
    if not sample.Dx.eq("Control").all() or (sample.Age < 18).any():
        raise ValueError("Expected adult control aging bundle")
    qc = inputs.frame(f"06_switch_mechanism/_m/ase_junction_switch/{region}/count_qc.parquet")
    if qc.sample_id.duplicated().any():
        raise ValueError("Duplicate count QC sample")
    completed = sample.merge(qc[["sample_id", "donor_id"]], on="sample_id", how="inner",
                             suffixes=("", "_count"), validate="one_to_one")
    if not completed.donor_id.eq(completed.donor_id_count).all():
        raise ValueError("Bundle/count donor identifiers disagree")
    return sample, completed.drop(columns="donor_id_count")


def donor_grid(counts, completed, pairs):
    """Explicit zeros only for completed samples, never for missing BAM/count runs."""
    base = pairs[["count_pair_id", "pair_key", "t1", "t2", "gene", "gene_name"]].drop_duplicates()
    grid = base.merge(completed[["sample_id", "donor_id"]], how="cross")
    if grid.empty:
        return grid.assign(n1=pd.Series(dtype=float), n2=pd.Series(dtype=float))
    if not counts.empty:
        expected = completed.set_index("sample_id").donor_id
        if not counts.donor_id.eq(counts.sample_id.map(expected)).all():
            raise ValueError("Count sample/donor mismatch")
        c = counts.pivot(index=["sample_id", "donor_id", "pair_id"],
                         columns="isoform", values="n_frag").fillna(0)
        c = c.reindex(columns=[1, 2], fill_value=0).rename(columns={1: "n1", 2: "n2"}).reset_index()
        grid = grid.merge(c, left_on=["sample_id", "donor_id", "count_pair_id"],
                          right_on=["sample_id", "donor_id", "pair_id"], how="left",
                          validate="one_to_one").drop(columns="pair_id")
    else:
        grid["n1"], grid["n2"] = 0, 0
    grid[["n1", "n2"]] = grid[["n1", "n2"]].fillna(0).astype("int64")
    grid["total"] = grid.n1 + grid.n2
    grid["junction_fraction_1"] = grid.n1 / grid.total.replace(0, np.nan)
    return grid


def rank_concordance(frame, x, y, min_n=30):
    """Descriptive partial Spearman; rank both outcomes, adjust common covariates.

    No nominal P value: both measurements use the same reads and nominated pairs.
    """
    missing = [c for c in COVARIATES if c not in frame]
    if missing:
        return {"concordance_status": "missing_covariates:" + ",".join(missing),
                "n_concordance": 0, "rho_raw": np.nan, "rho_adjusted": np.nan}
    f = frame[[x, y] + COVARIATES].replace([np.inf, -np.inf], np.nan).dropna()
    out = {"n_concordance": len(f), "rho_raw": np.nan, "rho_adjusted": np.nan,
           "concordance_status": "too_few_complete_donors"}
    if len(f) < min_n:
        return out
    if f[x].nunique() < 3 or f[y].nunique() < 3:
        out["concordance_status"] = "insufficient_variation"
        return out
    out["rho_raw"] = float(stats.spearmanr(f[x], f[y]).statistic)
    numeric = [c for c in COVARIATES if c not in ("Sex", "MoD", "Age")]
    blocks = [np.ones((len(f), 1))]
    # A smooth age adjustment (centered/scaled quadratic and cubic), chosen before run.
    for col in ["Age"] + numeric:
        v = pd.to_numeric(f[col], errors="raise").to_numpy(float)
        sd = np.std(v)
        if sd > 0:
            z = (v - np.mean(v)) / sd
            blocks.append(z[:, None])
            if col == "Age":
                blocks.extend([z[:, None] ** 2, z[:, None] ** 3])
    for col in ("Sex", "MoD"):
        blocks.append(pd.get_dummies(f[col].astype(str), drop_first=True, dtype=float).to_numpy())
    design = np.column_stack(blocks)
    rank = np.linalg.matrix_rank(design)
    if len(f) - rank < 10:
        out["concordance_status"] = "insufficient_residual_degrees_of_freedom"
        return out
    ranks = np.column_stack([stats.rankdata(f[x]), stats.rankdata(f[y])])
    residuals = ranks - design @ np.linalg.lstsq(design, ranks, rcond=None)[0]
    if np.min(np.std(residuals, axis=0)) < 1e-8:
        out["concordance_status"] = "constant_residual"
        return out
    out["rho_adjusted"] = float(np.corrcoef(residuals.T)[0, 1])
    out["concordance_status"] = "estimated_descriptive"
    return out


def add_quantification(grid, inputs, region, sample):
    """Use the original transcript estimates; length-normalize before pair ratio."""
    if grid.empty:
        return grid.assign(transcript_fraction_1=pd.Series(dtype=float))
    path = inputs.path(f"inputs/processed/brainseq/{region}/tx_counts.parquet")
    names = set(pq.ParquetFile(path).schema.names)
    ids = sorted(set(grid.sample_id) & names)
    tx = sorted(set(grid.t1) | set(grid.t2))
    q = pd.read_parquet(path, columns=["Name", "EffectiveLength"] + ids,
                       filters=[("Name", "in", tx)]).set_index("Name")
    if q.index.duplicated().any():
        raise ValueError("Duplicate transcript quantification rows")
    # Shared effective-length field is the release's available normalization.
    abundance = q[ids].div(q.EffectiveLength.where(q.EffectiveLength > 0), axis=0)
    values = []
    for row in grid.itertuples():
        a = abundance.at[row.t1, row.sample_id] if row.t1 in abundance.index and row.sample_id in ids else np.nan
        b = abundance.at[row.t2, row.sample_id] if row.t2 in abundance.index and row.sample_id in ids else np.nan
        values.append(a / (a + b) if np.isfinite(a + b) and a + b > 0 else np.nan)
    out = grid.copy()
    out["transcript_fraction_1"] = values
    return out.merge(sample[["sample_id"] + COVARIATES], on="sample_id", validate="many_to_one")


def pair_statistics(grid, depths, min_donors, min_per_form):
    rows = []
    for key, g in grid.groupby("pair_key", sort=True):
        for depth in depths:
            covered = g[g.total >= depth]
            both = (covered.n1 >= min_per_form) & (covered.n2 >= min_per_form)
            n = len(covered)
            row = {"pair_key": key, "min_total_fragments": depth,
                   "min_fragments_per_form": min_per_form,
                   "n_completed_donors": len(g), "n_covered": n,
                   "n_both": int(both.sum()),
                   "fraction_both": float(both.mean()) if n else np.nan,
                   "n_positive_total": int((g.total > 0).sum()),
                   "median_total_fragments": float(covered.total.median()) if n else np.nan,
                   "median_minor_usage": float(np.minimum(covered.junction_fraction_1,
                                                             1-covered.junction_fraction_1).median()) if n else np.nan,
                   "measurement_status": "measured" if n >= min_donors else "insufficient_coverage"}
            row.update(rank_concordance(covered, "junction_fraction_1", "transcript_fraction_1", min_donors))
            rows.append(row)
    return pd.DataFrame(rows)


def collapse_targets(targets):
    """One pair-region outcome, preserving all event provenance in targets.csv."""
    rows = []
    for (region, key), group in targets[targets.region != ""].groupby(["region", "pair_key"]):
        # Use strongest anatomical match, then prefer an eligible anchor *within* it.
        g = group.assign(match_rank=group.match.map(MATCH_ORDER),
                         ineligible=group.structural_status.ne("eligible"))
        first = g.sort_values(["match_rank", "ineligible", "event_id"]).iloc[0].to_dict()
        first["traits"] = ",".join(sorted(group.trait.unique()))
        first["qtl_tissues"] = ",".join(sorted(group.qtl_tissue.unique()))
        first["n_nominating_events"] = group.event_id.nunique()
        first["is_primary"] = first["match"] in PRIMARY
        rows.append(first)
    return pd.DataFrame(rows).drop(columns=["match_rank", "ineligible"], errors="ignore")


def save_frame(frame, out, name):
    frame.to_csv(out / f"{name}.csv", index=False)
    frame.to_parquet(out / f"{name}.parquet", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--output-dir", type=Path)
    ap.add_argument("--min-total-fragments", type=int, default=10)
    ap.add_argument("--min-donors", type=int, default=30)
    ap.add_argument("--depth-sensitivity", type=int, nargs="+", default=[5, 20])
    ap.add_argument("--batch-size", type=int, default=250_000)
    args = ap.parse_args()
    if min(args.min_total_fragments, args.min_donors, args.batch_size, *args.depth_sensitivity) < 1:
        ap.error("Thresholds and batch size must be positive")
    start = time.monotonic()
    root = args.root.resolve()
    out = args.output_dir or root / OUT_PATH
    out.mkdir(parents=True, exist_ok=True)
    depths = sorted(set([args.min_total_fragments] + args.depth_sensitivity))
    params = {"primary_min_total_fragments": args.min_total_fragments,
              "min_covered_donors": args.min_donors, "depth_sensitivity": depths,
              "primary_min_fragments_per_form": 1, "per_form_sensitivity": 2,
              "primary_matches": sorted(PRIMARY), "match_mappings": MAPPINGS,
              "covariates": COVARIATES, "age_adjustment": "centered_scaled_cubic",
              "transcript_identity": "exact_versioned_unordered_pair",
              "junction_coordinates": "LeafCutter_exon_flanks_to_closed_intron_plus1_minus1",
              "batch_size": args.batch_size, "no_allelic_significance_filter": True}
    (out / "params.json").write_text(json.dumps(params, indent=2) + "\n")
    inputs = Inputs(root)
    events, targets = build_targets(inputs)
    targets = resolve_recount(targets, inputs)
    save_frame(events, out, "nominations")
    save_frame(targets, out, "target_event_pairs")
    collapsed = collapse_targets(targets)
    all_stats, all_donors, memberships = [], [], []
    for region in REGIONS:
        sample, completed = sample_membership(inputs, region)
        memberships.append({"region": region, "n_bundle_donors": len(sample),
                            "n_count_completed": len(completed)})
        selected = collapsed[(collapsed.region == region) & (collapsed.structural_status == "eligible")]
        print(f"{region}: {len(selected)} exact nominated pairs; {len(completed)}/{len(sample)} count-completed donors", flush=True)
        if selected.empty:
            continue
        path = inputs.path(f"06_switch_mechanism/_m/ase_junction_switch/{region}/junction_allelic_counts.parquet")
        counts = aggregate_counts(path, set(completed.sample_id), set(selected.count_pair_id), args.batch_size)
        donors = donor_grid(counts, completed, selected)
        donors = add_quantification(donors, inputs, region, sample)
        donors["region"] = region
        all_donors.append(donors)
        for min_per_form in (1, 2):
            s = pair_statistics(donors, depths, args.min_donors, min_per_form)
            s["region"] = region
            all_stats.append(s)
    if not all_stats:
        raise ValueError("No countable nominated pairs; inspect target_event_pairs.csv")
    donor_data = pd.concat(all_donors, ignore_index=True)
    save_frame(donor_data, out, "donor_pair_counts")
    s = pd.concat(all_stats, ignore_index=True)
    # Cross every target with every threshold, retaining structural failures explicitly.
    conditions = pd.DataFrame([(d, f) for d in depths for f in (1, 2)],
                              columns=["min_total_fragments", "min_fragments_per_form"])
    outcomes = collapsed.merge(conditions, how="cross").merge(
        s, on=["region", "pair_key", "min_total_fragments", "min_fragments_per_form"],
        how="left", validate="one_to_one")
    outcomes["measurement_status"] = outcomes.measurement_status.fillna(outcomes.structural_status)
    outcomes["n_covered"] = outcomes.n_covered.fillna(0).astype(int)
    outcomes["n_both"] = outcomes.n_both.fillna(0).astype(int)
    save_frame(outcomes, out, "pair_results_all_thresholds")
    primary = outcomes[(outcomes.min_total_fragments == args.min_total_fragments)
                       & (outcomes.min_fragments_per_form == 1) & outcomes.is_primary].copy()
    save_frame(primary, out, "primary_pair_results")
    measured = primary[primary.measurement_status == "measured"]
    gene_rows = []
    for (gene, region), group in primary.groupby(["gene_name", "region"], sort=True):
        m = group[group.measurement_status == "measured"]
        gene_rows.append({"gene_name": gene, "region": region, "n_target_pairs": len(group),
                          "n_measured_pairs": len(m),
                          "median_fraction_both": m.fraction_both.median(),
                          "min_fraction_both": m.fraction_both.min(),
                          "max_fraction_both": m.fraction_both.max(),
                          "min_covered_donors": m.n_covered.min(),
                          "max_covered_donors": m.n_covered.max(),
                          "median_adjusted_rho": m.rho_adjusted.median(),
                          "n_concordance_pairs": m.rho_adjusted.notna().sum()})
    save_frame(pd.DataFrame(gene_rows), out, "primary_gene_region_summary")
    status = primary.groupby("measurement_status").agg(
        pair_regions=("pair_key", "size"), genes=("gene", "nunique")).reset_index()
    save_frame(status, out, "primary_coverage")
    sens = []
    for label, mask in [("primary", outcomes.is_primary), ("strict_anatomy", outcomes.match.eq("exact")),
                        ("adjacent_only", outcomes.match.eq("adjacent")),
                        ("secondary_only", outcomes.match.eq("secondary")),
                        ("primary_anchor_specific", outcomes.is_primary & outcomes.anchor_specific)]:
        for (depth, per_form), g in outcomes[mask].groupby(["min_total_fragments", "min_fragments_per_form"]):
            m = g[g.measurement_status == "measured"]
            # Equal-weight gene summary; pairs and regions are not independent replicates.
            gene_means = m.groupby("gene").fraction_both.mean()
            sens.append({"scope": label, "min_total_fragments": depth, "min_fragments_per_form": per_form,
                         "n_genes_targeted": g.gene.nunique(), "n_pair_regions_targeted": len(g),
                         "n_genes_measured": m.gene.nunique(), "n_pair_regions_measured": len(m),
                         "gene_balanced_mean_fraction_both": gene_means.mean()})
    save_frame(pd.DataFrame(sens), out, "sensitivity_summary")
    genes_mapped = set(primary.gene)
    genes_measured = set(measured.gene)
    gene_audit = events[["gene", "gene_name"]].drop_duplicates().copy()
    gene_audit["primary_region_available"] = gene_audit.gene.isin(genes_mapped)
    gene_audit["has_measured_primary_pair"] = gene_audit.gene.isin(genes_measured)
    save_frame(gene_audit, out, "gene_coverage")
    summary = {"n_nominated_genes": events.gene.nunique(), "n_nominating_events": len(events),
               "n_primary_region_genes": len(genes_mapped), "n_primary_pair_regions": len(primary),
               "n_structurally_eligible_primary_pair_regions": int(primary.structural_status.eq("eligible").sum()),
               "n_measured_genes": len(genes_measured), "n_measured_pair_regions": len(measured),
               "n_measured_distinct_pairs": measured.pair_key.nunique(),
               "both_form_donor_threshold": args.min_donors,
               "n_pairs_both_in_min_donors": int((measured.n_both >= args.min_donors).sum()),
               "fraction_both_min": measured.fraction_both.min(), "fraction_both_max": measured.fraction_both.max(),
               "covered_donors_min": int(measured.n_covered.min()) if len(measured) else 0,
               "covered_donors_max": int(measured.n_covered.max()) if len(measured) else 0,
               "n_descriptive_concordances": int(measured.rho_adjusted.notna().sum()),
               "median_descriptive_adjusted_rho": measured.rho_adjusted.median(),
               "sample_membership": memberships,
               "primary_status_counts": primary.measurement_status.value_counts().to_dict(),
               "nominated_genes_without_primary_region": sorted(gene_audit.loc[~gene_audit.primary_region_available, "gene_name"]),
               "primary_region_genes_without_measured_pairs": sorted(gene_audit.loc[gene_audit.primary_region_available & ~gene_audit.has_measured_primary_pair, "gene_name"])}
    def clean(obj):
        if isinstance(obj, dict): return {k: clean(v) for k, v in obj.items()}
        if isinstance(obj, list): return [clean(v) for v in obj]
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (float, np.floating)): return float(obj) if np.isfinite(obj) else None
        return obj
    summary = clean(summary)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    caption = (
        "**Pair-discriminating short-read junction support for disease-prioritized transcript pairs.** "
        f"Each small point is one exact nominated transcript pair; large diamonds show gene–region medians "
        f"and horizontal lines span the observed pair range, not a confidence interval. The x axis is the "
        f"fraction of adult control donors with at least one fragment assigned to each alternative, among donors "
        f"with at least {args.min_total_fragments} total pair-discriminating fragments. Pairs require at least "
        f"{args.min_donors} covered donors. Right-hand labels give measured/targeted pairs and the range of "
        f"covered donor counts. Of {summary['n_nominated_genes']} signal-colocalization-prioritized genes, "
        f"{summary['n_primary_region_genes']} have a primary BrainSEQ regional counterpart; "
        f"{summary['n_measured_pair_regions']} pair–region combinations in {summary['n_measured_genes']} genes "
        f"are measurable with existing counts. Remaining regional nominations are retained in the coverage "
        f"ledger rather than treated as negative. Primary mappings include BA9–DLPFC, hippocampus–hippocampus "
        f"and generic cortex–DLPFC (an anatomical proxy); adjacent and secondary regions are separate sensitivities. "
        f"Counts are pooled across haplotypes including unassigned fragments, without allelic-significance or "
        f"heterozygosity selection. Support is for the pair-discriminating junction models, not necessarily direct "
        f"observation of the particular colocalizing junction or full-length transcript. The measurements reuse "
        f"BrainSEQ reads and are complementary corroboration, not independent phenotype replication.\n")
    (out / "FIGURE_CAPTION.md").write_text(caption)
    (out / "RESULTS.md").write_text(
        "# Junction-pair corroboration\n\n" + caption + "\n```json\n" + json.dumps(summary, indent=2) + "\n```\n")
    scripts = [Path(__file__).resolve(), root / "manuscript/_h/junction_pair_corroboration_figure.R",
               root / "06_switch_mechanism/_h/02g.junction_pair_corroboration.sh",
               root / "tests/test_junction_pair_corroboration.py",
               root / "06_switch_mechanism/JUNCTION_PAIR_CORROBORATION.md"]
    for p in scripts:
        if p.exists(): inputs.paths.add(p)
    versions = {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "pyarrow", "scipy"]}
    git = lambda *a: subprocess.check_output(["git", "-C", str(root), *a], text=True).strip()
    manifest = {"parameters": params, "python": platform.python_version(), "packages": versions,
                "git_commit": git("rev-parse", "HEAD"), "git_status": git("status", "--short"),
                "inputs_and_scripts": [{"path": str(p.relative_to(root)), "bytes": p.stat().st_size,
                                         "sha256": sha256(p)} for p in sorted(inputs.paths)],
                "elapsed_seconds": round(time.monotonic() - start, 2)}
    (out / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    print(f"Wrote {out}; {manifest['elapsed_seconds']:.1f} seconds", flush=True)


if __name__ == "__main__":
    main()
