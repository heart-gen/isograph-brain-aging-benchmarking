"""Stage 3: independent RBP-binding support for the predicted co-switch regulons.

`rbp_regulon.py` calls candidate regulons from ATtRACT *motif* gain/loss — a prediction of
altered RBP binding, not measured binding. This module tests those predictions against
independent **ENCODE eCLIP** peaks: for each switch gene we ask whether the RBP actually
binds the *switched* sequence (the exonic interval present in one switch-pair isoform but not
the other) preferentially over the gene's constitutive (shared) exons. The headline unit is
**per RBP** (`_rbp_binding_test`): one two-sided exact McNemar over the RBP's unique nominated
regulon genes, BH-corrected across the testable RBPs — because a gene's binding is
region-invariant, the per-(region x module x RBP) regulon rows are NOT independent and are kept
descriptive only. The switched-vs-constitutive skew is small and near-universal, so this is
binding *capacity* for a subset of factors, not selective regulon occupancy.

Coverage is honest and partial: ENCODE eCLIP profiles ~168 RBPs in HepG2/K562 (not brain), so
only the subset of significant regulon RBPs with an eCLIP experiment is testable (~39 of 130),
and binding is cross-cell-type evidence of physical capacity, not brain-specific occupancy.
Both caveats are stated in the summary.

Subcommands:
  fetch  — download IDR-thresholded eCLIP peaks (GRCh38) for the testable RBPs into
           inputs/raw/rbp_binding/<RBP>.bed.gz (union of cell lines; chrom/start/end).
  run    — per region, overlap peaks with each gene's switched exon intervals (GENCODE v47
           gtf cache) and write rbp_binding_calls.parquet; then the per-RBP support table
           rbp_binding_support.parquet (headline), the descriptive rbp_binding_regulon.parquet,
           and RBP_BINDING_SUMMARY.md.
"""
from __future__ import annotations

import argparse
import gzip
import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from isograph_benchmark.paths import ensure_dir, rel
from isograph_benchmark.real_data.interpret_modules import DEFAULT_GTF_CACHE
from isograph_benchmark.real_data.qtl_anchoring import _bare
from isograph_benchmark.real_data.rbp_regulon import _REGIONS, _gene_tags, _TREE_OF

_BIND_DIR = rel("inputs", "raw", "rbp_binding")
_OUT_DIR = rel("real_data", "_m", "rbp")
_ENCODE = "https://www.encodeproject.org"
# ENCODE eCLIP peak beds carry output_type "peaks" (IDR-reproducible peaks from the
# eCLIP pipeline); we take the GRCh38 bed from each of an RBP's eCLIP experiments.
_PEAK_OUTPUT_TYPE = "peaks"


# --------------------------------------------------------------------------- #
# fetch
# --------------------------------------------------------------------------- #
def _md_table(df: pd.DataFrame) -> str:
    """GitHub-flavored markdown table without the optional `tabulate` dependency."""
    cols = list(df.columns)
    head = "| " + " | ".join(map(str, cols)) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def _encode_get(path: str) -> dict:
    req = urllib.request.Request(_ENCODE + path, headers={"Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def testable_rbps(regulon_path: Path, q: float = 0.05) -> list[str]:
    rr = pd.read_parquet(regulon_path)
    return sorted(set(rr[rr["q"] <= q]["rbp"]))


def eclip_peak_files(rbp: str) -> list[str]:
    """hrefs of released GRCh38 eCLIP peak beds for an RBP, across its experiments/cell lines."""
    try:
        d = _encode_get(f"/search/?type=Experiment&assay_title=eCLIP&target.label={rbp}"
                        f"&status=released&format=json&limit=all&field=accession")
    except Exception as e:
        print(f"  [{rbp}] experiment search failed: {e}", flush=True)
        return []
    hrefs = []
    for g in d.get("@graph", []):
        try:
            exp = _encode_get(f"/experiments/{g['accession']}/?format=json")
        except Exception as e:
            print(f"  [{rbp}] {g['accession']} fetch failed: {e}", flush=True)
            continue
        for f in exp.get("files", []):
            if (f.get("file_format") == "bed" and f.get("output_type") == _PEAK_OUTPUT_TYPE
                    and f.get("assembly") == "GRCh38" and f.get("status") == "released"
                    and f.get("href")):
                hrefs.append(f["href"])
    return sorted(set(hrefs))


def fetch(regulon_path: Path, out_dir: Path) -> None:
    ensure_dir(out_dir)
    rbps = testable_rbps(regulon_path)
    manifest = []
    for rbp in rbps:
        dst = out_dir / f"{rbp}.bed.gz"
        hrefs = eclip_peak_files(rbp)
        if not hrefs:
            print(f"[{rbp}] no GRCh38 IDR eCLIP peak file — skipping", flush=True)
            continue
        rows = []
        for h in hrefs:
            try:
                raw = urllib.request.urlopen(_ENCODE + h, timeout=300).read()
                text = gzip.decompress(raw).decode() if h.endswith(".gz") else raw.decode()
            except Exception as e:
                print(f"  [{rbp}] download failed {h}: {e}", flush=True)
                continue
            for line in text.splitlines():
                p = line.split("\t")
                if len(p) >= 3:
                    rows.append((p[0], int(p[1]), int(p[2])))
        if not rows:
            continue
        bed = pd.DataFrame(rows, columns=["chrom", "start", "end"]).drop_duplicates()
        bed.to_csv(dst, sep="\t", index=False, header=False, compression="gzip")
        manifest.append({"rbp": rbp, "n_peaks": len(bed), "n_files": len(hrefs)})
        print(f"[{rbp}] {len(bed)} peaks from {len(hrefs)} file(s) -> {dst.name}", flush=True)
    pd.DataFrame(manifest).to_csv(out_dir / "eclip_manifest.tsv", sep="\t", index=False)
    print(f"fetched {len(manifest)}/{len(rbps)} testable RBPs -> {out_dir}", flush=True)


# --------------------------------------------------------------------------- #
# run: switched-interval x peak overlap
# --------------------------------------------------------------------------- #
def _exon_intervals(gtf: pd.DataFrame, tx_ids: set[str]) -> dict[str, list[tuple[int, int]]]:
    e = gtf[(gtf["feature"] == "exon") & (gtf["transcript_id"].isin(tx_ids))]
    out: dict[str, list[tuple[int, int]]] = {}
    for tid, g in e.groupby("transcript_id"):
        out[tid] = list(zip(g["start"].to_numpy(int), g["end"].to_numpy(int)))
    return out


def _switched_intervals(sp: pd.DataFrame, gtf: pd.DataFrame) -> pd.DataFrame:
    """Per (gene) and per switch pair, the genomic intervals covered by exactly ONE isoform
    (`switched` = symmetric difference = the switched sequence) and by BOTH (`constitutive` =
    intersection = the shared sequence, the within-gene binding null). Returns gene_id, chrom,
    start, end, interval_type."""
    tx_ids = set(sp["transcript_id_1"]) | set(sp["transcript_id_2"])
    ex = _exon_intervals(gtf, tx_ids)
    chrom_of = gtf.drop_duplicates("transcript_id").set_index("transcript_id")["chrom"].to_dict()
    rows = []
    for r in sp.itertuples(index=False):
        t1, t2 = r.transcript_id_1, r.transcript_id_2
        if t1 not in ex or t2 not in ex:
            continue
        chrom = chrom_of.get(t1)
        a, b = _merge(ex[t1]), _merge(ex[t2])
        for s, e in _sym_diff(a, b):
            rows.append((r.gene_id, chrom, s, e, "switched"))
        for s, e in _intersect(a, b):
            rows.append((r.gene_id, chrom, s, e, "constitutive"))
    return pd.DataFrame(rows, columns=["gene_id", "chrom", "start", "end", "interval_type"])


def _intersect(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Intervals present in BOTH merged interval lists a and b."""
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo <= hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _merge(iv: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def _sym_diff(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Intervals in exactly one of merged interval lists a,b (a XOR b)."""
    return _subtract(a, b) + _subtract(b, a)


def _subtract(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out = []
    for s, e in a:
        cur = [(s, e)]
        for bs, be in b:
            nxt = []
            for cs, ce in cur:
                if be < cs or bs > ce:
                    nxt.append((cs, ce))
                else:
                    if cs < bs:
                        nxt.append((cs, bs - 1))
                    if ce > be:
                        nxt.append((be + 1, ce))
            cur = nxt
        out.extend(cur)
    return out


def _overlap_flags(intervals: pd.DataFrame, peaks_by_chrom: dict) -> np.ndarray:
    """Boolean per interval row: does any peak overlap it? (peaks pre-sorted per chrom)."""
    flags = np.zeros(len(intervals), dtype=bool)
    for i, (chrom, s, e) in enumerate(zip(intervals["chrom"], intervals["start"], intervals["end"])):
        arr = peaks_by_chrom.get(chrom)
        if arr is None:
            continue
        starts, ends = arr
        # any peak with start<=e and end>=s
        lo = np.searchsorted(starts, e, side="right")
        if lo == 0:
            continue
        if (ends[:lo] >= s).any():
            flags[i] = True
    return flags


def _overlap_bp(intervals: pd.DataFrame, peaks_by_chrom: dict) -> np.ndarray:
    """Peak-covered bp per interval row.

    The binary any-overlap flag makes a long switched interval beat a short constitutive one
    by construction: the longer the sequence, the likelier *some* peak touches it. Covered bp
    divided by interval kb gives a binding density that is comparable across intervals of
    different length, which is what the switched-vs-constitutive contrast needs.
    """
    covered = np.zeros(len(intervals), dtype=np.int64)
    for i, (chrom, s, e) in enumerate(zip(intervals["chrom"], intervals["start"],
                                          intervals["end"])):
        arr = peaks_by_chrom.get(chrom)
        if arr is None:
            continue
        starts, ends = arr
        lo = np.searchsorted(starts, e, side="right")
        if lo == 0:
            continue
        ps, pe = starts[:lo], ends[:lo]
        hit = pe >= s
        if not hit.any():
            continue
        # Peaks may overlap each other; union them so shared bases are not double-counted.
        a = np.maximum(ps[hit], s)
        b = np.minimum(pe[hit], e)
        order = np.argsort(a)
        a, b = a[order], b[order]
        total, cur_s, cur_e = 0, a[0], b[0]
        for x, y in zip(a[1:], b[1:]):
            if x > cur_e + 1:
                total += cur_e - cur_s + 1
                cur_s, cur_e = x, y
            else:
                cur_e = max(cur_e, y)
        covered[i] = total + (cur_e - cur_s + 1)
    return covered


def _load_peaks(rbp: str, bind_dir: Path) -> dict | None:
    p = bind_dir / f"{rbp}.bed.gz"
    if not p.exists():
        return None
    bed = pd.read_csv(p, sep="\t", header=None, names=["chrom", "start", "end"])
    by = {}
    for chrom, g in bed.groupby("chrom"):
        g = g.sort_values("start")
        by[chrom] = (g["start"].to_numpy(int), g["end"].to_numpy(int))
    return by


def _binding_calls_region(region_tree: str, region: str, gtf: pd.DataFrame,
                          rbps: list[str], bind_dir: Path, fdr: float) -> pd.DataFrame:
    art = rel("real_data", region_tree, region, "_m", "isograph_vae")
    sp_path = art / "module_interpret" / "structure_switch_pairs.parquet"
    if not sp_path.exists():
        return pd.DataFrame()
    tag, pool_source = _gene_tags(art, fdr)
    if tag.empty:
        return pd.DataFrame()
    sp = pd.read_parquet(sp_path)
    sp["gene"] = _bare(sp["gene_id"])
    sw = _switched_intervals(sp, gtf)
    if sw.empty:
        return pd.DataFrame()
    sw["gene"] = _bare(sw["gene_id"])
    tag = tag.groupby("gene").agg(go_invisible=("go_invisible", "max"),
                                  module_id=("module_id", "first")).reset_index()
    sw = sw.merge(tag, on="gene", how="inner")
    if sw.empty:
        return pd.DataFrame()

    rows = []
    for rbp in rbps:
        peaks = _load_peaks(rbp, bind_dir)
        if peaks is None:
            continue
        sw = sw.assign(_hit=_overlap_flags(sw, peaks),
                       _bp=_overlap_bp(sw, peaks),
                       _len=sw["end"] - sw["start"] + 1)
        # per gene: is the RBP bound at the SWITCHED interval and/or the CONSTITUTIVE interval?
        wide = (sw.pivot_table(index=["gene", "module_id", "go_invisible"],
                               columns="interval_type", values="_hit", aggfunc="max")
                  .reindex(columns=["switched", "constitutive"], fill_value=False)
                  .reset_index())
        wide["bound_switched"] = wide["switched"] == True   # noqa: E712 (NaN-safe cast)
        wide["bound_constitutive"] = wide["constitutive"] == True  # noqa: E712
        # Length-normalised binding density alongside the binary flag.
        span = (sw.pivot_table(index="gene", columns="interval_type",
                               values=["_bp", "_len"], aggfunc="sum", fill_value=0))
        for kind in ("switched", "constitutive"):
            bp = (span[("_bp", kind)] if ("_bp", kind) in span.columns
                  else pd.Series(0, index=span.index))
            ln = (span[("_len", kind)] if ("_len", kind) in span.columns
                  else pd.Series(0, index=span.index))
            wide[f"bp_{kind}"] = wide["gene"].map(bp).fillna(0).astype(int)
            wide[f"len_{kind}"] = wide["gene"].map(ln).fillna(0).astype(int)
            wide[f"density_{kind}"] = np.where(
                wide[f"len_{kind}"] > 0,
                wide[f"bp_{kind}"] / (wide[f"len_{kind}"] / 1000.0), np.nan)
        wide["rbp"] = rbp
        wide["region"] = region
        wide["pool_source"] = pool_source
        rows.append(wide[["region", "gene", "module_id", "go_invisible", "rbp",
                          "bound_switched", "bound_constitutive",
                          "bp_switched", "bp_constitutive",
                          "len_switched", "len_constitutive",
                          "density_switched", "density_constitutive", "pool_source"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


N_BOOT = 2000
SEED = 13
_Z95 = 1.959963984540054


def _family_note(rbp_support: pd.DataFrame) -> str:
    """How many independent motif families the supported RBPs represent.

    "17 of 39 RBPs" overstates the evidence when several of those RBPs share a
    near-identical motif, so the family-collapsed count is reported next to it.
    """
    fam_path = rel("real_data", "_m", "rbp", "rbp_motif_families.parquet")
    if not fam_path.exists() or rbp_support.empty:
        return ""
    fam = pd.read_parquet(fam_path, columns=["rbp", "family_id"]).drop_duplicates()
    by_rbp = fam.groupby("rbp")["family_id"].apply(set)
    sup = rbp_support.loc[rbp_support["binding_supported"], "rbp"]
    tested = rbp_support["rbp"]
    n_sup_fam = len(set().union(*[by_rbp.get(r, set()) for r in sup])) if len(sup) else 0
    n_test_fam = len(set().union(*[by_rbp.get(r, set()) for r in tested])) if len(tested) else 0
    if not n_test_fam:
        return ""
    return (f", which represent **{n_sup_fam}** of **{n_test_fam}** independent motif "
            f"families (see `rbp_motif_families.py`)")


def _matched_or(b: int, c: int) -> dict[str, float]:
    """Haldane-corrected matched odds ratio (b/c) with a Wald CI on the log scale.

    The 0.5 correction keeps the estimate and its interval finite when one discordant cell
    is empty, which happens for the RBPs with few nominated genes.
    """
    if b + c == 0:
        return {"matched_or": np.nan, "matched_or_ci_low": np.nan,
                "matched_or_ci_high": np.nan}
    bh, ch = b + 0.5, c + 0.5
    log_or = float(np.log(bh / ch))
    se = float(np.sqrt(1.0 / bh + 1.0 / ch))
    return {"matched_or": float(np.exp(log_or)),
            "matched_or_ci_low": float(np.exp(log_or - _Z95 * se)),
            "matched_or_ci_high": float(np.exp(log_or + _Z95 * se))}


def _density_ratio(g: pd.DataFrame, n_boot: int, seed: int) -> dict[str, float]:
    """Peak-bp-per-kb in switched vs constitutive intervals, with a gene-level bootstrap CI.

    Pooled over the RBP's nominated genes (total covered bp over total interval kb) rather
    than averaged over per-gene ratios, so genes with tiny intervals cannot dominate. Genes
    are the resampling unit because a gene's switched and constitutive intervals are not
    independent of each other.
    """
    bp_s = g["bp_switched"].to_numpy(float)
    bp_c = g["bp_constitutive"].to_numpy(float)
    ln_s = g["len_switched"].to_numpy(float) / 1000.0
    ln_c = g["len_constitutive"].to_numpy(float) / 1000.0

    def ratio(idx):
        s = bp_s[idx].sum() / ln_s[idx].sum() if ln_s[idx].sum() > 0 else np.nan
        c = bp_c[idx].sum() / ln_c[idx].sum() if ln_c[idx].sum() > 0 else np.nan
        return s, c, (s / c if c and np.isfinite(c) and c > 0 else np.nan)

    all_idx = np.arange(len(g))
    d_s, d_c, obs = ratio(all_idx)
    lo = hi = np.nan
    if len(g) >= 3 and np.isfinite(obs):
        rng = np.random.default_rng(seed)
        draws = np.empty(n_boot)
        for i in range(n_boot):
            draws[i] = ratio(rng.integers(0, len(g), len(g)))[2]
        draws = draws[np.isfinite(draws)]
        if draws.size > 1:
            lo, hi = (float(v) for v in np.quantile(draws, [0.025, 0.975]))
    return {"density_switched": float(d_s) if np.isfinite(d_s) else np.nan,
            "density_constitutive": float(d_c) if np.isfinite(d_c) else np.nan,
            "density_ratio": float(obs) if np.isfinite(obs) else np.nan,
            "density_ci_low": lo, "density_ci_high": hi}


def _rbp_binding_test(calls: pd.DataFrame, regulon_path: Path, fdr: float) -> pd.DataFrame:
    """Independence-respecting per-RBP binding support — the headline unit.

    A gene's switched/constitutive binding for an RBP is (near-)invariant across the regions
    and modules the gene appears in, so the per-(region, module, RBP) regulon rows are NOT
    independent and must not be pooled into one FDR family (that would count the same
    gene x RBP binding fact once per region). We instead restrict to the genes actually
    nominated as that RBP's regulon members (member of a module where the RBP is a significant
    motif regulon), collapse to the unique gene x RBP unit (bound if bound in ANY of the gene's
    region-specific switch definitions), and run ONE two-sided exact McNemar per RBP over its
    unique nominated genes, BH-correcting across the testable RBPs. `binding_supported` =
    preferential (switched > constitutive bound rate) AND BH q<=0.05."""
    from scipy.stats import binomtest
    if calls.empty:
        return pd.DataFrame()
    rr = pd.read_parquet(regulon_path)
    sig = rr[rr["q"] <= fdr][["region", "module_id", "rbp"]].drop_duplicates()
    nominated = calls.merge(sig, on=["region", "module_id", "rbp"], how="inner")
    if nominated.empty:
        return pd.DataFrame()
    agg = {"bound_switched": ("bound_switched", "max"),
           "bound_constitutive": ("bound_constitutive", "max")}
    for c in ("bp_switched", "bp_constitutive", "len_switched", "len_constitutive"):
        if c in nominated.columns:
            agg[c] = (c, "max")
    per_gene = nominated.groupby(["rbp", "gene"]).agg(**agg).reset_index()
    has_bp = "bp_switched" in per_gene.columns
    out = []
    for rbp, g in per_gene.groupby("rbp"):
        n = len(g)
        n_sw = int(g["bound_switched"].sum())
        n_co = int(g["bound_constitutive"].sum())
        b = int((g["bound_switched"] & ~g["bound_constitutive"]).sum())   # switched-only
        c = int((~g["bound_switched"] & g["bound_constitutive"]).sum())   # constitutive-only
        p = binomtest(b, b + c, 0.5, alternative="two-sided").pvalue if (b + c) > 0 else np.nan
        row = {"rbp": rbp, "n_genes": n,
               "rate_switched": round(n_sw / n, 3), "rate_constitutive": round(n_co / n, 3),
               "rate_diff": round(n_sw / n - n_co / n, 3),
               "n_switched_only": b, "n_constitutive_only": c, "mcnemar_p": p}
        # Haldane-corrected matched odds ratio with a Wald CI, so the McNemar p is
        # accompanied by an effect size rather than a bare significance call.
        row.update(_matched_or(b, c))
        if has_bp:
            row.update(_density_ratio(g, n_boot=N_BOOT, seed=SEED))
        out.append(row)
    res = pd.DataFrame(out)
    ok = res["mcnemar_p"].notna()
    res.loc[ok, "mcnemar_fdr"] = multipletests(res.loc[ok, "mcnemar_p"], method="fdr_bh")[1]
    res["preferential"] = res["rate_switched"] > res["rate_constitutive"]
    res["binding_supported"] = res["preferential"] & (res["mcnemar_fdr"] <= 0.05)
    return res.sort_values("mcnemar_fdr").reset_index(drop=True)


def _regulon_binding_join(calls: pd.DataFrame, regulon_path: Path, fdr: float,
                          rbp_support: pd.DataFrame) -> pd.DataFrame:
    """Annotate each significant module x RBP regulon with a DESCRIPTIVE within-gene binding
    contrast: switched vs constitutive bound rates, the discordant counts, and a two-sided
    exact McNemar p. These per-regulon rows are NOT an independent FDR family (the same
    gene x RBP binding recurs across regions), so the significance gate is inherited from the
    per-RBP `_rbp_binding_test`: a regulon is `binding_supported` when it shows the preferential
    direction AND its RBP has independent binding support."""
    from scipy.stats import binomtest
    rr = pd.read_parquet(regulon_path)
    sig = rr[(rr["q"] <= fdr)].copy()
    if calls.empty or sig.empty:
        return pd.DataFrame()
    out = []
    for (region, mid, rbp), g in calls.groupby(["region", "module_id", "rbp"]):
        n = len(g)
        n_sw = int(g["bound_switched"].sum())
        n_co = int(g["bound_constitutive"].sum())
        # discordant genes: bound at exactly one interval type
        b = int((g["bound_switched"] & ~g["bound_constitutive"]).sum())   # switched-only
        c = int((~g["bound_switched"] & g["bound_constitutive"]).sum())   # constitutive-only
        p = binomtest(b, b + c, 0.5, alternative="two-sided").pvalue if (b + c) > 0 else np.nan
        out.append({"region": region, "module_id": mid, "rbp": rbp, "n_genes": n,
                    "rate_switched": round(n_sw / n, 3), "rate_constitutive": round(n_co / n, 3),
                    "n_switched_only": b, "n_constitutive_only": c,
                    "mcnemar_p": p, "go_invisible": bool(g["go_invisible"].max())})
    res = pd.DataFrame(out)
    if res.empty:
        return res
    res["preferential"] = res["rate_switched"] > res["rate_constitutive"]
    supported = (set(rbp_support.loc[rbp_support["binding_supported"], "rbp"])
                 if not rbp_support.empty else set())
    fdr_map = (rbp_support.set_index("rbp")["mcnemar_fdr"].to_dict()
               if not rbp_support.empty else {})
    res["rbp_binding_fdr"] = res["rbp"].map(fdr_map)
    res["binding_supported"] = res["preferential"] & res["rbp"].isin(supported)
    merged = sig.merge(res, on=["region", "module_id", "rbp"], how="inner", suffixes=("", "_bind"))
    return merged


def run(fdr: float, bind_dir: Path, gtf_cache: Path) -> None:
    regulon_path = _OUT_DIR / "rbp_regulon.parquet"
    rbps = testable_rbps(regulon_path, fdr)
    gtf = pd.read_parquet(gtf_cache, columns=["transcript_id", "gene_id", "chrom", "feature",
                                              "start", "end"])
    all_calls = []
    for tree, region in _REGIONS:
        c = _binding_calls_region(tree, region, gtf, rbps, bind_dir, fdr)
        if not c.empty:
            all_calls.append(c)
            print(f"[{tree}/{region}] {c['gene'].nunique()} genes x {c['rbp'].nunique()} RBPs; "
                  f"switched-bound {int(c['bound_switched'].sum())} / "
                  f"constitutive-bound {int(c['bound_constitutive'].sum())}", flush=True)
    calls = pd.concat(all_calls, ignore_index=True) if all_calls else pd.DataFrame()
    out = ensure_dir(_OUT_DIR)
    calls.to_parquet(out / "rbp_binding_calls.parquet", index=False, compression="zstd")
    rbp_support = _rbp_binding_test(calls, regulon_path, fdr)
    rbp_support.to_parquet(out / "rbp_binding_support.parquet", index=False, compression="zstd")
    joined = _regulon_binding_join(calls, regulon_path, fdr, rbp_support)
    joined.to_parquet(out / "rbp_binding_regulon.parquet", index=False, compression="zstd")
    _write_summary(calls, joined, rbp_support, rbps, bind_dir, out)


def _write_summary(calls, joined, rbp_support, testable, bind_dir, out) -> None:
    have = sorted({p.stem.replace(".bed", "") for p in bind_dir.glob("*.bed.gz")})
    lines = ["# RBP binding-evidence support for predicted co-switch regulons", ""]
    lines += [f"- Testable regulon RBPs (sig ∩ ENCODE eCLIP): **{len(testable)}**; "
              f"peak files present: **{len(have)}**.", ""]
    if not rbp_support.empty:
        nt = len(rbp_support)
        npref = int(rbp_support["preferential"].sum())
        nsup = int(rbp_support["binding_supported"].sum())
        med = rbp_support["rate_diff"].median()
        lines += [
            "## Per-RBP binding support (independence-respecting headline)", "",
            "Unit = one **two-sided exact McNemar per RBP** over its **unique nominated regulon "
            "genes** (a gene member of a module where the RBP is a significant motif regulon), "
            "deduplicated across regions because a gene's switched/constitutive binding is "
            f"region-invariant; BH-corrected across the {nt} testable RBPs. The earlier "
            "per-(region×module×RBP) regulon rows are **not** an independent FDR family — the "
            "same gene×RBP binding recurs across regions — so they are reported descriptively "
            "below, not as the significance count.", "",
            f"- RBPs with preferential switched-interval binding: **{npref}/{nt}**; "
            f"**binding-supported** (preferential AND BH q≤0.05): **{nsup}/{nt}**"
            f"{_family_note(rbp_support)}.",
            f"- Median switched−constitutive bound-rate gap across RBPs: **{med:.3f}** — small "
            "and near-universal. This is binding *capacity* at alternative vs constitutive "
            "exons, **not** factor-specific occupancy; it does not establish that the predicted "
            "regulon factors selectively bind the switched sequence beyond a generic "
            "alternative-exon skew.", "",
            "### Binding-supported RBPs (switched > constitutive)", "",
            "`matched_or` is the Haldane-corrected discordant-pair odds ratio behind the "
            "McNemar p; `density_ratio` is peak-covered bp per kb in switched vs constitutive "
            "intervals, with a gene-level bootstrap CI. The binary rates alone cannot "
            "distinguish a real preference from a longer switched interval having more chance "
            "to be touched by some peak.", "",
            _md_table(rbp_support[rbp_support["binding_supported"]]
                      [[c for c in ("rbp", "n_genes", "rate_switched", "rate_constitutive",
                                    "rate_diff", "matched_or", "matched_or_ci_low",
                                    "matched_or_ci_high", "density_ratio", "density_ci_low",
                                    "density_ci_high", "mcnemar_fdr")
                        if c in rbp_support.columns]].round(4)),
            "",
        ]
    else:
        lines += ["- No per-RBP binding support (check that `fetch` ran).", ""]
    if not joined.empty:
        nreg = len(joined)
        nreg_sup = int(joined["binding_supported"].sum())
        gi_sup = int(joined[joined["go_invisible"] == True]["binding_supported"].sum())
        lines += [
            "## Regulon annotation (descriptive, not an FDR family)", "",
            f"- Significant module×RBP regulons with a testable RBP: **{nreg}**.",
            f"- Regulons whose RBP is binding-supported and that show the preferential "
            f"direction: **{nreg_sup}**; of those GO-invisible: **{gi_sup}**.", ""]
    lines += ["", "_Metric: the RBP's eCLIP peaks are overlapped with each member gene's SWITCHED "
              "exon interval (present in one switch isoform) vs its CONSTITUTIVE interval (shared "
              "by both); an RBP is scored over its unique nominated regulon genes (two-sided exact "
              "McNemar on discordant genes, BH across RBPs). This within-gene null neutralizes "
              "peak-dense RBPs. Caveats: ENCODE eCLIP is HepG2/K562, not brain — binding capacity, "
              "not brain-specific occupancy; only regulon RBPs with an eCLIP experiment are "
              "testable; the switched-vs-constitutive skew is small and shared by most RBPs, so "
              "this supports a subset of factors' binding *capacity*, not selective regulon "
              "occupancy._"]
    (out / "RBP_BINDING_SUMMARY.md").write_text("\n".join(lines))
    print(f"wrote {out/'RBP_BINDING_SUMMARY.md'}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    pf = sub.add_parser("fetch", help="download eCLIP peaks for testable RBPs")
    pf.add_argument("--bind-dir", type=Path, default=_BIND_DIR)
    pr = sub.add_parser("run", help="overlap switched intervals with peaks + join regulons")
    pr.add_argument("--bind-dir", type=Path, default=_BIND_DIR)
    pr.add_argument("--gtf-cache", type=Path, default=DEFAULT_GTF_CACHE)
    pr.add_argument("--fdr", type=float, default=0.05)
    args = p.parse_args()
    if args.cmd == "fetch":
        fetch(_OUT_DIR / "rbp_regulon.parquet", args.bind_dir)
    else:
        run(args.fdr, args.bind_dir, args.gtf_cache)


if __name__ == "__main__":
    main()
