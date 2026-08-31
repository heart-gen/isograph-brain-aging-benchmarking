"""A2: independent long-read confirmation of IsoGraph aging switch transcript pairs.

IsoGraph discovers co-switching modules from short-read (Salmon) isoform composition.
A reviewer's fair question is whether the transcript switches that define those modules
are *real* at the isoform level, or artefacts of short-read transcript quantification.
This CLI answers that on a fully independent platform, lab, and quantifier:

  **Aguzzoli-Heberle et al., Nat Biotechnol 2024** -- deep Oxford Nanopore (PromethION)
  long-read RNA-seq of aged human dorsolateral prefrontal cortex (Brodmann area 9/46),
  n=12 (6 AD / 6 control), quantified with **Bambu**. Processed transcript-count matrix
  is open on Zenodo (10.5281/zenodo.8180677, ``counts_transcript.txt``).

Region matching: the long-read tissue is DLPFC BA9/46, so the region-matched IsoGraph
switches come from **GTEx frontal_cortex_ba9 + cortex aging switch modules**. (BrainSEQ
DLPFC aging is vacuous -- its switch modules are caudate-only -- so it is not used here.)

Tier-1 scope (this module): confirm, on the long-read Bambu matrix, that
  (i)  the transcripts constituting IsoGraph's consequential (coding) switch pairs are
       independently *detected*  -- transcript-pair existence; and
  (ii) each switch gene is genuinely multi-isoform in long-read and its two dominant
       isoforms are *anti-correlated in usage* across samples -- the platform-independent
       signature of an isoform switch (one isoform up implies the other down).
Both are reported as *rates over all tested switches* (not cherry-picked examples).

What Tier-1 deliberately does NOT claim: a powered re-association of usage with age.
The public matrix carries no per-sample age, and n=12 is under-powered for a directional
age test -- so a signed usage-vs-age concordance is the Tier-1.5 extension (supply a
per-sample age table via ``--sample-meta``) / Tier-2 (raw-read realignment). This module
computes the direction-agnostic switch signature that the open data *can* adjudicate, and
says so.

This is confirmation of discovered structure on an orthogonal modality -- exactly the
frontier-modality anchor the short-read critique asks for -- not a rediscovery.

Outputs land in ``06_switch_mechanism/_m/longread_switch_confirm/``:
  * ``pair_confirmation.parquet`` -- per switch pair: transcripts, detection, mean IF,
    cross-sample usage anti-correlation, passed-through structural consequence flags.
  * ``gene_confirmation.parquet`` -- per switch gene: expressed / multi-isoform /
    top-2 anti-correlated / any consequential pair detected.
  * ``summary.json`` -- per-region + pooled confirmed rates, thresholds, provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data.validate_switch_splicing import (
    _artifact_dir,
    _module_switch_genes,
    _strip_ver,
)

# --------------------------------------------------------------------------- #
# Dataset (Aguzzoli-Heberle 2024, Zenodo 10.5281/zenodo.8180677)
# --------------------------------------------------------------------------- #
_ZENODO_RECORD = "8180677"
_COUNTS_FILE = "counts_transcript.txt"
_COUNTS_URL = f"https://zenodo.org/records/{_ZENODO_RECORD}/files/{_COUNTS_FILE}?download=1"
# Region-matched to the long-read tissue (DLPFC BA9/46). BrainSEQ DLPFC aging is
# vacuous (caudate-only switch modules) so only these GTEx cortical regions are used.
_DEFAULT_REGIONS = ("frontal_cortex_ba9", "cortex")


def _default_data_dir():
    return ensure_dir(rel("inputs", "_m", "longread_aged_dlpfc"))


def _out_dir():
    return ensure_dir(stage_out("mechanism", "longread_switch_confirm"))


# --------------------------------------------------------------------------- #
# Fetch (Tier-1: processed Bambu matrix only -- avoids the ~936 GB raw pull)
# --------------------------------------------------------------------------- #
def _md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _download(url: str, tmp: Path, retries: int = 4) -> None:
    """Fetch ``url`` to ``tmp``, tolerating PSC's flaky TLS to Zenodo.

    urllib's handshake intermittently fails against Zenodo from the cluster
    (``SSL: UNEXPECTED_EOF_WHILE_READING``), while ``curl``/``wget`` (which retry
    the TLS negotiation and support resume) succeed. Try urllib with backoff,
    then fall back to whichever external downloader is available.
    """
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            urllib.request.urlretrieve(url, tmp)  # noqa: S310 (fixed https Zenodo URL)
            return
        except Exception as err:  # noqa: BLE001 -- retry any transport failure
            last_err = err
            print(f"[fetch] urllib attempt {attempt}/{retries} failed: {err}")
            time.sleep(min(2 ** attempt, 30))
    for tool, cmd in (
        ("curl", ["curl", "-fSL", "--retry", "5", "--retry-all-errors",
                  "-C", "-", "-o", str(tmp), url]),
        ("wget", ["wget", "-c", "-O", str(tmp), url]),
    ):
        if shutil.which(tool) is None:
            continue
        print(f"[fetch] falling back to {tool}")
        try:
            subprocess.run(cmd, check=True)  # noqa: S603 -- fixed argv, trusted URL
            return
        except subprocess.CalledProcessError as err:
            last_err = err
            print(f"[fetch] {tool} failed: {err}")
    raise RuntimeError(f"could not download {url}") from last_err


def fetch(data_dir: Path, force: bool = False) -> Path:
    """Download the processed Bambu transcript-count matrix from Zenodo (~92 MB)."""
    data_dir = ensure_dir(data_dir)
    dest = data_dir / _COUNTS_FILE
    if dest.exists() and not force:
        print(f"[fetch] {dest} already present ({dest.stat().st_size/1e6:.1f} MB); "
              "use --force to re-download.")
    else:
        print(f"[fetch] downloading {_COUNTS_URL}\n        -> {dest}")
        tmp = dest.with_suffix(dest.suffix + ".part")
        _download(_COUNTS_URL, tmp)
        tmp.replace(dest)
        print(f"[fetch] done ({dest.stat().st_size/1e6:.1f} MB)")
    prov = {
        "source": "Aguzzoli-Heberle et al., Nat Biotechnol 2024 (DOI 10.1038/s41587-024-02245-9)",
        "zenodo_record": _ZENODO_RECORD,
        "url": _COUNTS_URL,
        "file": _COUNTS_FILE,
        "bytes": dest.stat().st_size,
        "md5": _md5(dest),
        "quantifier": "Bambu (transcript-level counts, intronic reads excluded)",
        "assay": "ONT PromethION cDNA (SQK-PCS111), DLPFC BA9/46, n=12 (6 AD / 6 control)",
    }
    (data_dir / "provenance.json").write_text(json.dumps(prov, indent=2))
    return dest


# --------------------------------------------------------------------------- #
# Long-read matrix -> per-transcript detection + per-gene isoform fractions
# --------------------------------------------------------------------------- #
def _load_longread(counts_path: Path, keep_genes: set[str]):
    """Return (per-transcript table, sample columns) restricted to ``keep_genes``.

    Bambu ``counts_transcript.txt``: TXNAME | GENEID | <one column per sample>. Known
    transcripts carry unversioned ENST ids, novel ones ``BambuTxN`` (kept, since they
    contribute to a gene's usage denominator). GENEID is unversioned ENSG.
    """
    lr = pd.read_csv(counts_path, sep="\t")
    sample_cols = [c for c in lr.columns if c not in ("TXNAME", "GENEID")]
    lr["gene"] = _strip_ver(lr["GENEID"].astype(str))
    lr["tx"] = _strip_ver(lr["TXNAME"].astype(str))
    lr = lr[lr["gene"].isin(keep_genes)].copy()
    lr[sample_cols] = lr[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return lr, sample_cols


def _isoform_fractions(lr: pd.DataFrame, sample_cols: list[str]):
    """Per-sample within-gene isoform fractions; NaN in a sample with zero gene counts."""
    gene_tot = lr.groupby("gene")[sample_cols].transform("sum")
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = lr[sample_cols].to_numpy() / gene_tot.to_numpy()
    frac = pd.DataFrame(frac, index=lr.index, columns=sample_cols)
    return frac


def _tx_stats(lr: pd.DataFrame, frac: pd.DataFrame, sample_cols: list[str],
              min_count: float, min_samples: int) -> pd.DataFrame:
    """Per-transcript long-read summary: detection + mean isoform fraction."""
    counts = lr[sample_cols]
    n_det = (counts >= min_count).sum(axis=1)
    out = pd.DataFrame({
        "gene": lr["gene"].values,
        "tx": lr["tx"].values,
        "total_count": counts.sum(axis=1).values,
        "n_detected_samples": n_det.values,
        "detected": (n_det >= min_samples).values,
        "mean_if": frac.mean(axis=1, skipna=True).values,
    })
    return out


# --------------------------------------------------------------------------- #
# Switch pairs (consequential/coding pairs from switch_consequence)
# --------------------------------------------------------------------------- #
def _switch_pairs(cohort: str, region: str, variant: str, trait: str,
                  switch_genes: set[str], coding_only: bool) -> pd.DataFrame:
    """Transcript pairs from ``switch_consequence/pair_consequence.parquet`` restricted
    to trait-associated switch-module genes (and, by default, coding-consequence pairs --
    the consequential switches the manuscript foregrounds)."""
    path = _artifact_dir(cohort, region, variant, trait) / "switch_consequence" / "pair_consequence.parquet"
    pc = pd.read_parquet(path)
    pc["gene"] = _strip_ver(pc["gene"].astype(str))
    pc = pc[pc["gene"].isin(switch_genes)].copy()
    if coding_only and "coding_consequence" in pc.columns:
        pc = pc[pc["coding_consequence"].fillna(False).astype(bool)].copy()
    pc["t1"] = _strip_ver(pc["transcript_id_1"].astype(str))
    pc["t2"] = _strip_ver(pc["transcript_id_2"].astype(str))
    return pc


# --------------------------------------------------------------------------- #
# Confirmation
# --------------------------------------------------------------------------- #
def _anticorr(frac: pd.DataFrame, idx1: int, idx2: int):
    """Spearman rho of two isoforms' usage across samples; <0 => switch-like."""
    a = frac.loc[idx1].to_numpy(dtype=float)
    b = frac.loc[idx2].to_numpy(dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4 or np.nanstd(a[ok]) == 0 or np.nanstd(b[ok]) == 0:
        return np.nan
    return float(stats.spearmanr(a[ok], b[ok]).correlation)


def confirm_region(cohort: str, region: str, variant: str, trait: str,
                   lr_all: pd.DataFrame, sample_cols: list[str],
                   alpha: float, coding_only: bool,
                   min_count: float, min_samples: int, min_if: float):
    switch_genes = _module_switch_genes(cohort, region, variant, trait, alpha)
    if not switch_genes:
        return None, None, {"region": region, "n_switch_genes": 0,
                            "note": "no trait-associated switch-module genes (vacuous)"}

    pairs = _switch_pairs(cohort, region, variant, trait, switch_genes, coding_only)

    lr = lr_all[lr_all["gene"].isin(switch_genes)].copy()
    frac = _isoform_fractions(lr, sample_cols)
    txs = _tx_stats(lr, frac, sample_cols, min_count, min_samples)
    det = dict(zip(txs["tx"], txs["detected"]))
    mif = dict(zip(txs["tx"], txs["mean_if"]))
    # row index per (gene, tx) for anti-correlation lookups
    lr_row = {(g, t): i for i, (g, t) in enumerate(zip(lr["gene"], lr["tx"]))}
    frac = frac.reset_index(drop=True)

    # ---- pair-level existence + switch-likeness ------------------------------
    rec = []
    for r in pairs.itertuples(index=False):
        d1, d2 = det.get(r.t1, False), det.get(r.t2, False)
        i1, i2 = mif.get(r.t1, np.nan), mif.get(r.t2, np.nan)
        rho = np.nan
        p1, p2 = lr_row.get((r.gene, r.t1)), lr_row.get((r.gene, r.t2))
        if p1 is not None and p2 is not None:
            rho = _anticorr(frac, p1, p2)
        coexpr = bool(d1 and d2 and (i1 >= min_if) and (i2 >= min_if))
        rec.append({
            "region": region, "gene": r.gene, "t1": r.t1, "t2": r.t2,
            "t1_detected": bool(d1), "t2_detected": bool(d2),
            "pair_detected": bool(d1 and d2),
            "t1_mean_if": i1, "t2_mean_if": i2,
            "coexpressed": coexpr, "usage_spearman": rho,
            "switch_like": bool(coexpr and np.isfinite(rho) and rho < 0),
            "coding_consequence": bool(getattr(r, "coding_consequence", False)),
            "nmd_switch": bool(getattr(r, "nmd_switch", False)),
        })
    pair_df = pd.DataFrame(rec)

    # ---- gene-level confirmation --------------------------------------------
    grec = []
    for gene in sorted(switch_genes):
        g_tx = txs[txs["gene"] == gene]
        expressed = bool((g_tx["total_count"] > 0).any())
        top2 = g_tx.sort_values("mean_if", ascending=False).head(2)
        multi_iso = bool(len(top2) == 2 and top2["mean_if"].iloc[1] >= min_if)
        top2_rho = np.nan
        if multi_iso:
            r1 = lr_row.get((gene, top2["tx"].iloc[0]))
            r2 = lr_row.get((gene, top2["tx"].iloc[1]))
            if r1 is not None and r2 is not None:
                top2_rho = _anticorr(frac, r1, r2)
        gp = pair_df[pair_df["gene"] == gene] if len(pair_df) else pair_df
        any_pair_detected = bool(len(gp) and gp["pair_detected"].any())
        any_pair_switchlike = bool(len(gp) and gp["switch_like"].any())
        grec.append({
            "region": region, "gene": gene, "expressed_longread": expressed,
            "multi_isoform": multi_iso, "top2_usage_spearman": top2_rho,
            "top2_anticorr": bool(np.isfinite(top2_rho) and top2_rho < 0),
            "n_switch_pairs": int(len(gp)),
            "any_pair_detected": any_pair_detected,
            "any_pair_switchlike": any_pair_switchlike,
            # confirmed = independently a real, switch-capable isoform locus
            "confirmed": bool(expressed and multi_iso
                              and np.isfinite(top2_rho) and top2_rho < 0),
        })
    gene_df = pd.DataFrame(grec)

    expr = gene_df[gene_df["expressed_longread"]]
    summ = {
        "region": region,
        "n_switch_genes": int(len(gene_df)),
        "n_genes_expressed_longread": int(len(expr)),
        "n_genes_multi_isoform": int(gene_df["multi_isoform"].sum()),
        "n_genes_confirmed": int(gene_df["confirmed"].sum()),
        "confirmed_rate_of_expressed": (float(gene_df["confirmed"].sum() / len(expr))
                                        if len(expr) else float("nan")),
        "n_switch_pairs": int(len(pair_df)),
        "n_pairs_detected": int(pair_df["pair_detected"].sum()) if len(pair_df) else 0,
        "pair_detection_rate": (float(pair_df["pair_detected"].mean())
                                if len(pair_df) else float("nan")),
        "n_pairs_switch_like": int(pair_df["switch_like"].sum()) if len(pair_df) else 0,
        "pair_switch_like_rate": (float(pair_df["switch_like"].mean())
                                  if len(pair_df) else float("nan")),
        "coding_pairs_only": coding_only,
    }
    return pair_df, gene_df, summ


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _run_confirm(args) -> None:
    data_dir = Path(args.data_dir) if args.data_dir else _default_data_dir()
    counts_path = data_dir / _COUNTS_FILE
    if not counts_path.exists():
        counts_path = fetch(data_dir, force=False)

    regions = tuple(args.regions) if args.regions else _DEFAULT_REGIONS
    # union of switch genes across regions -> one restricted long-read load
    all_switch = set()
    for region in regions:
        all_switch |= _module_switch_genes(args.cohort, region, args.variant, args.trait, args.alpha)
    print(f"[confirm] {len(all_switch)} region-matched switch genes across {regions}")
    lr_all, sample_cols = _load_longread(counts_path, all_switch)
    print(f"[confirm] long-read matrix: {len(lr_all)} transcript rows, {len(sample_cols)} samples")

    pair_frames, gene_frames, summaries = [], [], []
    for region in regions:
        pdf, gdf, summ = confirm_region(
            args.cohort, region, args.variant, args.trait, lr_all, sample_cols,
            args.alpha, args.coding_only, args.min_count, args.min_samples, args.min_if)
        summaries.append(summ)
        if pdf is not None:
            pair_frames.append(pdf)
            gene_frames.append(gdf)
        cr = summ.get("confirmed_rate_of_expressed", float("nan"))
        pr = summ.get("pair_detection_rate", float("nan"))
        sr = summ.get("pair_switch_like_rate", float("nan"))
        print(f"[{args.cohort}/{region}/{args.trait}] switch genes {summ['n_switch_genes']} | "
              f"confirmed {summ.get('n_genes_confirmed', 0)}/"
              f"{summ.get('n_genes_expressed_longread', 0)} expressed "
              f"(rate {cr:.3f}) | pairs {summ.get('n_switch_pairs', 0)} "
              f"detected-rate {pr:.3f} switch-like-rate {sr:.3f}")

    out = _out_dir()
    if pair_frames:
        pd.concat(pair_frames, ignore_index=True).to_parquet(out / "pair_confirmation.parquet")
        pd.concat(gene_frames, ignore_index=True).to_parquet(out / "gene_confirmation.parquet")

    # pooled over non-vacuous regions
    gf = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    pf = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    pooled = None
    if len(gf):
        expr = gf[gf["expressed_longread"]]
        pooled = {
            "regions": [s["region"] for s in summaries if s.get("n_switch_genes")],
            "n_genes_confirmed": int(gf["confirmed"].sum()),
            "n_genes_expressed_longread": int(len(expr)),
            "confirmed_rate_of_expressed": (float(gf["confirmed"].sum() / len(expr))
                                            if len(expr) else float("nan")),
            "pair_detection_rate": (float(pf["pair_detected"].mean()) if len(pf) else float("nan")),
            "pair_switch_like_rate": (float(pf["switch_like"].mean()) if len(pf) else float("nan")),
        }
    summary = {
        "analysis": "A2 long-read (Bambu, ONT DLPFC BA9/46) confirmation of IsoGraph "
                    "aging switch transcript pairs",
        "tier": "Tier-1 (processed matrix): isoform-level existence + switch-like "
                "anti-correlated usage; signed usage-vs-age is Tier-1.5/Tier-2",
        "dataset": "Aguzzoli-Heberle 2024 NBT; Zenodo 8180677; n=12 (6 AD/6 control)",
        "cohort": args.cohort, "trait": args.trait, "variant": args.variant,
        "thresholds": {"alpha": args.alpha, "min_count": args.min_count,
                       "min_samples": args.min_samples, "min_if": args.min_if,
                       "coding_pairs_only": args.coding_only},
        "per_region": summaries,
        "pooled": pooled,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[confirm] wrote {out}/summary.json")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pf = sub.add_parser("fetch", help="download the Bambu long-read matrix from Zenodo")
    pf.add_argument("--data-dir", default=None)
    pf.add_argument("--force", action="store_true")

    pc = sub.add_parser("confirm", help="confirm switch pairs against the long-read matrix")
    pc.add_argument("--data-dir", default=None)
    pc.add_argument("--cohort", choices=["gtex", "brainseq"], default="gtex")
    pc.add_argument("--region", action="append", dest="regions",
                    help="repeatable; default = region-matched frontal_cortex_ba9 + cortex")
    pc.add_argument("--trait", choices=["age", "dx"], default="age")
    pc.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    pc.add_argument("--alpha", type=float, default=0.05,
                    help="module trait-association FDR for switch-gene selection")
    pc.add_argument("--coding-only", dest="coding_only", action="store_true", default=True,
                    help="restrict to coding-consequence switch pairs (default on)")
    pc.add_argument("--all-pairs", dest="coding_only", action="store_false",
                    help="use all transcript pairs, not just coding-consequence ones")
    pc.add_argument("--min-count", type=float, default=5.0,
                    help="min Bambu count for a transcript to be 'detected' in a sample")
    pc.add_argument("--min-samples", type=int, default=3,
                    help="min samples meeting --min-count for transcript detection")
    pc.add_argument("--min-if", type=float, default=0.05,
                    help="min mean isoform fraction for an isoform to count as used")

    args = p.parse_args()
    if args.cmd == "fetch":
        fetch(Path(args.data_dir) if args.data_dir else _default_data_dir(), force=args.force)
    else:
        _run_confirm(args)


if __name__ == "__main__":
    main()
