"""Biology gate for phenotype-significant, GO-invisible IsoGraph switch modules.

Some IsoGraph switch modules are phenotype-associated yet carry no canonical GO:BP
enrichment. This script tests whether those modules are coherent isoform-switch
biology or noise, using only the saved module_enrichment + module_interpret outputs
(no re-fit). For each phenotype-significant module it reports whether members carry
real isoform switches (anticorrelated transcript pairs) and how functionally
consequential those switches are, contrasts GO-invisible vs GO-visible modules
against a pooled background, and names the top switch-driver genes.

Per analysis it writes, under <artifact-parent>/_m/:
  go_invisible_gate.parquet — one row per phenotype-significant module: size,
      phenotype fdr, go_invisible flag, switch metrics, functional-consequence
      fractions, and named top switch drivers.
  GO_INVISIBLE_GATE.md — the verdict writeup.

Reads:
  module_enrichment/isograph_modules.parquet (module selection: pheno_fdr, n_go_terms)
  isograph_vae[...]/module_interpret/<module>/{gene_driver,transcript_polarity}.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.interpret_modules import DEFAULT_GTF_PATH
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

FUNCTIONAL_COLS = ("cds_changed", "coding_status_change", "biotype_switch", "utr_changed")
_ABS_R = 0.3
_Q = 0.05
_SYMBOL_CACHE = stage_out("tmp", "gene_id_symbol.parquet")


def gene_symbol_map(gtf_path: Path, cache: Path = _SYMBOL_CACHE) -> dict[str, str]:
    """Unversioned gene_id -> gene_name from the GTF gene lines (cached)."""
    if cache.exists():
        c = pd.read_parquet(cache)
        return dict(zip(c["gene_id"], c["gene_name"]))
    rows = []
    if gtf_path and gtf_path.exists():
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                attrs = f[8]
                gid = attrs.split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
                name = attrs.split('gene_name "', 1)[1].split('"', 1)[0] if 'gene_name "' in attrs else gid
                rows.append((gid, name))
    out = pd.DataFrame(rows, columns=["gene_id", "gene_name"]).drop_duplicates("gene_id")
    if not out.empty:
        out.to_parquet(ensure_dir(cache.parent) / cache.name, index=False, compression="zstd")
    return dict(zip(out["gene_id"], out["gene_name"]))


def _sig_switch_tx(interpret_dir: Path, module_id: str) -> pd.DataFrame:
    t = pd.read_parquet(interpret_dir / module_id / "transcript_polarity_table.parquet")
    t["absr"] = pd.to_numeric(t["r"], errors="coerce").abs()
    t["ss"] = pd.to_numeric(t["switch_strength"], errors="coerce")
    t["q"] = pd.to_numeric(t["qvalue"], errors="coerce")
    return t


def _functional_fracs(sig: pd.DataFrame) -> dict[str, float]:
    out = {}
    for col in FUNCTIONAL_COLS:
        v = sig[col].dropna() if col in sig.columns else pd.Series(dtype=float)
        out[f"frac_{col}"] = round(float(v.mean()), 3) if len(v) else np.nan
    return out


def _real_switch_genes(t: pd.DataFrame) -> int:
    """Genes with both an up- and a down-correlated transcript at q<0.05 (true DTU)."""
    q = t[t["q"] < _Q]
    if q.empty:
        return 0
    r = pd.to_numeric(q["r"], errors="coerce")
    g = q.assign(_pos=r > 0, _neg=r < 0).groupby("gene_id")[["_pos", "_neg"]].any()
    return int((g["_pos"] & g["_neg"]).sum())


def _top_switch_genes(t: pd.DataFrame, symbols: dict[str, str], n: int) -> str:
    sig = t[(t["absr"] >= _ABS_R) & (t["q"] < _Q)]
    if sig.empty:
        return ""
    gg = sig.groupby("gene_id").agg(ss=("ss", "max"), mr=("absr", "max"))
    gg["score"] = gg["ss"] * gg["mr"]
    top = gg.sort_values("score", ascending=False).head(n).index
    return ", ".join(symbols.get(g.split(".")[0], g.split(".")[0]) for g in top)


def run_gate(analysis: str, region: str | None, variant: str, fdr: float, top_n: int,
             gtf_path: Path | None) -> pd.DataFrame:
    iso_dir = _artifact_dir(analysis, region, variant)
    enrich = pd.read_parquet(iso_dir.parent / "module_enrichment" / "isograph_modules.parquet")
    enrich["module_id"] = enrich["module_id"].astype(str)
    interpret_dir = iso_dir / "module_interpret"

    disease = enrich[enrich["pheno_fdr"] <= fdr].copy()
    disease["go_invisible"] = disease["n_go_terms"] == 0
    symbols = gene_symbol_map(gtf_path) if gtf_path else {}

    # pooled background: every interpreted module's significant switch transcripts
    bg = pd.concat(
        [_sig_switch_tx(interpret_dir, d.name) for d in interpret_dir.iterdir() if d.is_dir()],
        ignore_index=True,
    )
    bg_sig = bg[(bg["absr"] >= _ABS_R) & (bg["q"] < _Q)]
    background = {"module": "_background", "go_invisible": None, **_functional_fracs(bg_sig),
                 "n_sig_switch_tx": int(len(bg_sig))}

    rows = []
    for _, m in disease.sort_values(["go_invisible", "pheno_fdr"]).iterrows():
        mid = m["module_id"]
        if not (interpret_dir / mid).is_dir():
            continue
        t = _sig_switch_tx(interpret_dir, mid)
        sig = t[(t["absr"] >= _ABS_R) & (t["q"] < _Q)]
        rows.append({
            "module": mid, "n_genes": int(m["n_genes"]), "pheno_fdr": round(float(m["pheno_fdr"]), 4),
            "go_invisible": bool(m["go_invisible"]), "n_go_terms": int(m["n_go_terms"]),
            "genes_with_real_switch": _real_switch_genes(t),
            "max_switch_strength": round(float(t["ss"].max()), 3),
            "n_sig_switch_tx": int(len(sig)), **_functional_fracs(sig),
            "top_switch_genes": _top_switch_genes(t, symbols, top_n),
        })
    gate = pd.DataFrame(rows)

    out_dir = ensure_dir(iso_dir.parent)
    gate.to_parquet(out_dir / "go_invisible_gate.parquet", index=False, compression="zstd")
    (out_dir / "go_invisible_gate_background.json").write_text(json.dumps(background, indent=2))
    _write_report(out_dir, analysis, region, fdr, gate, background)
    return gate


def _markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    head = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def _write_report(out_dir: Path, analysis: str, region: str | None, fdr: float,
                  gate: pd.DataFrame, background: dict) -> None:
    inv = gate[gate["go_invisible"]]
    vis = gate[~gate["go_invisible"]]
    bg = ", ".join(f"{c.replace('frac_', '')} {background[c]}" for c in
                   [f"frac_{x}" for x in FUNCTIONAL_COLS] if not pd.isna(background.get(c)))
    lines = [
        f"# Biology gate — GO-invisible switch modules ({analysis}"
        + (f"/{region}" if region else "") + ")",
        "",
        f"Phenotype-significant IsoGraph switch modules (pheno_fdr <= {fdr}): "
        f"**{len(gate)}** — {len(vis)} GO-enriched, **{len(inv)} GO-invisible** "
        f"({', '.join(inv['module'])}).",
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.go_invisible_gate "
        f"--analysis {analysis}" + (f" --region {region}" if region else "") + "`",
        "",
        "## Per-module switch coherence",
        "",
        _markdown_table(gate[["module", "go_invisible", "n_genes", "pheno_fdr",
                              "genes_with_real_switch", "max_switch_strength",
                              "n_sig_switch_tx", "top_switch_genes"]]),
        "",
        f"Pooled background functional-consequence fractions: {bg}.",
        "",
        "## Reading",
        "",
        "- `genes_with_real_switch` near `n_genes` = members carry genuine DTU "
        "(both up- and down-correlated transcripts), not abundance shifts.",
        "- Functional-consequence fractions at/above background and indistinguishable "
        "from GO-visible disease modules => GO-invisible modules are not lower quality; "
        "GO-invisibility reflects GO's gene-level/abundance bias, blind to DTU.",
        "- `top_switch_genes` heterogeneous within a module (shared switch axis, not a "
        "shared GO process) => frame as a complementary DTU-without-DGE layer, NOT "
        "pathway discovery WGCNA misses.",
    ]
    (out_dir / "GO_INVISIBLE_GATE.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="GO-invisible switch-module biology gate.")
    p.add_argument("--analysis", required=True,
                   help="e.g. brainseq-sczd, brainseq-aging, gtex-aging")
    p.add_argument("--region", default=None)
    p.add_argument("--variant", default="standard")
    p.add_argument("--fdr", type=float, default=0.10)
    p.add_argument("--top-n", type=int, default=8)
    p.add_argument("--gtf", default=str(DEFAULT_GTF_PATH))
    p.add_argument("--no-gtf", action="store_true", help="skip gene-symbol naming")
    args = p.parse_args()
    gtf = None if args.no_gtf else Path(args.gtf)
    gate = run_gate(args.analysis, args.region, args.variant, args.fdr, args.top_n, gtf)
    print(f"Wrote gate for {len(gate)} phenotype-significant modules "
          f"({int(gate['go_invisible'].sum())} GO-invisible)")


if __name__ == "__main__":
    main()
