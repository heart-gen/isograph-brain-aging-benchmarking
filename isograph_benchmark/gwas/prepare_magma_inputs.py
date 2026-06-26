from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from isograph_benchmark.paths import ensure_dir, rel


DEFAULT_CONFIG = rel("configs", "gwas_magma.yaml")
DEFAULT_OUTPUT_DIR = rel("real_data", "gwas", "_m", "tmp", "magma_pvals")


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def _sep_arg(value: str | None) -> str:
    if value in (None, "", "whitespace"):
        return r"\s+"
    if value == "tab":
        return "\t"
    return value


def _find_column(columns: list[str], spec: Any, label: str) -> str:
    candidates = [str(item) for item in _as_list(spec)]
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"Missing {label} column; tried {candidates}, found {columns}")


def _chrom_to_numeric(series: pd.Series) -> pd.Series:
    value = series.astype("string").str.replace("^chr", "", regex=True)
    return pd.to_numeric(value, errors="coerce")


def _chunk_reader(path: Path, trait_cfg: dict[str, Any], chunksize: int):
    return pd.read_csv(
        path,
        sep=_sep_arg(trait_cfg.get("sep")),
        compression="infer",
        chunksize=chunksize,
        dtype="string",
        engine="python",
    )


def _build_output_frame(
    chunk: pd.DataFrame,
    trait_cfg: dict[str, Any],
    mhc_cfg: dict[str, Any],
    seen_snps: set[str],
) -> pd.DataFrame:
    columns = list(chunk.columns)
    snp_col = _find_column(columns, trait_cfg.get("snp_col"), "SNP")
    p_col = _find_column(columns, trait_cfg.get("p_col"), "p-value")

    out = pd.DataFrame()
    raw_snp = chunk[snp_col].astype("string").str.strip()
    if trait_cfg.get("snp_rsid_extract"):
        # SNP ids carry the rsID embedded in a compound token (e.g.
        # "rs11106131:78045649:G:T" or bare "chr:pos:a1:a2" with no rsID). The
        # MAGMA reference panel is rsID-keyed, so pull the rsID; rows without one
        # become NA and are dropped by the validity filter below.
        out["SNP"] = raw_snp.str.extract(r"(rs\d+)", expand=False).astype("string")
    else:
        out["SNP"] = raw_snp
    out["P"] = pd.to_numeric(chunk[p_col], errors="coerce")

    if "fixed_n" in trait_cfg:
        out["N"] = float(trait_cfg["fixed_n"])
    elif trait_cfg.get("n_col"):
        n_col = _find_column(columns, trait_cfg.get("n_col"), "N")
        out["N"] = pd.to_numeric(chunk[n_col], errors="coerce")
    elif trait_cfg.get("n_sum_cols"):
        total = None
        for n_col in trait_cfg["n_sum_cols"]:
            col = _find_column(columns, n_col, "N component")
            values = pd.to_numeric(chunk[col], errors="coerce")
            total = values if total is None else total + values
        out["N"] = total
    else:
        raise ValueError("Trait config must specify fixed_n, n_col, or n_sum_cols")

    valid = (
        out["SNP"].notna()
        & (out["SNP"] != "")
        & (out["SNP"] != "NA")
        & out["P"].between(0, 1, inclusive="both")
        & out["N"].gt(0)
    )

    if trait_cfg.get("chr_col") and trait_cfg.get("pos_col"):
        chr_col = _find_column(columns, trait_cfg.get("chr_col"), "chromosome")
        pos_col = _find_column(columns, trait_cfg.get("pos_col"), "position")
        chrom = _chrom_to_numeric(chunk[chr_col])
        pos = pd.to_numeric(chunk[pos_col], errors="coerce")
        mhc_chrom = int(mhc_cfg.get("chromosome", 6))
        mhc_start = int(mhc_cfg.get("start", 25_000_000))
        mhc_end = int(mhc_cfg.get("end", 34_000_000))
        in_mhc = chrom.eq(mhc_chrom) & pos.between(mhc_start, mhc_end, inclusive="both")
        valid &= ~in_mhc

    out = out.loc[valid, ["SNP", "P", "N"]].copy()
    out["N"] = out["N"].round().astype("Int64")

    duplicate = out["SNP"].isin(seen_snps)
    out = out.loc[~duplicate]
    seen_snps.update(out["SNP"].tolist())
    return out


def prepare_trait(
    trait: str,
    trait_cfg: dict[str, Any],
    output_dir: Path,
    mhc_cfg: dict[str, Any],
    chunksize: int,
    force: bool = False,
) -> Path:
    source = Path(str(trait_cfg["path"]))
    if not source.exists():
        raise FileNotFoundError(f"{trait}: GWAS source not found: {source}")

    out_path = output_dir / f"{trait}_noMHC.tsv"
    if out_path.exists() and not force:
        print(f"{trait}: exists, skipping {out_path}")
        return out_path

    ensure_dir(output_dir)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    seen_snps: set[str] = set()
    rows_written = 0

    with tmp_path.open("w") as handle:
        handle.write("SNP\tP\tN\n")
        for chunk in _chunk_reader(source, trait_cfg, chunksize):
            normalized = _build_output_frame(chunk, trait_cfg, mhc_cfg, seen_snps)
            if normalized.empty:
                continue
            normalized.to_csv(handle, sep="\t", index=False, header=False)
            rows_written += len(normalized)

    tmp_path.replace(out_path)
    print(f"{trait}: wrote {rows_written:,} SNPs to {out_path}")
    return out_path


def prepare_from_config(
    config_path: Path = DEFAULT_CONFIG,
    traits: list[str] | None = None,
    force: bool = False,
) -> list[Path]:
    config = _read_yaml(config_path)
    trait_map = config.get("traits") or {}
    selected = traits or list(trait_map)
    output_dir = rel(str(config.get("output_dir", DEFAULT_OUTPUT_DIR)))
    chunksize = int(config.get("chunksize", 500_000))
    mhc_cfg = config.get("mhc") or {}

    paths = []
    for trait in selected:
        if trait not in trait_map:
            raise KeyError(f"Unknown trait {trait!r}; available: {sorted(trait_map)}")
        paths.append(prepare_trait(trait, trait_map[trait], output_dir, mhc_cfg, chunksize, force=force))
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare MAGMA SNP p-value inputs from GWAS summary statistics.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--trait", action="append", help="Trait id to prepare. May be repeated; default is all.")
    parser.add_argument("--force", action="store_true", help="Rebuild outputs that already exist.")
    args = parser.parse_args()
    prepare_from_config(args.config, traits=args.trait, force=args.force)


if __name__ == "__main__":
    main()
