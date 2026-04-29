from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, rel


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_file(src: Path, dst: Path) -> dict[str, object]:
    ensure_dir(dst.parent)
    if not dst.exists() or src.stat().st_size != dst.stat().st_size:
        shutil.copy2(src, dst)
    return {"source": str(src), "target": str(dst), "bytes": dst.stat().st_size, "sha256": sha256(dst)}


def copy_brainseq(cfg: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    raw = rel("inputs", "raw", "brainseq")
    source_root = Path(cfg["brainseq"]["source_root"])
    for region in cfg["brainseq"]["regions"]:
        if region in set(cfg["brainseq"].get("exclude_regions", [])):
            continue
        for subdir, names in {
            "counts": ["tx-counts.tsv", "gene-counts.tsv", "psi-events.tsv.gz"],
            "tpm": ["tx-tpm.tsv", "gene-tpm.csv"],
        }.items():
            for name in names:
                rows.append(copy_file(source_root / subdir / region / name, raw / subdir / region / name))
    rows.append(copy_file(Path(cfg["brainseq"]["metadata"]["rnaseq"]), raw / "metadata" / "libd_rnaseq_metadata.tab"))
    for region, files in cfg["brainseq"]["metadata"]["metrics"].items():
        for file_name in files:
            src = Path(file_name)
            rows.append(copy_file(src, raw / "metadata" / region / src.name))
    for name in ["gene-annotation.tsv", "transcript-annotation.tsv"]:
        rows.append(copy_file(source_root / "annotations" / name, raw / "annotations" / name))
    return rows


def copy_gtex(cfg: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    raw = rel("inputs", "raw", "gtex_v11")
    counts_root = Path(cfg["gtex_v11"]["counts_root"])
    metadata_root = Path(cfg["gtex_v11"]["metadata_root"])
    for name in cfg["gtex_v11"]["files"].values():
        base = counts_root if "Analysis_2025" in name else metadata_root
        rows.append(copy_file(base / name, raw / ("counts" if base == counts_root else "metadata") / name))
    return rows


def main() -> None:
    cfg = load_yaml("configs/data_sources.yaml")
    rows = copy_brainseq(cfg) + copy_gtex(cfg)
    manifest = pd.DataFrame(rows)
    out = rel("reports", "raw_copy_manifest.parquet")
    ensure_dir(out.parent)
    try:
        manifest.to_parquet(out, index=False, compression="zstd")
    except ImportError:
        out = out.with_suffix(".tsv")
        manifest.to_csv(out, sep="\t", index=False)
    print(out)


if __name__ == "__main__":
    main()
