"""Freeze the evidence base behind the PI review reports.

One row per ``ANALYSIS_MAP.md`` analysis, recording where each analysis's code,
wrapper, outputs and logs actually are, when its outputs last *changed* (not when
they were last copied), and which prose already interprets it.

Run recency deliberately does **not** come from filesystem mtimes. Everything under
the stage ``_m/`` trees carries a single checkout timestamp from the repository
reorganisation, so mtime says when the tree was materialised, not when an analysis
ran. The reliable signal is the last git commit that *changed the content* of an
output (renames excluded, so the stage reorg does not read as a re-run), with mtime
kept only as a fallback for untracked outputs.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import COHORTS, ensure_dir, root, stage_out

# --------------------------------------------------------------------------- #
# Where each stage keeps its outputs and logs
# --------------------------------------------------------------------------- #
STAGE_DIRS: dict[str, str] = {
    "01": "01_synthetic_benchmark",
    "02": "02_module_discovery",
    "03": "03_module_characterization",
    "04": "04_module_trust",
    "05": "05_genetic_anchoring",
    "06": "06_switch_mechanism",
    "07": "07_rbp_regulation",
    "08": "08_integration",
}

# ANALYSIS_MAP.md writes per-region outputs against a placeholder store. These
# expand over the cohort x region artifact stores that actually exist on disk.
STORE_PLACEHOLDERS = ("<store>", "<cohort>/<region>/_m", "<caudate_sczd store>")


def region_stores() -> list[Path]:
    """Every ``<modules>/<cohort>/<region>/_m`` directory present on disk."""
    stores: list[Path] = []
    for cohort in COHORTS:
        cdir = stage_out("modules", cohort)
        if not cdir.exists():
            continue
        for region in sorted(d for d in cdir.iterdir() if d.is_dir() and d.name != "_m"):
            store = region / "_m"
            if store.exists():
                stores.append(store)
    return stores

# Directories that never hold logs but do hold very many entries.
_WALK_SKIP = {
    "_h", ".git", "datasets", "datasets_residual_cost", "partitions", "_shards",
    "susie", "abf", "arms", "figures", "module_interpret", "beds", "ldscores",
    "sumstats", "genotypes", "results",
}

_MAP_ROW = re.compile(r"^\|(?P<cells>.+)\|\s*$")
_STAGE_HEAD = re.compile(r"^##\s+(?P<num>\d{2})\s+—\s+(?P<name>.+?)\s*$")
_CODE = re.compile(r"`([^`]+)`")


@dataclass
class Row:
    stage: str
    stage_name: str
    analysis: str
    cli_raw: str
    wrapper_raw: str
    outputs_raw: str
    display_raw: str
    cli_paths: list[str] = field(default_factory=list)
    wrapper_paths: list[str] = field(default_factory=list)
    output_paths: list[str] = field(default_factory=list)


def _cells(line: str) -> list[str] | None:
    m = _MAP_ROW.match(line.rstrip())
    if m is None:
        return None
    return [c.strip() for c in m.group("cells").split("|")]


def parse_map(map_path: Path) -> list[Row]:
    """One Row per data row of every stage table in ANALYSIS_MAP.md."""
    rows: list[Row] = []
    stage = stage_name = ""
    for line in map_path.read_text().splitlines():
        head = _STAGE_HEAD.match(line)
        if head:
            stage, stage_name = head.group("num"), head.group("name")
            continue
        cells = _cells(line)
        if cells is None or not stage or len(cells) < 5:
            continue
        if cells[0].startswith("---") or cells[0] == "Analysis":
            continue
        rows.append(
            Row(
                stage=stage,
                stage_name=stage_name,
                analysis=cells[0],
                cli_raw=cells[1],
                wrapper_raw=cells[2],
                outputs_raw=cells[3],
                display_raw=cells[4],
            )
        )
    return rows


# --------------------------------------------------------------------------- #
# Path resolution
# --------------------------------------------------------------------------- #
def _candidate_cli(token: str) -> list[Path]:
    """Expand a CLI cell token into concrete module paths under the package."""
    token = token.strip().rstrip(",")
    if not token or token.startswith("("):
        return []
    # brace expansion: real_data/coloc_{prep,summary}.py
    brace = re.match(r"^(.*?)\{([^}]*)\}(.*)$", token)
    if brace:
        pre, opts, post = brace.groups()
        return [
            p
            for opt in opts.split(",")
            for p in _candidate_cli(f"{pre}{opt.strip()}{post}")
        ]
    if not token.endswith(".py") and "/" in token:
        token = f"{token}.py"
    cands = [root() / "isograph_benchmark" / token, root() / token]
    return [p for p in cands if p.exists()] or [root() / "isograph_benchmark" / token]


def _candidate_wrapper(token: str, stage: str) -> list[Path]:
    """Expand a wrapper cell token into concrete files under ``<stage>/_h/``."""
    token = token.strip().rstrip(",")
    if not token or token.startswith("(") or token.lower() == "—":
        return []
    stage_dir = STAGE_DIRS.get(stage)
    if stage_dir is None:
        return []
    hdir = root() / stage_dir / "_h"
    if not hdir.exists():
        return []
    # A cell is either an explicit path, or a numbering hint like "_h/01–02".
    if token.endswith(".sh") or token.endswith(".R"):
        p = root() / token if "/" in token else hdir / token
        return [p]
    nums = re.findall(r"(\d{2})", token.replace("–", "-"))
    if not nums:
        return []
    lo, hi = nums[0], nums[-1]
    wanted = {f"{n:02d}" for n in range(int(lo), int(hi) + 1)}
    return sorted(f for f in hdir.iterdir() if f.name[:2] in wanted)


def _expand(token: str) -> list[str]:
    """Expand the map's path notation: ``{a,b}`` alternation and ``[x]`` optionals."""
    brace = re.match(r"^(.*?)\{([^}]*)\}(.*)$", token)
    if brace:
        pre, opts, post = brace.groups()
        return [x for opt in opts.split(",") for x in _expand(f"{pre}{opt.strip()}{post}")]
    opt = re.match(r"^(.*?)\[([^\]]*)\](.*)$", token)
    if opt:
        pre, inner, post = opt.groups()
        return _expand(f"{pre}{post}") + _expand(f"{pre}{inner}{post}")
    return [token]


# Extensions to try when the map names an output without one, and the model
# subdirectory several per-region outputs actually live under.
_EXT_FALLBACKS = (".parquet", ".json", ".csv", ".md")
_MODEL_SUBDIR = "isograph_vae"


def _resolve_with_fallback(cand: Path) -> tuple[Path, str]:
    """(resolved path, discrepancy note). Empty note means the map was accurate.

    The map is prose maintained by hand, so a declared path can be right about the
    analysis and wrong about the file. Rather than reporting such an analysis as
    never run, resolve it and record precisely how the map disagrees with disk.
    """
    if cand.exists():
        return cand, ""
    for ext in _EXT_FALLBACKS:
        if cand.with_suffix(ext).exists():
            return cand.with_suffix(ext), f"{cand.name} -> {cand.name}{ext}"
    # Several per-region outputs sit inside the fit directory, not the store root.
    nested = cand.parent / _MODEL_SUBDIR / cand.name
    if nested.exists():
        return nested, f"{cand.name} -> {_MODEL_SUBDIR}/{cand.name}"
    stem_hits = sorted(
        h
        for h in ((cand.parent / _MODEL_SUBDIR).glob(f"{cand.stem}*") if (cand.parent / _MODEL_SUBDIR).exists() else [])
    )
    if stem_hits:
        return stem_hits[0], f"{cand.name} -> {_MODEL_SUBDIR}/{stem_hits[0].name}"
    return cand, ""


def _resolve_one(tok: str, stage_dir: str | None) -> list[Path]:
    """Concrete paths for a single (already brace-expanded) output token."""
    tok = tok.strip().rstrip("/").lstrip("/")
    if not tok:
        return []
    bases = [root()]
    if stage_dir and not tok.startswith(tuple(STAGE_DIRS.values())):
        bases.insert(0, root() / stage_dir)
    if "*" in tok:
        # The wildcard may be in ANY segment -- `_m/ase_junction_switch/<region>/x.parquet`
        # becomes `.../*/x.parquet` -- so glob the whole relative pattern from each base
        # rather than only the final name. Globbing `parent.name` under a parent that itself
        # contains `*` finds nothing, which is what made the per-region allelic outputs read
        # as `not_run` until 2026-09-20.
        hits: list[Path] = []
        for base in bases:
            if base.exists():
                hits.extend(sorted(base.glob(tok)))
        return hits
    existing = [b / tok for b in bases if (b / tok).exists()]
    return existing or [bases[0] / tok]


def _candidate_outputs(row: Row) -> tuple[list[Path], int, list[str]]:
    """(paths, n_stores_expanded, map_discrepancies) for a row's declared outputs.

    Region-store placeholders fan out over every cohort x region store on disk, so a
    per-region analysis is judged on how many of its stores actually carry the output
    rather than being called unrun. Where the map names a path that does not exist but
    an equivalent one does, the resolved path is used and the disagreement recorded.
    """
    out: list[Path] = []
    notes: set[str] = set()
    n_stores = 0
    stage_dir = STAGE_DIRS.get(row.stage)
    for token in _CODE.findall(row.outputs_raw):
        for tok in _expand(token.strip()):
            tok = tok.strip()
            placeholder = next((m for m in STORE_PLACEHOLDERS if m in tok), None)
            if placeholder:
                suffix = tok.split(placeholder, 1)[1].strip("/")
                stores = region_stores()
                if placeholder == "<caudate_sczd store>":
                    stores = [st for st in stores if st.parent.name == "caudate_sczd"]
                n_stores = max(n_stores, len(stores))
                for store in stores:
                    cand = store / suffix if suffix else store
                    if "*" in cand.name:
                        out.extend(sorted(cand.parent.glob(cand.name)))
                        continue
                    resolved, note = _resolve_with_fallback(cand)
                    out.append(resolved)
                    if note:
                        notes.add(note)
                continue
            # A remaining <placeholder> is a wildcard over that path segment.
            tok = re.sub(r"<[^>]+>", "*", tok)
            for cand in _resolve_one(tok, stage_dir):
                if "*" in str(cand):
                    continue          # _resolve_one already globbed; a leftover * found nothing
                resolved, note = _resolve_with_fallback(cand)
                out.append(resolved)
                if note:
                    notes.add(note)
    seen: set[Path] = set()
    uniq = [pp for pp in out if not (pp in seen or seen.add(pp))]
    return uniq, n_stores, sorted(notes)


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #
def _git(args: list[str]) -> str:
    try:
        return subprocess.run(
            ["git", *args],
            cwd=root(),
            capture_output=True,
            text=True,
            timeout=60,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):  # pragma: no cover
        return ""


_HISTORY: dict[str, tuple[str, str, str]] | None = None
_TRACKED: set[str] | None = None


def _build_history() -> dict[str, tuple[str, str, str]]:
    """Map repo-relative path -> (date, sha, subject) of its last content change.

    Built in two batched git passes instead of one call per file. ``--diff-filter=AM``
    with rename detection on means a pure move — such as the 2026-08-29 stage
    reorganisation — never registers as a content change; the rename pairs are
    collected separately so a file's pre-move history still reaches its new path.
    """
    hist: dict[str, tuple[str, str, str]] = {}
    out = _git(
        [
            "log",
            "--diff-filter=AM",
            "-M",
            "--name-only",
            "--date=short",
            "--format=\x01%ad\t%h\t%s",
        ]
    )
    cur: tuple[str, str, str] | None = None
    for line in out.splitlines():
        if line.startswith("\x01"):
            parts = line[1:].split("\t")
            while len(parts) < 3:
                parts.append("")
            cur = (parts[0], parts[1], parts[2])
        elif line.strip() and cur is not None:
            hist.setdefault(line.strip(), cur)  # log is newest-first

    # Chain renames so a moved file inherits the history of its former path.
    renames = _git(["log", "--diff-filter=R", "-M", "--name-status", "--format="])
    pairs: list[tuple[str, str]] = []
    for line in renames.splitlines():
        cells = line.split("\t")
        if len(cells) == 3 and cells[0].startswith("R"):
            pairs.append((cells[1], cells[2]))
    for old, new_path in pairs:
        if new_path not in hist and old in hist:
            hist[new_path] = hist[old]
    return hist


def history() -> dict[str, tuple[str, str, str]]:
    global _HISTORY
    if _HISTORY is None:
        _HISTORY = _build_history()
    return _HISTORY


def tracked_files() -> set[str]:
    global _TRACKED
    if _TRACKED is None:
        _TRACKED = set(_git(["ls-files"]).splitlines())
    return _TRACKED


def content_change(path: Path) -> tuple[str, str, str]:
    """(date, sha, subject) of the last commit that changed this path's content.

    A directory takes the newest content change among the files inside it.
    """
    rp = str(path.relative_to(root()))
    hist = history()
    if rp in hist:
        return hist[rp]
    prefix = rp + "/"
    inside = [v for k, v in hist.items() if k.startswith(prefix)]
    if not inside:
        return "", "", ""
    return max(inside, key=lambda t: t[0])


def is_tracked(path: Path) -> bool:
    rp = str(path.relative_to(root()))
    if rp in tracked_files():
        return True
    prefix = rp + "/"
    return any(f.startswith(prefix) for f in tracked_files())


def newest_mtime(paths: list[Path]) -> str:
    stamps = []
    for p in paths:
        if p.is_dir():
            stamps.extend(c.stat().st_mtime for c in p.iterdir() if c.is_file())
        elif p.exists():
            stamps.append(p.stat().st_mtime)
    if not stamps:
        return ""
    return datetime.fromtimestamp(max(stamps)).strftime("%Y-%m-%d %H:%M")


def prose_docs(paths: list[Path]) -> list[Path]:
    """Markdown sitting beside an analysis's outputs — its existing interpretation."""
    docs: list[Path] = []
    for p in paths:
        d = p if p.is_dir() else p.parent
        if not d.exists():
            continue
        docs.extend(sorted(f for f in d.iterdir() if f.suffix == ".md"))
    seen: set[Path] = set()
    return [x for x in docs if not (x in seen or seen.add(x))]


def _log_dirs(stage: str) -> list[Path]:
    """Every ``logs/`` directory beneath a stage, wherever it ended up.

    Wrappers written before the 2026-08-29 reorganisation pointed sbatch at a bare
    ``_m/logs``, so a stage's logs can sit under a sibling tree (notably the cohort
    stores under ``02_module_discovery``) rather than under the stage that owns the
    analysis. Walking the stage tree finds the nested ones; the displaced ones are
    recorded as a caveat in the reports rather than silently attributed.
    """
    stage_dir = STAGE_DIRS.get(stage)
    if stage_dir is None:
        return []
    base = root() / stage_dir
    if not base.exists():
        return []
    # Every real logs/ dir sits within four levels of a stage root; the walk is
    # depth-capped and prunes the bulk data trees, which hold thousands of
    # directories and no logs.
    found: list[Path] = []
    stack = [(base, 0)]
    while stack:
        cur, depth = stack.pop()
        if depth > 4:
            continue
        try:
            children = [c for c in cur.iterdir() if c.is_dir()]
        except OSError:  # pragma: no cover - unreadable dir
            continue
        for child in children:
            if child.name == "logs":
                found.append(child)
            elif child.name not in _WALK_SKIP:
                stack.append((child, depth + 1))
    return sorted(found)


_LOG_STATS: dict[str, tuple[int, str]] = {}


def log_stats(stage: str) -> tuple[int, str]:
    """(number of SLURM logs, newest log timestamp) for a stage, computed once.

    Logs are gitignored, so unlike the committed ``_m/`` outputs their mtimes were
    written by the job itself rather than by the checkout that materialised the
    tree — this is the one filesystem timestamp in the repository that is a genuine
    run date. ``scandir`` carries the stat with the entry, which matters because the
    synthetic stage alone holds ~5,000 log files on a slow shared filesystem.
    """
    if stage in _LOG_STATS:
        return _LOG_STATS[stage]
    count = 0
    newest = 0.0
    for d in _log_dirs(stage):
        try:
            with os.scandir(d) as entries:
                for e in entries:
                    if not e.is_file() or not e.name.endswith((".log", ".out")):
                        continue
                    count += 1
                    newest = max(newest, e.stat().st_mtime)
        except OSError:  # pragma: no cover - unreadable dir
            continue
    stamp = datetime.fromtimestamp(newest).strftime("%Y-%m-%d %H:%M") if newest else ""
    _LOG_STATS[stage] = (count, stamp)
    return _LOG_STATS[stage]


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def classify(outputs: list[Path], n_stores: int) -> tuple[str, str]:
    """(run_status, note) for an analysis, from what is actually on disk.

    For a per-region analysis, partial coverage of the cohort x region stores is a
    real finding (several are deliberate spot-checks, several are gaps), so it is
    reported as coverage rather than collapsed into pass/fail.
    """
    if not outputs:
        return "unresolved", "no output path could be resolved from the map"
    present = [p for p in outputs if p.exists()]
    if not present:
        return "not_run", "no declared output exists"
    empty = [p for p in present if p.is_dir() and not any(p.iterdir())]
    if empty:
        return "partially_run", f"{len(empty)} declared output dir(s) empty"
    if len(present) < len(outputs):
        missing = len(outputs) - len(present)
        if n_stores:
            return (
                "partial_coverage",
                f"present in {len(present)} of {len(outputs)} cohort x region stores",
            )
        return "partially_run", f"{missing} of {len(outputs)} declared outputs missing"
    return "run", ""


def build() -> pd.DataFrame:
    rows = parse_map(root() / "ANALYSIS_MAP.md")
    stages = sorted({r.stage for r in rows})
    stage_log = {st: log_stats(st) for st in stages}
    records = []
    for row in rows:
        cli = [p for tok in _CODE.findall(row.cli_raw) for p in _candidate_cli(tok)]
        wrappers = [
            p for tok in _CODE.findall(row.wrapper_raw) for p in _candidate_wrapper(tok, row.stage)
        ]
        outputs, n_stores, map_notes = _candidate_outputs(row)
        status, note = classify(outputs, n_stores)
        present = [p for p in outputs if p.exists()]

        changes = [content_change(p) for p in present if is_tracked(p)]
        changes = [c for c in changes if c[0]]
        last_change = max((c[0] for c in changes), default="")
        change_commit = next((c[1] for c in changes if c[0] == last_change), "")
        change_subject = next((c[2] for c in changes if c[0] == last_change), "")
        tracked = sum(1 for p in present if is_tracked(p))

        docs = prose_docs(present)

        records.append(
            {
                "stage": row.stage,
                "stage_name": row.stage_name,
                "analysis": row.analysis,
                "run_status": status,
                "status_note": note,
                "cli_paths": ";".join(str(p.relative_to(root())) for p in cli),
                "cli_missing": ";".join(
                    str(p.relative_to(root())) for p in cli if not p.exists()
                ),
                "wrapper_paths": ";".join(str(p.relative_to(root())) for p in wrappers),
                "wrapper_missing": ";".join(
                    str(p.relative_to(root())) for p in wrappers if not p.exists()
                ),
                "n_stores_expected": n_stores,
                "map_discrepancy": "; ".join(map_notes),
                "n_outputs_declared": len(outputs),
                "n_outputs_present": len(present),
                "output_paths": ";".join(str(p.relative_to(root())) for p in present),
                "outputs_missing": ";".join(
                    str(p.relative_to(root())) for p in outputs if not p.exists()
                ),
                "n_outputs_tracked": tracked,
                "last_content_change": last_change,
                "change_commit": change_commit,
                "change_subject": change_subject,
                "newest_mtime": newest_mtime(present),
                "prose_docs": ";".join(str(p.relative_to(root())) for p in docs),
                "n_stage_logs": stage_log[row.stage][0],
                "stage_newest_log": stage_log[row.stage][1],
                "display_items": row.display_raw,
            }
        )
    return pd.DataFrame.from_records(records)


def to_markdown(df: pd.DataFrame) -> str:
    """A readable companion to the parquet — the status board, per stage."""
    counts = df["run_status"].value_counts().to_dict()
    lines = [
        "# Evidence inventory",
        "",
        "One row per `ANALYSIS_MAP.md` analysis. Generated by "
        "`python -m isograph_benchmark.reporting.evidence_inventory`; do not edit by hand.",
        "",
        "**Run recency is the last git commit that changed an output's content** "
        "(`--diff-filter=AM` with rename detection, so the 2026-08-29 stage "
        "reorganisation does not read as a re-run). Filesystem mtimes are *not* run "
        "dates — every `_m/` tree carries a single checkout timestamp. The one "
        "trustworthy filesystem timestamp is the SLURM logs, which are gitignored and "
        "so were written by the jobs themselves; those are the `newest log` column.",
        "",
        f"**Status board.** {len(df)} analyses: "
        + ", ".join(f"{v} {k}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1]))
        + ".",
        "",
    ]
    for stage, grp in df.groupby("stage"):
        lines += [
            f"## {stage} — {grp['stage_name'].iloc[0]}",
            "",
            "| Analysis | Status | Last content change | Outputs (present/declared) | Tracked | Display |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for _, r in grp.iterrows():
            change = r["last_content_change"] or "—"
            if r["change_commit"]:
                change = f"{change} `{r['change_commit']}`"
            lines.append(
                f"| {r['analysis']} | {r['run_status']} | {change} | "
                f"{r['n_outputs_present']}/{r['n_outputs_declared']} | "
                f"{r['n_outputs_tracked']} | {r['display_items']} |"
            )
        lines.append("")
        flagged = grp[grp["run_status"] != "run"]
        for _, r in flagged.iterrows():
            note = r["status_note"] or "—"
            lines.append(f"- **{r['analysis']}** — {r['run_status']}: {note}.")
            if r["outputs_missing"]:
                lines.append(f"  Missing: `{r['outputs_missing'].replace(';', '`, `')}`")
        if len(flagged):
            lines.append("")

    disc = df[df["map_discrepancy"] != ""]
    if len(disc):
        lines += [
            "## ANALYSIS_MAP.md vs disk",
            "",
            "Where the map names an output path that does not exist but an equivalent one "
            "does. The analysis ran; the index is wrong about where its output lands.",
            "",
            "| Stage | Analysis | Declared -> actual |",
            "| --- | --- | --- |",
        ]
        for _, r in disc.iterrows():
            lines.append(f"| {r['stage']} | {r['analysis']} | `{r['map_discrepancy']}` |")
        lines.append("")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check",
        action="store_true",
        help="rebuild and compare against the committed inventory instead of writing",
    )
    args = ap.parse_args(argv)

    df = build().sort_values(["stage", "analysis"]).reset_index(drop=True)
    out_dir = stage_out("reports.evidence")
    pq = out_dir / "inventory.parquet"
    md = out_dir / "inventory.md"

    if args.check:
        if not pq.exists():
            print(f"FAIL: {pq} does not exist; run without --check first")
            return 1
        prev = pd.read_parquet(pq)
        same = prev.equals(df)
        print("inventory stable" if same else "inventory CHANGED since last write")
        if not same:
            print(f"  rows before/after: {len(prev)}/{len(df)}")
        return 0 if same else 1

    ensure_dir(out_dir)
    df.to_parquet(pq, index=False)
    md.write_text(to_markdown(df))
    summary = df["run_status"].value_counts().to_dict()
    print(f"wrote {pq.relative_to(root())} ({len(df)} analyses)")
    print(f"wrote {md.relative_to(root())}")
    print("status:", json.dumps(summary))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
