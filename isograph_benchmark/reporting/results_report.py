"""Assemble the collaborator-facing results report built from the PI review.

The prose in ``reports/results/`` is hand-written. What this module owns is the
mechanical half, so that it is reproducible rather than a sequence of ad-hoc copies:

``figures``
    stage the display items the report cites into ``reports/results/figures/``, so the
    report tree is self-contained and survives conversion to Word or PDF.
``build``
    compile ``results_briefing.tex`` with pdflatex (twice, for frame numbers).
``sync``
    copy the whole tree into the manuscript repository, which is where the shareable
    copy lives. The destination is an argument, never a baked-in path.

The figure manifest below is the single definition of which display items the report is
allowed to show. A figure named here must exist; a missing one is an error rather than a
silently absent image, because a report that cites a figure it does not carry is worse
than one that fails to build.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from isograph_benchmark.paths import ensure_dir, root, stage_out

# --------------------------------------------------------------------------- #
# Figure manifest
# --------------------------------------------------------------------------- #
# staged name -> source path, repository-relative. Staged names keep the display-item
# identity used in FIGURE_ORDERING.md so a reader can trace a figure back to its claim.
FIGURE_MANIFEST: dict[str, tuple[str, ...]] = {
    # synthetic benchmark (question 1)
    "fig1_benchmark_overview.png": (
        "01_synthetic_benchmark",
        "03_metrics",
        "figures",
        "fig1_benchmark_overview.png",
    ),
    # real-data display items (questions 2-6)
    **{
        f"{name}.png": ("manuscript", "_m", "figures", f"{name}.png")
        for name in (
            "figConceptOverview",
            "figSeparation",
            "figGoInvisible",
            "figIsaConcordance",
            "figTrustFunnel",
            "figCrossCohortReplication",
            "figBaselineRates",
            "figOrthogonalConfirm",
            "figSwitchConsequence",
            "figCompositionRobustness",
            "figQtlSpecificity",
            "figGeneticAnchoring",
            "figRbpRegulon",
        )
    },
}

# Files that make up the report tree, in the order they are listed to the user.
REPORT_FILES: tuple[str, ...] = (
    "results_report.md",
    "appendix_genetics.md",
    "appendix_mechanism.md",
    "results_briefing.tex",
    "results_briefing.pdf",
)

DEFAULT_DEST_ENV = "ISOGRAPH_MANUSCRIPT_REPO"
DEFAULT_DEST_RELATIVE = ("..", "..", "manuscript", "isograph-brain-manuscript")
DEST_SUBDIR = ("drafts", "results-report")


def report_dir() -> Path:
    return stage_out("reports.results")


def figures_dir() -> Path:
    return stage_out("reports.results.figures")


def default_dest() -> Path:
    """Manuscript repository root: the environment override, else the sibling checkout."""
    env = os.environ.get(DEFAULT_DEST_ENV)
    if env:
        return Path(env).expanduser()
    return root().joinpath(*DEFAULT_DEST_RELATIVE)


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def stage_figures() -> list[Path]:
    dest = ensure_dir(figures_dir())
    missing = [
        name for name, parts in FIGURE_MANIFEST.items() if not root().joinpath(*parts).exists()
    ]
    if missing:
        raise FileNotFoundError(
            "figure manifest names files that do not exist: " + ", ".join(sorted(missing))
        )
    staged = []
    for name, parts in FIGURE_MANIFEST.items():
        target = dest / name
        shutil.copy2(root().joinpath(*parts), target)
        staged.append(target)
    return staged


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def build_briefing() -> Path:
    tex = report_dir() / "results_briefing.tex"
    if not tex.exists():
        raise FileNotFoundError(f"{tex} does not exist; write the briefing first")
    for _ in range(2):  # second pass resolves frame numbers
        proc = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
            cwd=tex.parent,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            tail = "\n".join(proc.stdout.strip().splitlines()[-25:])
            raise RuntimeError(f"pdflatex failed:\n{tail}")
    for suffix in (".aux", ".log", ".nav", ".out", ".snm", ".toc"):
        tex.with_suffix(suffix).unlink(missing_ok=True)
    return tex.with_suffix(".pdf")


# --------------------------------------------------------------------------- #
# sync
# --------------------------------------------------------------------------- #
def sync(dest_repo: Path) -> Path:
    src = report_dir()
    if not src.exists():
        raise FileNotFoundError(f"{src} does not exist; nothing to sync")
    if not dest_repo.exists():
        raise FileNotFoundError(f"manuscript repository not found at {dest_repo}")
    dest = dest_repo.joinpath(*DEST_SUBDIR)
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(src, dest, ignore=shutil.ignore_patterns("*.aux", "*.log", "*.nav", "*.out", "*.snm", "*.toc"))
    return dest


# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)

    sub.add_parser("figures", help="stage the cited display items into reports/results/figures/")
    sub.add_parser("build", help="compile results_briefing.tex with pdflatex")

    p_sync = sub.add_parser("sync", help="copy the report tree into the manuscript repository")
    p_sync.add_argument(
        "--dest",
        type=Path,
        default=None,
        help=(
            "manuscript repository root; defaults to $"
            f"{DEFAULT_DEST_ENV}, else the sibling checkout"
        ),
    )

    args = ap.parse_args(argv)

    if args.command == "figures":
        staged = stage_figures()
        print(f"staged {len(staged)} figures into {figures_dir().relative_to(root())}")
        return 0

    if args.command == "build":
        pdf = build_briefing()
        print(f"built {pdf.relative_to(root())}")
        return 0

    dest = sync(args.dest or default_dest())
    present = [name for name in REPORT_FILES if (dest / name).exists()]
    print(f"synced {report_dir().relative_to(root())} -> {dest}")
    print("  " + ", ".join(present))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
