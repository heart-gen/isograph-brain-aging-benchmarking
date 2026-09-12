from __future__ import annotations

from pathlib import Path

try:
    from pyhere import here as _pyhere
except ImportError:  # pragma: no cover - local bootstrap before env creation
    _pyhere = None


def root() -> Path:
    if _pyhere is not None:
        return Path(_pyhere()).resolve()
    current = Path.cwd().resolve()
    for candidate in [current, *current.parents]:
        if (candidate / ".here").exists():
            return candidate
    return current


def rel(*parts: str) -> Path:
    return root().joinpath(*parts)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


# --------------------------------------------------------------------------- #
# Layout registry
# --------------------------------------------------------------------------- #
# The physical layout of the repository is defined here and nowhere else.
# Analysis code addresses a stage by its logical bucket name, so a stage can be
# relocated by editing one tuple below.
OUTPUT_DIRS: dict[str, tuple[str, ...]] = {
    # 01 — synthetic benchmark (stage root; sub-stages 00_design .. 03_metrics)
    "synthetic": ("01_synthetic_benchmark",),
    # 02 — module discovery; the cohort x region artifact store lives beneath it
    "modules": ("02_module_discovery",),
    # 03 — module trust
    "trust.stability": ("03_module_trust", "_m", "stability"),
    "trust.replication": ("03_module_trust", "_m", "replication"),
    # 04 — module characterization
    "characterize": ("04_module_characterization", "_m"),
    # 05 — genetic anchoring
    "anchoring": ("05_genetic_anchoring", "_m"),
    "anchoring.coloc": ("05_genetic_anchoring", "_m", "coloc"),
    # per-gene sQTL-vs-eQTL coloc.abf contrast, and its signal-level (coloc.susie) upgrade
    "anchoring.coloc_modality": (
        "05_genetic_anchoring", "_m", "coloc_modality_contrast"),
    "anchoring.coloc_signal": (
        "05_genetic_anchoring", "_m", "coloc_signal_susie"),
    # BrainSEQ switch-QTL (S_g) / abundance-QTL (A_g) mapping and its junction layer
    "anchoring.brainseq_qtl": ("05_genetic_anchoring", "_m", "brainseq_switch_qtl"),
    # signal-level coloc of the BrainSEQ S_g / A_g QTLs (EA-only), in-sample QTL LD
    "anchoring.coloc_brainseq": ("05_genetic_anchoring", "_m", "coloc_brainseq"),
    "anchoring.smr": ("05_genetic_anchoring", "_m", "smr_heidi"),
    "anchoring.ldsc": ("05_genetic_anchoring", "_m", "ldsc"),
    "anchoring.gwas": ("05_genetic_anchoring", "_m", "gwas"),
    # 06 — switch mechanism
    "mechanism": ("06_switch_mechanism", "_m"),
    # 07 — RBP regulation
    "regulation": ("07_rbp_regulation", "_m"),
    # manuscript display items
    "manuscript": ("manuscript", "_m"),
    # PI-facing review reports and their frozen evidence inventory
    "reports": ("reports", "pi"),
    "reports.evidence": ("reports", "pi", "_evidence"),
    # shared gitignored scratch (GTF parse cache, id maps)
    "tmp": ("04_module_characterization", "_m", "tmp"),
}

COHORTS: tuple[str, ...] = ("brainseq", "gtex")

# Analysis name -> (cohort, fixed region). ``None`` means the region is supplied
# by the caller. This is the single definition of the analysis namespace used by
# every ``--analysis`` CLI flag.
ANALYSIS_COHORT: dict[str, tuple[str, str | None]] = {
    "brainseq-sczd": ("brainseq", "caudate_sczd"),
    "brainseq-aging": ("brainseq", None),
    "gtex-aging": ("gtex", None),
}


def stage_out(bucket: str, *parts: str) -> Path:
    """Path inside a stage output bucket, e.g. ``stage_out("mechanism", "switch_validation")``."""
    try:
        base = OUTPUT_DIRS[bucket]
    except KeyError:  # pragma: no cover - programming error
        raise KeyError(
            f"Unknown output bucket {bucket!r}; known: {sorted(OUTPUT_DIRS)}"
        ) from None
    return rel(*base, *parts)


def stage_dir(stage: str, *parts: str) -> Path:
    """Path inside a stage directory (not its ``_m``), e.g. the synthetic sub-stages."""
    return stage_out(stage, *parts)


def cohort_dir(cohort: str, *parts: str) -> Path:
    """Path inside a cohort tree of the module-discovery artifact store."""
    return stage_out("modules", cohort, *parts)


def region_store(cohort: str, region: str, *parts: str) -> Path:
    """Path inside ``<modules>/<cohort>/<region>/_m`` — the per-region artifact store."""
    return cohort_dir(cohort, region, "_m", *parts)


def resolve_analysis(analysis: str, region: str | None = None) -> tuple[str, str]:
    """Map an ``--analysis`` name (+ region) onto its ``(cohort, region)`` store key."""
    try:
        cohort, fixed = ANALYSIS_COHORT[analysis]
    except KeyError:
        raise ValueError(f"Unknown analysis: {analysis!r}") from None
    if fixed is not None:
        return cohort, fixed
    if region is None:
        raise ValueError(f"Analysis {analysis!r} requires a region")
    return cohort, region


def analysis_store(analysis: str, region: str | None = None, *parts: str) -> Path:
    """Per-region artifact store addressed by analysis name."""
    cohort, resolved = resolve_analysis(analysis, region)
    return region_store(cohort, resolved, *parts)


def region_artifact_dir(
    analysis: str, region: str | None = None, variant: str = "standard"
) -> Path:
    """The IsoGraph fit directory for an analysis.

    ``variant="with-abundance"`` addresses the abundance-channel refit, which is
    swept and resolution-selected separately from the switch-only fit.
    """
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    return analysis_store(analysis, region, subdir)
