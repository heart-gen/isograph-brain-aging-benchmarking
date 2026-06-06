"""
GO enrichment density for data-driven Leiden parameter selection.

Downloads and caches three files:
  - go-basic.obo                Gene Ontology hierarchy  (OBO library)
  - goa_human.gaf.gz            UniProt → GO terms for human  (EBI, ~10 MB)
  - Homo_sapiens.gene_info.gz   ENSEMBL ↔ gene symbol lookup  (NCBI, ~5 MB)

ID mapping chain: ENSEMBL ID → gene symbol (gene_info) → GO terms (GAF).
ENSEMBL IDs are used directly as population/study keys in GOEnrichmentStudy,
so no Entrez integer conversion is needed.

Optional dependency: goatools (≥1.4).  If not installed the helper degrades
gracefully; callers receive NaN density and fall back to the trait-association
criterion.  Install via:
    pip install 'isograph[go-enrichment]'
"""
from __future__ import annotations

import gzip
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from goatools.go_enrichment import GOEnrichmentStudy
    from goatools.obo_parser import GODag

    HAS_GOATOOLS = True
except ImportError:
    HAS_GOATOOLS = False

GO_FDR_ALPHA = 0.05
# Modules with fewer mapped genes than this are skipped.
_MIN_STUDY_GENES = 5

_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
_GAF_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
_GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
    "Homo_sapiens.gene_info.gz"
)


def _download(url: str, dest: Path) -> None:
    print(f"  [GO] Downloading {dest.name} ...", flush=True)
    urllib.request.urlretrieve(url, dest)
    print(f"  [GO] {dest.name}: {dest.stat().st_size / 1e6:.1f} MB", flush=True)


def _parse_gaf_symbol2gos(gaf_path: Path) -> dict[str, set[str]]:
    """Parse goa_human.gaf.gz → {gene_symbol: set[GO_IDs]} (BP aspect only).

    GAF 2.x column layout (1-indexed in spec, 0-indexed here):
      1  DB
      2  DB_Object_ID      (UniProt accession)
      3  DB_Object_Symbol  (gene symbol, e.g. DRD2)
      5  GO_ID
      9  Aspect            P=biological_process, F=molecular_function, C=cell_component
    """
    symbol2gos: dict[str, set[str]] = {}
    with gzip.open(gaf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            symbol = cols[2]
            go_id  = cols[4]
            aspect = cols[8]
            if aspect == "P" and go_id:
                symbol2gos.setdefault(symbol, set()).add(go_id)
    return symbol2gos


def _parse_ensembl_symbol(gene_info_path: Path) -> dict[str, str]:
    """Parse Homo_sapiens.gene_info.gz → {ENSEMBL_ID: gene_symbol}.

    Column layout (tab-separated, 0-indexed):
      1  GeneID   (Entrez — not used here)
      2  Symbol   (HGNC gene symbol)
      5  dbXrefs  (pipe-separated, e.g. "Ensembl:ENSG00000149295|HGNC:3020")
    """
    mapping: dict[str, str] = {}
    with gzip.open(gene_info_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            symbol = cols[2]
            for xref in cols[5].split("|"):
                if xref.startswith("Ensembl:"):
                    mapping[xref[8:]] = symbol
    return mapping


class GoAnnotations:
    """Loads GO annotations once; serves enrichment queries for many module tables.

    Typical usage::

        ga = GoAnnotations(cache_dir=Path("inputs/go_annotations"))
        ga.prepare(background_gene_ids=all_gene_ids)   # once per dataset
        density, n_enriched = ga.enrichment_density(module_table)

    The background population is fixed at the dataset level so the enrichment
    statistic is comparable across resolutions in the same sweep.
    """

    def __init__(self, cache_dir: Path) -> None:
        if not HAS_GOATOOLS:
            raise ImportError(
                "goatools is required for GO enrichment.\n"
                "Install via:  pip install goatools\n"
                "           or pip install 'isograph[go-enrichment]'"
            )
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        obo_path       = self._cache_dir / "go-basic.obo"
        gaf_path       = self._cache_dir / "goa_human.gaf.gz"
        gene_info_path = self._cache_dir / "Homo_sapiens.gene_info.gz"

        for url, path in [
            (_OBO_URL,       obo_path),
            (_GAF_URL,       gaf_path),
            (_GENE_INFO_URL, gene_info_path),
        ]:
            if not path.exists():
                _download(url, path)

        print("[GO] Parsing GO DAG ...", flush=True)
        self._godag = GODag(str(obo_path), load_obsolete=False, prt=None)

        print("[GO] Parsing EBI GAF (human BP terms) ...", flush=True)
        symbol2gos = _parse_gaf_symbol2gos(gaf_path)
        print(f"[GO] {len(symbol2gos):,} gene symbols with BP annotations", flush=True)

        print("[GO] Building ENSEMBL → symbol mapping ...", flush=True)
        ensembl_to_symbol = _parse_ensembl_symbol(gene_info_path)

        # Compose ENSEMBL → set[GO_IDs] via the symbol lookup.
        self._ensembl_to_gos: dict[str, set[str]] = {
            eid: symbol2gos[sym]
            for eid, sym in ensembl_to_symbol.items()
            if sym in symbol2gos
        }
        print(
            f"[GO] {len(self._ensembl_to_gos):,} ENSEMBL IDs with GO annotations",
            flush=True,
        )

        self._study_obj: GOEnrichmentStudy | None = None
        self._bg_set: frozenset[str] = frozenset()

    def prepare(self, background_gene_ids: list[str]) -> None:
        """Build GOEnrichmentStudy for the given background gene population.

        Uses bare ENSEMBL IDs (version suffix stripped) as population/study
        keys — no integer Entrez conversion.  Idempotent for the same background.
        """
        # Strip ENSEMBL version suffix (e.g. ENSG00000000003.16 → ENSG00000000003)
        # so that versioned gene IDs from feature_scores match the unversioned
        # IDs stored in _ensembl_to_gos (sourced from NCBI gene_info).
        bg_assoc = {
            g.split(".")[0]: self._ensembl_to_gos[g.split(".")[0]]
            for g in background_gene_ids
            if g.split(".")[0] in self._ensembl_to_gos
        }
        bg_set = frozenset(bg_assoc)
        if bg_set == self._bg_set and self._study_obj is not None:
            return
        self._bg_set = bg_set
        n_covered = len(bg_set)
        n_total   = len(background_gene_ids)
        print(
            f"[GO] Background: {n_covered:,}/{n_total:,} genes have annotations "
            f"({n_covered / n_total:.1%})",
            flush=True,
        )
        # log=None suppresses goatools' verbose startup messages.
        self._study_obj = GOEnrichmentStudy(
            list(bg_set),
            bg_assoc,
            self._godag,
            propagate_counts=True,
            alpha=GO_FDR_ALPHA,
            methods=["fdr_bh"],
            log=None,
        )

    def enrichment_density(
        self, module_table: pd.DataFrame
    ) -> tuple[float, int]:
        """Return (density, n_enriched) for the given module partition.

        density    = n_enriched / n_tested  (NaN if no module passes size filter)
        n_enriched = modules with ≥1 GO:BP term at BH FDR ≤ GO_FDR_ALPHA
        """
        if self._study_obj is None:
            raise RuntimeError("Call prepare() before enrichment_density().")

        module_ids = module_table["module_id"].unique()
        n_enriched = 0
        n_tested   = 0
        for mod_id in module_ids:
            genes = module_table.loc[
                module_table["module_id"] == mod_id, "gene_id"
            ].tolist()
            # Strip version suffix to match bare ENSEMBL IDs in _bg_set.
            study = {g.split(".")[0] for g in genes if g.split(".")[0] in self._bg_set}
            if len(study) < _MIN_STUDY_GENES:
                continue
            n_tested += 1
            results = self._study_obj.run_study(study, prt=None)
            if any(r.p_fdr_bh < GO_FDR_ALPHA for r in results):
                n_enriched += 1

        density = n_enriched / n_tested if n_tested > 0 else np.nan
        return float(density), int(n_enriched)

    def enrich_gene_set(
        self, gene_ids: list[str], max_terms: int | None = None
    ) -> pd.DataFrame:
        """GO:BP enrichment for a single gene set against the prepared background.

        Returns a DataFrame (term_id, term_name, p_value, p_fdr_bh,
        study_count, study_n) of enriched BP terms (FDR ≤ GO_FDR_ALPHA), sorted
        by FDR.  Empty DataFrame if the set is too small or nothing is enriched.
        """
        if self._study_obj is None:
            raise RuntimeError("Call prepare() before enrich_gene_set().")
        cols = ["term_id", "term_name", "p_value", "p_fdr_bh", "study_count", "study_n"]
        study = {g.split(".")[0] for g in gene_ids if g.split(".")[0] in self._bg_set}
        if len(study) < _MIN_STUDY_GENES:
            return pd.DataFrame(columns=cols)
        results = self._study_obj.run_study(study, prt=None)
        rows = [
            {
                "term_id": r.GO,
                "term_name": r.name,
                "p_value": float(r.p_uncorrected),
                "p_fdr_bh": float(r.p_fdr_bh),
                "study_count": int(r.ratio_in_study[0]),
                "study_n": int(r.ratio_in_study[1]),
            }
            for r in results
            if r.enrichment == "e" and r.NS == "BP" and r.p_fdr_bh <= GO_FDR_ALPHA
        ]
        out = pd.DataFrame(rows, columns=cols).sort_values("p_fdr_bh").reset_index(drop=True)
        return out.head(max_terms) if max_terms else out
