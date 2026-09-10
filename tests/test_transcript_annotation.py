"""The BrainSEQ bundles feed every IsoGraph fit, so their input must be regenerable.

`build_bundles` reads `inputs/raw/brainseq/annotations/transcript-annotation.tsv` to
attach `gene_id`, `transcript_name` and `transcript_type` to the RSEM transcript matrix.
That file is derived, gitignored and was absent -- so a clean checkout could not rebuild
the bundles at all, and nothing noticed because the failure was a bare FileNotFoundError
deep inside bundle construction.

These tests pin the two things that make the regeneration trustworthy: that it reproduces
the annotation the committed bundles were actually built with, and that it carries real
TRANSCRIPT names, because the DRD2 sanity check degrades silently rather than failing if
it is handed gene symbols instead.
"""
from __future__ import annotations

import pandas as pd
import pytest

from isograph_benchmark.inputs.build_transcript_annotation import COLUMNS
from isograph_benchmark.paths import rel

_ANNOT = rel("inputs", "raw", "brainseq", "annotations", "transcript-annotation.tsv")
_BUNDLE = rel("inputs", "bundles", "brainseq_v1", "caudate", "transcripts.parquet")

pytestmark = pytest.mark.skipif(
    not _ANNOT.exists() or not _BUNDLE.exists(),
    reason="regenerate with `python -m isograph_benchmark.inputs.build_transcript_annotation`")


@pytest.fixture(scope="module")
def annot() -> pd.DataFrame:
    return pd.read_csv(_ANNOT, sep="\t")


@pytest.fixture(scope="module")
def bundle() -> pd.DataFrame:
    return pd.read_parquet(_BUNDLE)


def test_carries_the_columns_build_bundles_selects(annot):
    assert {"transcript_id", "gene_id", "transcript_name", "transcript_type"} <= set(annot.columns)
    assert list(annot.columns) == COLUMNS


def test_reproduces_the_committed_bundle_exactly(annot, bundle):
    """The bundle was written by the ORIGINAL file, so agreement proves the regeneration
    recovers that file's content rather than merely something plausible from the same
    GENCODE release."""
    m = bundle.merge(annot, on="transcript_id", how="left", suffixes=("_b", "_n"))
    assert m["transcript_name_n"].isna().sum() == 0, "bundle transcripts missing from annotation"
    for c in ("gene_id", "transcript_name", "transcript_type"):
        diff = (m[f"{c}_b"].astype(str) != m[f"{c}_n"].astype(str)).sum()
        assert diff == 0, f"{diff} rows differ on {c}"


def test_transcript_name_is_a_transcript_name_not_a_gene_symbol(annot):
    """`run_models._drd2_gene_id` matches `transcript_name.str.startswith('DRD2-')`. Handed
    gene symbols it would match nothing and print a warning instead of failing, so the
    DRD2 check would quietly stop checking. Guard the distinction explicitly."""
    drd2 = annot[annot["transcript_name"].str.startswith("DRD2-", na=False)]
    assert len(drd2) > 0
    assert drd2["gene_id"].str.startswith("ENSG00000149295").all()
    # The failure mode being guarded is substituting the shared r-variables annotation,
    # whose `gene_name` column is the GENE symbol. In a real transcript annotation the two
    # are never equal: a transcript is either `SYMBOL-NNN` or, for novel/TAGENE models
    # with no symbol-based name, a bare ENST accession (about 29% of GENCODE v47). Neither
    # is ever the bare gene symbol, so equality is the exact discriminator -- and it does
    # not depend on guessing what fraction carries the -NNN form.
    same_as_gene = (annot["transcript_name"].astype(str)
                    == annot["gene_name"].astype(str))
    assert same_as_gene.sum() == 0, (
        f"{int(same_as_gene.sum()):,} transcript_name values equal the gene symbol; "
        "this looks like a gene-level annotation, which would silently disable the "
        "DRD2 check")


def test_the_drd2_check_can_still_resolve_its_gene(annot):
    from isograph_benchmark.real_data.run_models import _drd2_gene_id
    assert _drd2_gene_id(annot) is not None


def test_transcript_ids_are_unique(annot):
    assert not annot["transcript_id"].duplicated().any()


def test_build_bundles_loader_raises_a_useful_error_when_absent(tmp_path, monkeypatch):
    """The original failure mode was a bare FileNotFoundError from inside bundle
    construction. It must now say what to run."""
    import isograph_benchmark.inputs.build_bundles as bb
    monkeypatch.setattr(bb, "rel", lambda *p: tmp_path / "nope.tsv")
    with pytest.raises(SystemExit) as e:
        bb._load_transcript_annotation()
    assert "build_transcript_annotation" in str(e.value)
