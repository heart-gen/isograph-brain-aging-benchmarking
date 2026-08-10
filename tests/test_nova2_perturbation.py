from __future__ import annotations

import gzip

import pandas as pd

from isograph_benchmark.real_data.nova2_perturbation import (
    _event_flanks,
    _localization_effect,
    _overlap_event_flanks,
    _read_bed12,
    _read_event_statistics,
)


def _write_statistics(path, rows: list[list[object]]) -> None:
    header = [
        "Ensemble_gene//Gene_symbol",
        "chr",
        "chr.Start",
        "chr.End",
        "AS.name",
        "strand",
        "AS.type",
        "isoformIDs",
        "coverage",
        "Exon.inclusion_rate.group1",
        "Exon.inclusion_rate.group2",
        "delta.exon.inclusion.rate(dI)",
        "pvalue",
        "FDR",
    ]
    with gzip.open(path, "wt", newline="") as handle:
        handle.write("\t".join(header) + "\r\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\r\n")


def test_event_statistics_handle_crlf_and_duplicate_ids(tmp_path) -> None:
    path = tmp_path / "events.txt.gz"
    common = [
        "GENE//Gene",
        "chr1",
        1,
        100,
        "event",
        "+",
        "cass",
        "INC/SKIP",
        20,
        0.2,
        0.4,
    ]
    _write_statistics(
        path,
        [
            [*common, 0.20, 0.01, 0.04],
            [*common, -0.30, 0.001, 0.02],
            [
                "G2//G2",
                "chr2",
                1,
                50,
                "below",
                "-",
                "cass",
                "INC/SKIP",
                20,
                0.2,
                0.25,
                0.05,
                0.5,
                0.5,
            ],
        ],
    )
    frame = _read_event_statistics(path, "emx1", 0.05, 0.10)
    assert list(frame["event_id"]) == ["below", "event"]
    selected = frame.set_index("event_id").loc["event"]
    assert selected["delta_inclusion"] == -0.30
    assert selected["source_stat_rows"] == 2
    assert bool(selected["passes_event_thresholds"])
    assert not bool(frame.set_index("event_id").loc["below", "passes_event_thresholds"])


def test_bed12_flanks_preserve_strand_semantics_and_short_introns(tmp_path) -> None:
    stat_path = tmp_path / "events.txt.gz"
    rows = []
    for gene, chrom, event, strand in [
        ("G1//G1", "chr1", "plus", "+"),
        ("G2//G2", "chr2", "minus", "-"),
        ("G3//G3", "chr3", "short", "+"),
    ]:
        rows.append(
            [
                gene,
                chrom,
                1,
                500,
                event,
                strand,
                "cass",
                "INC/SKIP",
                30,
                0.2,
                0.5,
                0.3,
                0.001,
                0.01,
            ]
        )
    _write_statistics(stat_path, rows)
    bed_path = tmp_path / "events.bed.gz"
    with gzip.open(bed_path, "wt") as handle:
        handle.write("chr1\t100\t320\tplus\t0\t+\t100\t320\t0\t2\t10,10\t0,210\n")
        handle.write("chr2\t100\t320\tminus\t0\t-\t100\t320\t0\t2\t10,10\t0,210\n")
        handle.write("chr3\t100\t160\tshort\t0\t+\t100\t160\t0\t2\t10,10\t0,50\n")
    statistics = _read_event_statistics(stat_path, "ctx", 0.05, 0.10)
    bed12 = _read_bed12(bed_path, "ctx")
    flanks = _event_flanks(statistics, bed12, [50])
    plus = flanks[flanks["event_id"].eq("plus")].set_index("structural_class")
    minus = flanks[flanks["event_id"].eq("minus")].set_index("structural_class")
    short = flanks[flanks["event_id"].eq("short")]
    assert (plus.loc["donor_flank", ["start", "end"]].tolist()) == [110, 160]
    assert (plus.loc["acceptor_flank", ["start", "end"]].tolist()) == [260, 310]
    assert (minus.loc["acceptor_flank", ["start", "end"]].tolist()) == [110, 160]
    assert (minus.loc["donor_flank", ["start", "end"]].tolist()) == [260, 310]
    assert short.iloc[0]["structural_class"] == "short_intron"
    assert short.iloc[0][["start", "end"]].tolist() == [110, 150]


def test_event_overlap_requires_width_strand_and_structural_class() -> None:
    human = pd.DataFrame(
        {
            "window_id": ["hit", "wrong_strand", "wrong_class", "wrong_width"],
            "candidate_id": ["c"] * 4,
            "gene_id": ["g"] * 4,
            "match_id": ["m1", "m2", "m3", "m4"],
            "window_role": ["case"] * 4,
            "match_status": ["matched"] * 4,
            "chrom": ["chr1"] * 4,
            "start": [100] * 4,
            "end": [200] * 4,
            "strand": ["+", "-", "+", "+"],
            "window_width": [100, 100, 100, 50],
            "structural_class": [
                "donor_flank",
                "donor_flank",
                "acceptor_flank",
                "donor_flank",
            ],
        }
    )
    event = pd.DataFrame(
        {
            "window_id": ["event_flank"],
            "event_id": ["event"],
            "event_coordinate_id": ["coord"],
            "window_width": [100],
            "structural_class": ["donor_flank"],
            "fdr": [0.01],
            "delta_inclusion": [0.2],
            "absolute_delta_inclusion": [0.2],
            "reciprocal_mapped": [True],
            "forward_mapped_chrom": ["chr1"],
            "forward_mapped_start": [150],
            "forward_mapped_end": [250],
            "forward_mapped_strand": ["+"],
        }
    )
    calls, overlaps = _overlap_event_flanks(human, event, "ctx")
    localized = calls.set_index("window_id")["localized"]
    assert bool(localized["hit"])
    assert not bool(localized["wrong_strand"])
    assert not bool(localized["wrong_class"])
    assert not bool(localized["wrong_width"])
    assert overlaps.iloc[0]["overlap_nt"] == 50


def test_localization_effect_uses_exact_paired_discordance() -> None:
    frame = pd.DataFrame(
        {
            "case_callable": [True] * 12,
            "control_callable": [True] * 12,
            "case_localized": [True] * 10 + [False] * 2,
            "control_localized": [False] * 10 + [True] * 2,
        }
    )
    effect = _localization_effect(frame, "gene")
    assert effect["case_only_discordant"] == 10
    assert effect["control_only_discordant"] == 2
    assert effect["informative_discordant"] == 12
    assert effect["matched_odds_ratio_haldane"] == 10.5 / 2.5
