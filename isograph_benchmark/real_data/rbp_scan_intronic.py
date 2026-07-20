"""Stage 1 (intronic scope): scan switch-isoform intronic splice-site flanks for RBP motifs.

The mature-transcript scan (`rbp_scan.py`) only sees exonic + UTR sequence, so it is blind to
the splicing-regulatory RBPs that act from *intronic* splice-site flanks (the canonical binding
niche for the spliceosome-adjacent regulators: U2AF2, PTBP, RBFOX, NOVA, hnRNP families, ...).
This script closes that gap. For every transcript that appears in any region's IsoGraph switch
pairs, it derives the introns from the GENCODE v47 GTF, extracts a strand-aware window reaching
into the intron at each splice site from the genome FASTA (pre-mRNA sense; minus-strand genes are
reverse-complemented), and scans those flanks against the same ATtRACT human PWMs with identical
MOODS settings as the mature scan. Per-(transcript, RBP) hit counts are aggregated over all of a
transcript's intronic flanks so the output schema is identical to `rbp_counts.parquet`, letting
Stage 2 (`rbp_regulon.py --scope intronic`) consume it unchanged.

Runs in the dedicated `motif` env (MOODS + pyfaidx + pandas).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyfaidx
import MOODS.scan
import MOODS.tools

from isograph_benchmark.real_data.rbp_scan import (
    _load_pwms, _needed_transcripts, _P_THRESH,
)

_REPO = Path(__file__).resolve().parents[2]
_GENOME_FA = Path(
    "/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/"
    "fasta/GRCh38.primary_assembly.genome.fa")
_GTF = Path(
    "/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/"
    "gtf/gencode.v47.annotation.gtf")
_OUT = _REPO / "real_data" / "_m" / "rbp" / "rbp_counts_intronic.parquet"
_FLANK = 100              # nt reaching into the intron from each splice site

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _transcript_id(attr: str) -> str | None:
    """Pull the versioned transcript_id out of a GTF attribute column."""
    i = attr.find('transcript_id "')
    if i < 0:
        return None
    i += len('transcript_id "')
    j = attr.find('"', i)
    return attr[i:j] if j > i else None


def _transcript_exons(needed: set[str]) -> dict[str, dict]:
    """{transcript_id: {chrom, strand, exons:[(start,end),...]}} for needed transcripts.

    Coordinates are GTF 1-based inclusive; exons are collected unsorted and sorted later.
    """
    tx: dict[str, dict] = {}
    with open(_GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = _transcript_id(f[8])
            if tid is None or tid not in needed:
                continue
            rec = tx.get(tid)
            if rec is None:
                rec = {"chrom": f[0], "strand": f[6], "exons": []}
                tx[tid] = rec
            rec["exons"].append((int(f[3]), int(f[4])))
    return tx


def _intron_flanks(rec: dict) -> list[tuple[str, int, int]]:
    """Genomic (chrom, start, end) 1-based-inclusive windows reaching into each intron.

    For each intron we take up to `_FLANK` nt at the 5' end and `_FLANK` nt at the 3' end;
    short introns (<= 2*_FLANK) are emitted once as the whole intron to avoid double-counting.
    Orientation is handled downstream (minus strand reverse-complemented), so donor/acceptor
    identity is irrelevant here.
    """
    exons = sorted(rec["exons"])
    chrom = rec["chrom"]
    out: list[tuple[str, int, int]] = []
    for (_, e1_end), (e2_start, _) in zip(exons[:-1], exons[1:]):
        istart, iend = e1_end + 1, e2_start - 1
        L = iend - istart + 1
        if L <= 0:
            continue
        if L <= 2 * _FLANK:
            out.append((chrom, istart, iend))
        else:
            out.append((chrom, istart, istart + _FLANK - 1))
            out.append((chrom, iend - _FLANK + 1, iend))
    return out


def run() -> None:
    needed = _needed_transcripts()
    print(f"switch transcripts to scan (intronic): {len(needed):,}")

    tx = _transcript_exons(needed)
    multi = {t: r for t, r in tx.items() if len(r["exons"]) >= 2}
    print(f"  parsed GTF exons for {len(tx):,} transcripts "
          f"({len(multi):,} multi-exon with introns)")

    matrices, labels = _load_pwms()
    print(f"ATtRACT human PWMs: {len(matrices)} matrices, "
          f"{len({r for _, r in labels})} RBPs")

    bg = MOODS.tools.flat_bg(4)
    thresholds = [MOODS.tools.threshold_from_p(m, bg, _P_THRESH) for m in matrices]
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    genome = pyfaidx.Fasta(str(_GENOME_FA), sequence_always_upper=True)

    rows = []
    n = 0
    for tid, rec in multi.items():
        strand = rec["strand"]
        per_rbp: dict[str, int] = {}
        for chrom, gstart, gend in _intron_flanks(rec):
            if chrom not in genome:
                continue
            seq = str(genome[chrom][gstart - 1:gend])  # 1-based inclusive -> 0-based half-open
            if strand == "-":
                seq = _revcomp(seq)
            seq = seq.replace("U", "T")               # PWM U column already mapped to T
            for (_, rbp), matches in zip(labels, scanner.scan(seq)):
                c = len(matches)
                if c:
                    per_rbp[rbp] = per_rbp.get(rbp, 0) + c
        for rbp, c in per_rbp.items():
            rows.append((tid, rbp, c))
        n += 1
        if n % 5000 == 0:
            print(f"  scanned {n:,} transcripts")

    out = pd.DataFrame(rows, columns=["transcript_id", "rbp", "count"])
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(_OUT, index=False)
    print(f"scanned {n:,} transcripts -> {_OUT} ({len(out):,} (tx,RBP) hit rows)")


def main() -> None:
    argparse.ArgumentParser(
        description="Scan switch-isoform intronic splice-site flanks for ATtRACT RBP motifs."
    ).parse_args()
    run()


if __name__ == "__main__":
    main()
