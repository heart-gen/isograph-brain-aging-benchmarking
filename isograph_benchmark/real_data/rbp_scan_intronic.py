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

It inherits the mature scan's composition-matched background (reviewer item 6a) rather than a
flat 0.25 one, which matters more here than for mature transcripts: intronic sequence is
markedly AT-rich, so a flat background inflates exactly the AU-binding regulators this scope
exists to find. Backgrounds are binned over each transcript's pooled flank composition, and
the flat-background counts are emitted alongside under `bg_mode='flat'`. Hits are also tallied
per motif family. The `region` column is constant (`intron_flank`), so the schema matches.

Runs in the dedicated `motif` env (MOODS + pyfaidx + pandas).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyfaidx
import MOODS.scan
import MOODS.tools

from isograph_benchmark.paths import stage_out
from isograph_benchmark.real_data.rbp_scan import (
    _PSEUDOCOUNT,
    _P_THRESH,
    _bin_background,
    _load_pwms,
    _needed_transcripts,
    composition_bin,
)

_REPO = Path(__file__).resolve().parents[2]
_GENOME_FA = Path(
    "/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/"
    "fasta/GRCh38.primary_assembly.genome.fa")
_GTF = Path(
    "/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/"
    "gtf/gencode.v47.annotation.gtf")
_RBP_DIR = stage_out("regulation", "rbp")
_OUT = _RBP_DIR / "rbp_counts_intronic.parquet"
_FAM_OUT = _RBP_DIR / "rbp_family_counts_intronic.parquet"
_FAMILIES = _RBP_DIR / "rbp_motif_families.parquet"
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


def run(flat_too: bool = True) -> None:
    needed = _needed_transcripts()
    print(f"switch transcripts to scan (intronic): {len(needed):,}")

    tx = _transcript_exons(needed)
    multi = {t: r for t, r in tx.items() if len(r["exons"]) >= 2}
    print(f"  parsed GTF exons for {len(tx):,} transcripts "
          f"({len(multi):,} multi-exon with introns)")

    matrices, labels = _load_pwms()
    print(f"ATtRACT human PWMs: {len(matrices)} matrices, "
          f"{len({r for _, r in labels})} RBPs")

    if not _FAMILIES.exists():
        raise SystemExit(f"{_FAMILIES} not found; run rbp_motif_families.py first.")
    fam = pd.read_parquet(_FAMILIES, columns=["matrix_id", "family_id"])
    family_of = dict(zip(fam["matrix_id"], fam["family_id"]))

    genome = pyfaidx.Fasta(str(_GENOME_FA), sequence_always_upper=True)

    # Pass 1: gather each transcript's flank sequences and bin by their pooled composition.
    by_bin: dict[tuple[int, int], dict[str, list[str]]] = {}
    for tid, rec in multi.items():
        strand = rec["strand"]
        seqs = []
        for chrom, gstart, gend in _intron_flanks(rec):
            if chrom not in genome:
                continue
            seq = str(genome[chrom][gstart - 1:gend])  # 1-based inclusive -> 0-based half-open
            if strand == "-":
                seq = _revcomp(seq)
            seqs.append(seq.replace("U", "T"))         # PWM U column already mapped to T
        if not seqs:
            continue
        gb, pb, _, _ = composition_bin("".join(seqs))
        by_bin.setdefault((gb, pb), {})[tid] = seqs
    n = sum(len(v) for v in by_bin.values())
    print(f"  flank sequences for {n:,} transcripts in {len(by_bin)} composition bins",
          flush=True)

    rbp_rows: list = []
    fam_rows: list = []

    def scan_bins(bg_of, bg_mode: str) -> None:
        for k, (key, group) in enumerate(sorted(by_bin.items()), start=1):
            bg = bg_of(group)
            log_odds = [MOODS.tools.log_odds(m, bg, _PSEUDOCOUNT) for m in matrices]
            thresholds = [MOODS.tools.threshold_from_p(m, bg, _P_THRESH) for m in log_odds]
            scanner = MOODS.scan.Scanner(7)
            scanner.set_motifs(log_odds, bg, thresholds)
            for tid, seqs in group.items():
                per_rbp: dict[str, int] = {}
                per_fam: dict[str, int] = {}
                for seq in seqs:
                    for (mid, rbp), matches in zip(labels, scanner.scan(seq)):
                        c = len(matches)
                        if not c:
                            continue
                        per_rbp[rbp] = per_rbp.get(rbp, 0) + c
                        f = family_of.get(mid)
                        if f is not None:
                            per_fam[f] = per_fam.get(f, 0) + c
                for rbp, c in per_rbp.items():
                    rbp_rows.append((tid, rbp, "intron_flank", c, bg_mode))
                for f, c in per_fam.items():
                    fam_rows.append((tid, f, "intron_flank", c, bg_mode))
            print(f"  [{bg_mode}] bin {key} ({len(group):,} tx) [{k}/{len(by_bin)}]",
                  flush=True)

    scan_bins(lambda g: _bin_background([s for v in g.values() for s in v]), "composition")
    if flat_too:
        flat = MOODS.tools.flat_bg(4)
        scan_bins(lambda g: flat, "flat")

    out = pd.DataFrame(rbp_rows, columns=["transcript_id", "rbp", "region",
                                          "count", "bg_mode"])
    fam_out = pd.DataFrame(fam_rows, columns=["transcript_id", "family_id", "region",
                                              "count", "bg_mode"])
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(_OUT, index=False, compression="zstd")
    fam_out.to_parquet(_FAM_OUT, index=False, compression="zstd")
    print(f"scanned {n:,} transcripts -> {_OUT} ({len(out):,} rbp rows, "
          f"{len(fam_out):,} family rows)")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Scan switch-isoform intronic splice-site flanks for ATtRACT RBP motifs.")
    p.add_argument("--no-flat", action="store_true",
                   help="skip the parallel flat-background scan (halves runtime)")
    args = p.parse_args()
    run(not args.no_flat)


if __name__ == "__main__":
    main()
