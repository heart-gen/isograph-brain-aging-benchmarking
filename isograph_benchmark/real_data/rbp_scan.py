"""Stage 1 of the RBP-regulon analysis: scan switch-isoform sequences for RBP motifs.

Runs in the dedicated `motif` env (MOODS + pyfaidx), NOT the isograph env, so it depends
only on pandas + MOODS + gzip. It collects every transcript that appears in any region's
IsoGraph switch pairs, pulls its mature sequence from the GENCODE v47 transcript FASTA, and
scans it against the ATtRACT human RBP position-weight matrices, writing per-(transcript, RBP)
motif hit counts. Stage 2 (`rbp_regulon.py`, isograph env) turns those counts into the
within-pair "switch alters RBP X" call and the per-module regulon enrichment.

ATtRACT motifs are RNA (A C G U); the transcript FASTA is DNA (A C G T) and RBP motifs are
single-stranded, so we scan the sense strand only and treat the PWM's U column as T.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

import MOODS.scan
import MOODS.tools

_REPO = Path(__file__).resolve().parents[2]
_MOTIF_DIR = _REPO / "inputs" / "rbp_motifs"
_FASTA = _REPO / "inputs" / "raw" / "gencode_v47" / "gencode.v47.transcripts.fa.gz"
_OUT = _REPO / "real_data" / "_m" / "rbp" / "rbp_counts.parquet"
_P_THRESH = 1e-4          # per-position match p-value for a motif hit
_PSEUDOCOUNT = 0.1


def _needed_transcripts() -> set[str]:
    """Every transcript in any region's structure_switch_pairs (versioned IDs)."""
    tx: set[str] = set()
    for tree in ("brainseq", "gtex"):
        for f in (_REPO / "real_data" / tree).glob(
                "*/_m/isograph_vae/module_interpret/structure_switch_pairs.parquet"):
            d = pd.read_parquet(f, columns=["transcript_id_1", "transcript_id_2"])
            tx.update(d["transcript_id_1"])
            tx.update(d["transcript_id_2"])
    return tx


def _load_pwms() -> tuple[list, list[tuple[str, str]]]:
    """ATtRACT human PWMs as MOODS log-odds matrices + (matrix_id, rbp_name) labels."""
    db = pd.read_csv(_MOTIF_DIR / "ATtRACT_db.txt", sep="\t", dtype=str)
    db = db[db["Organism"] == "Homo_sapiens"]
    matrix_to_rbp = dict(zip(db["Matrix_id"], db["Gene_name"]))

    freqs: dict[str, list[list[float]]] = {}
    cur, rows = None, []
    for line in (_MOTIF_DIR / "pwm.txt").read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None:
                freqs[cur] = rows
            cur, rows = line[1:].split()[0], []
        else:
            rows.append([float(x) for x in line.split()])
    if cur is not None:
        freqs[cur] = rows

    bg = MOODS.tools.flat_bg(4)
    matrices, labels = [], []
    for mid, rws in freqs.items():
        rbp = matrix_to_rbp.get(mid)
        if rbp is None or len(rws) == 0:
            continue
        cols = list(zip(*rws))               # 4 columns: A, C, G, U(=T)
        freq = [list(cols[0]), list(cols[1]), list(cols[2]), list(cols[3])]
        matrices.append(MOODS.tools.log_odds(freq, bg, _PSEUDOCOUNT))
        labels.append((mid, rbp))
    return matrices, labels


def _iter_fasta(path: Path, needed: set[str]):
    """Yield (versioned transcript_id, sequence) for needed transcripts from a gzip FASTA."""
    tx_id, chunks = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if tx_id is not None and tx_id in needed:
                    yield tx_id, "".join(chunks)
                tx_id = line[1:].split("|", 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if tx_id is not None and tx_id in needed:
            yield tx_id, "".join(chunks)


def run() -> None:
    needed = _needed_transcripts()
    print(f"switch transcripts to scan: {len(needed):,}")
    matrices, labels = _load_pwms()
    rbp_of = [rbp for _, rbp in labels]
    print(f"ATtRACT human PWMs: {len(matrices)} matrices, {len(set(rbp_of))} RBPs")

    bg = MOODS.tools.flat_bg(4)
    thresholds = [MOODS.tools.threshold_from_p(m, bg, _P_THRESH) for m in matrices]
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    rows = []
    n = 0
    for tx_id, seq in _iter_fasta(_FASTA, needed):
        seq = seq.upper().replace("U", "T")
        results = scanner.scan(seq)
        # sum matrix hits per RBP (an RBP can have several matrices)
        per_rbp: dict[str, int] = {}
        for (_, rbp), matches in zip(labels, results):
            c = len(matches)
            if c:
                per_rbp[rbp] = per_rbp.get(rbp, 0) + c
        for rbp, c in per_rbp.items():
            rows.append((tx_id, rbp, c))
        n += 1
        if n % 5000 == 0:
            print(f"  scanned {n:,} transcripts")
    out = pd.DataFrame(rows, columns=["transcript_id", "rbp", "count"])
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(_OUT, index=False)
    print(f"scanned {n:,} transcripts -> {_OUT} ({len(out):,} (tx,RBP) hit rows)")


def main() -> None:
    argparse.ArgumentParser(description="Scan switch isoforms for ATtRACT RBP motifs.").parse_args()
    run()


if __name__ == "__main__":
    main()
