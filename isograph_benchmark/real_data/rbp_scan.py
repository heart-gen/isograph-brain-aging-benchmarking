"""Stage 1 of the RBP-regulon analysis: scan switch-isoform sequences for RBP motifs.

Runs in the dedicated `motif` env (MOODS + pandas), NOT the isograph env. It collects every
transcript that appears in any region's IsoGraph switch pairs, pulls its mature sequence from
the GENCODE v47 transcript FASTA, and scans it against the ATtRACT human position-weight
matrices. Stage 2 (`rbp_regulon.py`, isograph env) turns those counts into the within-pair
"switch alters RBP X" call and the per-module regulon enrichment.

ATtRACT motifs are RNA (A C G U); the transcript FASTA is DNA (A C G T) and RBP motifs are
single-stranded, so we scan the sense strand only and treat the PWM's U column as T.

Two robustness corrections (reviewer item 6a):

**Composition-matched background.** A flat 0.25 background makes an AU-rich motif look
enriched in an AU-rich transcript for no reason other than base composition, which
systematically inflates the ELAVL / CPEB / hnRNPD-type binders. Scoring thresholds are
therefore solved against a background matched to the transcript's own composition. A strictly
per-transcript background would need ~1,200 matrices x ~80,000 transcripts threshold solves,
which is not affordable, so transcripts are binned by GC and purine content (pinned in
``configs/rbp_families.yaml``) and one threshold set is solved per bin. Flat-background counts
are emitted alongside under ``bg_mode='flat'`` so the before/after is auditable in one table.

**Region partition.** The GENCODE transcript FASTA headers carry no UTR/CDS spans, so CDS is
projected from the GTF into mature-transcript coordinates and every hit is labelled
``5utr`` / ``cds`` / ``3utr`` / ``noncoding``. Motif gain or loss can then be compared within
a region class rather than across the whole mature transcript, where a switch that merely
lengthens the 3'UTR would masquerade as regulatory rewiring.

Counts are emitted at **two** resolutions. The previous schema summed all of an RBP's
matrices into one number, which is what makes N near-identical motifs read as N independent
lines of evidence, so ``rbp_family_counts.parquet`` additionally tallies by motif family
(from ``rbp_motif_families.py``). Matrix-level counts are *not* persisted: at ~1,200 matrices
x 4 regions x ~80,000 transcripts they run to ~10^8 rows, against ~2x10^7 for either
aggregate. Changing the family cut therefore requires a re-scan, which the pinned cut in
``configs/rbp_families.yaml`` makes a deliberate act rather than a routine one.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

import MOODS.scan
import MOODS.tools

_REPO = Path(__file__).resolve().parents[2]
_MOTIF_DIR = _REPO / "inputs" / "rbp_motifs"
_FASTA = _REPO / "inputs" / "raw" / "gencode_v47" / "gencode.v47.transcripts.fa.gz"
_OUT_DIR = _REPO / "real_data" / "_m" / "rbp"
_OUT = _OUT_DIR / "rbp_counts.parquet"
_FAM_OUT = _OUT_DIR / "rbp_family_counts.parquet"
_FAMILIES = _OUT_DIR / "rbp_motif_families.parquet"
_OPPORTUNITY = _OUT_DIR / "rbp_scan_opportunity.parquet"
_REGIONS_OUT = _OUT_DIR / "transcript_regions.parquet"

_P_THRESH = 1e-4          # per-position match p-value for a motif hit
_PSEUDOCOUNT = 0.1
_GC_BINS = 20
_PURINE_BINS = 2

REGION_CLASSES = ("5utr", "cds", "3utr", "noncoding")


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def _needed_transcripts() -> set[str]:
    """Every transcript in any region's structure_switch_pairs (versioned IDs)."""
    tx: set[str] = set()
    for tree in ("brainseq", "gtex"):
        root = _REPO / "real_data" / tree
        if not root.exists():
            continue
        for region_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            f = (region_dir / "_m" / "isograph_vae" / "module_interpret"
                 / "structure_switch_pairs.parquet")
            if not f.exists():
                continue
            d = pd.read_parquet(f, columns=["transcript_id_1", "transcript_id_2"])
            tx.update(d["transcript_id_1"])
            tx.update(d["transcript_id_2"])
    return tx


def _load_pwms() -> tuple[list[list[list[float]]], list[tuple[str, str]]]:
    """ATtRACT human PWMs as raw frequency matrices + (matrix_id, rbp_name) labels.

    Frequencies rather than log-odds, because the log-odds transform depends on the
    background and there is now one background per composition bin.
    """
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
            if cur is not None and rows:
                freqs[cur] = rows
            cur, rows = line[1:].split()[0], []
        else:
            rows.append([float(x) for x in line.split()])
    if cur is not None and rows:
        freqs[cur] = rows

    matrices, labels = [], []
    for mid, rws in freqs.items():
        rbp = matrix_to_rbp.get(mid)
        if rbp is None or not rws:
            continue
        cols = list(zip(*rws))               # 4 columns: A, C, G, U(=T)
        matrices.append([list(cols[0]), list(cols[1]), list(cols[2]), list(cols[3])])
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


# --------------------------------------------------------------------------- #
# Region partition: project genomic CDS spans into mature-transcript coordinates
# --------------------------------------------------------------------------- #
def transcript_region_spans(exons: list[tuple[int, int]], cds: list[tuple[int, int]],
                            strand: str) -> list[tuple[int, int, str]]:
    """``[(start, end, region_class)]`` in 0-based half-open mature-transcript coordinates.

    Exons are concatenated in transcription order (reverse genomic order on the minus
    strand). With no annotated CDS the whole transcript is ``noncoding``; otherwise the
    projected CDS splits it into 5'UTR / CDS / 3'UTR, and the three spans tile the mature
    length exactly.
    """
    if not exons:
        return []
    ordered = sorted(exons, key=lambda e: e[0], reverse=(strand == "-"))
    length = sum(e - s for s, e in ordered)
    if not cds:
        return [(0, length, "noncoding")]

    # genomic -> transcript offset for the CDS boundaries
    cds_lo = min(s for s, _ in cds)
    cds_hi = max(e for _, e in cds)

    def to_tx(pos: int) -> int | None:
        off = 0
        for s, e in ordered:
            if s <= pos < e:
                return off + (pos - s if strand == "+" else e - 1 - pos)
            off += e - s
        return None

    a = to_tx(cds_lo if strand == "+" else cds_hi - 1)
    b = to_tx(cds_hi - 1 if strand == "+" else cds_lo)
    if a is None or b is None:
        return [(0, length, "noncoding")]
    start, end = min(a, b), max(a, b) + 1

    spans = []
    if start > 0:
        spans.append((0, start, "5utr"))
    spans.append((start, end, "cds"))
    if end < length:
        spans.append((end, length, "3utr"))
    return spans


def _region_lookup(spans: list[tuple[int, int, str]], length: int) -> np.ndarray:
    """Per-base region-class code, so a hit position maps to its class in O(1)."""
    codes = np.zeros(length, dtype=np.int8)
    index = {c: i for i, c in enumerate(REGION_CLASSES)}
    for s, e, cls in spans:
        codes[s:min(e, length)] = index[cls]
    return codes


def load_transcript_regions(needed: set[str], gtf_cache: Path) -> dict[str, list]:
    """transcript_id -> region spans, read from the GTF parquet cache.

    The cache is written by ``isograph.explain.structure.parse_gtf`` in the isograph env;
    this stage runs in the motif env and only reads it, so the heavy GTF parse is not
    repeated and the two envs stay decoupled.
    """
    if not gtf_cache.exists():
        raise SystemExit(
            f"GTF cache not found: {gtf_cache}\nGenerate it once from the isograph env:\n"
            f"  python -c 'from isograph.explain.structure import parse_gtf; "
            f"from isograph_benchmark.real_data.interpret_modules import "
            f"DEFAULT_GTF_PATH, DEFAULT_GTF_CACHE; "
            f"parse_gtf(DEFAULT_GTF_PATH, cache=DEFAULT_GTF_CACHE)'")
    df = pd.read_parquet(gtf_cache,
                         columns=["transcript_id", "strand", "feature", "start", "end"])
    df = df[df["transcript_id"].isin(needed) & df["feature"].isin(["exon", "CDS"])]

    strand = (df.drop_duplicates("transcript_id")
              .set_index("transcript_id")["strand"].astype(str).to_dict())
    # The GTF cache carries GTF coordinates: 1-based, inclusive on both ends. Everything
    # below is 0-based half-open, so a raw (start, end) would drop one base per exon and the
    # region spans would not tile the mature length.
    df = df.assign(start=df["start"] - 1)
    grouped = {
        feat: {tx: list(zip(g["start"], g["end"]))
               for tx, g in sub.groupby("transcript_id", observed=True)}
        for feat, sub in df.groupby("feature", observed=True)
    }
    exons_by_tx = grouped.get("exon", {})
    cds_by_tx = grouped.get("CDS", {})

    out: dict[str, list] = {}
    for tx, exons in exons_by_tx.items():
        out[tx] = transcript_region_spans(exons, cds_by_tx.get(tx, []), strand.get(tx, "+"))
    return out


# --------------------------------------------------------------------------- #
# Composition binning
# --------------------------------------------------------------------------- #
def composition_bin(seq: str, gc_bins: int = _GC_BINS,
                    purine_bins: int = _PURINE_BINS) -> tuple[int, int, float, float]:
    """(gc_bin, purine_bin, gc, purine) for one sequence."""
    n = len(seq) or 1
    counts = {b: seq.count(b) for b in "ACGT"}
    gc = (counts["G"] + counts["C"]) / n
    purine = (counts["A"] + counts["G"]) / n
    gb = int(np.clip(int(gc * gc_bins), 0, gc_bins - 1))
    pb = int(np.clip(int(purine * purine_bins), 0, purine_bins - 1))
    return gb, pb, gc, purine


def _bin_background(seqs: list[str]) -> tuple[float, float, float, float]:
    """Mean A/C/G/T frequency over a composition bin's sequences."""
    counts = np.zeros(4)
    for s in seqs:
        counts += [s.count("A"), s.count("C"), s.count("G"), s.count("T")]
    total = counts.sum()
    if total <= 0:
        return MOODS.tools.flat_bg(4)
    return tuple(counts / total)


# --------------------------------------------------------------------------- #
# Scan
# --------------------------------------------------------------------------- #
def _scan_group(seqs: dict[str, str], matrices, labels, bg, regions: dict[str, list],
                bg_mode: str, family_of: dict[str, str],
                rbp_rows: list, fam_rows: list) -> None:
    """Scan one composition bin with a single threshold set, tallying hits by region.

    Matrix-level hits are folded into the per-RBP and per-family tallies as they are
    produced, so the ~10^8-row matrix-level intermediate never has to exist.
    """
    log_odds = [MOODS.tools.log_odds(m, bg, _PSEUDOCOUNT) for m in matrices]
    thresholds = [MOODS.tools.threshold_from_p(m, bg, _P_THRESH) for m in log_odds]
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(log_odds, bg, thresholds)

    for tx_id, seq in seqs.items():
        spans = regions.get(tx_id) or [(0, len(seq), "noncoding")]
        codes = _region_lookup(spans, len(seq))
        results = scanner.scan(seq)
        by_rbp: dict[tuple[str, str], int] = {}
        by_family: dict[tuple[str, str], int] = {}
        for (mid, rbp), matches in zip(labels, results):
            if not matches:
                continue
            pos = np.fromiter((m.pos for m in matches), dtype=np.int64, count=len(matches))
            pos = np.clip(pos, 0, len(seq) - 1)
            fam = family_of.get(mid)
            for code, c in zip(*np.unique(codes[pos], return_counts=True)):
                region = REGION_CLASSES[code]
                key = (rbp, region)
                by_rbp[key] = by_rbp.get(key, 0) + int(c)
                if fam is not None:
                    fkey = (fam, region)
                    by_family[fkey] = by_family.get(fkey, 0) + int(c)
        for (rbp, region), c in by_rbp.items():
            rbp_rows.append((tx_id, rbp, region, c, bg_mode))
        for (fam, region), c in by_family.items():
            fam_rows.append((tx_id, fam, region, c, bg_mode))


def run(gtf_cache: Path, gc_bins: int, purine_bins: int, flat_too: bool,
        limit: int | None) -> None:
    needed = _needed_transcripts()
    print(f"switch transcripts to scan: {len(needed):,}", flush=True)
    matrices, labels = _load_pwms()
    print(f"ATtRACT human PWMs: {len(matrices)} matrices, "
          f"{len({r for _, r in labels})} RBPs", flush=True)

    regions = load_transcript_regions(needed, gtf_cache)
    print(f"region spans resolved for {len(regions):,} transcripts", flush=True)

    # Pass 1: read sequences, bin by composition, record the opportunity covariates.
    by_bin: dict[tuple[int, int], dict[str, str]] = {}
    opportunity = []
    for i, (tx_id, seq) in enumerate(_iter_fasta(_FASTA, needed)):
        if limit is not None and i >= limit:
            break
        seq = seq.upper().replace("U", "T")
        gb, pb, gc, purine = composition_bin(seq, gc_bins, purine_bins)
        by_bin.setdefault((gb, pb), {})[tx_id] = seq
        spans = regions.get(tx_id) or [(0, len(seq), "noncoding")]
        for s, e, cls in spans:
            sub = seq[s:e]
            opportunity.append((tx_id, cls, e - s, (sub.count("G") + sub.count("C")) / max(len(sub), 1),
                                gb, pb, gc, purine, len(seq)))
    n_tx = sum(len(v) for v in by_bin.values())
    print(f"read {n_tx:,} sequences into {len(by_bin)} composition bins", flush=True)

    if not _FAMILIES.exists():
        raise SystemExit(
            f"{_FAMILIES} not found; run rbp_motif_families.py (isograph env) first so hits "
            f"can be tallied per motif family as well as per RBP.")
    fam = pd.read_parquet(_FAMILIES, columns=["matrix_id", "family_id"])
    family_of = dict(zip(fam["matrix_id"], fam["family_id"]))

    rbp_rows: list = []
    fam_rows: list = []
    for k, (key, seqs) in enumerate(sorted(by_bin.items()), start=1):
        bg = _bin_background(list(seqs.values()))
        _scan_group(seqs, matrices, labels, bg, regions, "composition", family_of,
                    rbp_rows, fam_rows)
        print(f"  bin {key} ({len(seqs):,} tx) bg=({bg[0]:.3f},{bg[1]:.3f},"
              f"{bg[2]:.3f},{bg[3]:.3f})  [{k}/{len(by_bin)}]", flush=True)

    if flat_too:
        flat = MOODS.tools.flat_bg(4)
        for k, (key, seqs) in enumerate(sorted(by_bin.items()), start=1):
            _scan_group(seqs, matrices, labels, flat, regions, "flat", family_of,
                        rbp_rows, fam_rows)
            print(f"  flat-bg bin {key} [{k}/{len(by_bin)}]", flush=True)

    out = pd.DataFrame(rbp_rows, columns=["transcript_id", "rbp", "region",
                                          "count", "bg_mode"])
    fam_out = pd.DataFrame(fam_rows, columns=["transcript_id", "family_id", "region",
                                              "count", "bg_mode"])
    opp = pd.DataFrame(opportunity, columns=["transcript_id", "region", "length", "gc",
                                             "gc_bin", "purine_bin", "tx_gc", "tx_purine",
                                             "tx_length"])
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_parquet(_OUT, index=False, compression="zstd")
    fam_out.to_parquet(_FAM_OUT, index=False, compression="zstd")
    opp.to_parquet(_OPPORTUNITY, index=False, compression="zstd")
    pd.DataFrame(
        [(t, s, e, c) for t, sp in regions.items() for s, e, c in sp],
        columns=["transcript_id", "start", "end", "region"],
    ).to_parquet(_REGIONS_OUT, index=False)

    print(f"scanned {n_tx:,} transcripts -> {_OUT} ({len(out):,} rbp rows, "
          f"{len(fam_out):,} family rows)", flush=True)
    for mode, g in out.groupby("bg_mode"):
        print(f"  {mode}: {int(g['count'].sum()):,} hits over "
              f"{g['rbp'].nunique()} RBPs", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Scan switch isoforms for ATtRACT RBP motifs.")
    p.add_argument("--gtf-cache", default=str(
        _REPO / "real_data" / "_m" / "tmp"
        / "gencode.v47.primary_assembly.annotation.gtf_cache.parquet"))
    p.add_argument("--gc-bins", type=int, default=_GC_BINS)
    p.add_argument("--purine-bins", type=int, default=_PURINE_BINS)
    p.add_argument("--no-flat", action="store_true",
                   help="skip the parallel flat-background scan (halves runtime, but the "
                        "before/after comparison is then not reproducible from one table)")
    p.add_argument("--limit", type=int, default=None, help="debug: scan only N transcripts")
    args = p.parse_args()
    run(Path(args.gtf_cache), args.gc_bins, args.purine_bins, not args.no_flat, args.limit)


if __name__ == "__main__":
    main()
