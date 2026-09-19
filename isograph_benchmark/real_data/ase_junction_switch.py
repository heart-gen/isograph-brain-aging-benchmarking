"""Allele-aware junction recount of switch pairs, from the WASP BAMs on Quest (PI item 10a).

WHAT THIS IS FOR
----------------
The phASER screen (`ase_switch_direction.py --stage screen`) showed that exon-unique
sequence is a ONE-SIDED discriminator: for most switch pairs only one isoform has unique
exonic sequence carrying a heterozygous site, so there is nothing to contrast. An
isoform-specific splice JUNCTION is two-sided by construction -- T1 splices one intron, T2
another -- so a fragment that crosses a specific junction AND carries a phased
heterozygous site is informative for isoform and haplotype at once.

This module counts those fragments from the ASE-grade BAMs (STAR two-pass,
`--waspOutputMode SAMtag`) and applies the SAME pre-registered gate as the exonic screen.
It runs on Quest, where the BAMs are; only the count table and the report come back.

  --stage targets   Isoform-specific junctions per switch pair (GENCODE v47), the all_samples
                    and ea_only switch-QTL lead per gene, the donors' phased genotype at each
                    lead, the fetch regions and the sample list for the array.
  --stage count     One BAM: fragment-level isoform x haplotype counts per pair, plus
                    isoform-agnostic per-site haplotype counts for the concordance check.
  --stage concordance
                    One sample: per-site haplotype counts from this recount against the
                    phASER release's ASEReadCounter counts at the same sites.
  --stage screen    Pool the shards into `junction_allelic_counts.parquet` and apply the gate.

WHAT IS COUNTED
---------------
A fragment (both mates, grouped by read name) is kept when:
  * every alignment is primary, mapped, not QC-fail, and uniquely mapped (`NH:i:1`);
  * no mate carries a WASP tag other than `vW:i:1` (reference mapping bias is the dominant
    confound in a within-donor allelic contrast; a fragment with a failing mate is dropped
    whole), and every mate that covers a heterozygous site carries `vW:i:1`;
  * its junctions (CIGAR `N`, anchored by >= `--min-anchor` aligned bases on both sides)
    hit at least one junction specific to one isoform of the pair, none specific to the
    other, and every junction it carries inside the pair's span is an intron of the
    isoform it is assigned to. A fragment with a junction unknown to the assigned isoform
    is `incompatible` and not counted.
Haplotype comes from the base at each phased heterozygous SNV the fragment covers (base
quality >= 10, the release's ASEReadCounter setting), phased by phASER's genome-wide
genotype `PW`. When both mates cover a site they must agree; all sites on a fragment must
name the same haplotype, otherwise it is a `hap_conflict` (kept as an error-rate
estimate). A fragment with no heterozygous site is recorded as hap 0: it says which
isoform but not which allele, and is kept only for isoform totals.

STAGE 4 -- THE ALLELIC TEST -- LIVES IN `ase_junction_allelic.py`
----------------------------------------------------------------
It runs only where this module's gate passed. It orients each lead-heterozygous donor to
the haplotype carrying the lead's ALT allele (`PW` is anchored to the same population
phasing as the genotype VCF, so the lead's phased GT from `lead_genotypes.parquet` is in
the same frame). It then fits a beta-binomial GLMM with a donor random intercept per pair.
Donors homozygous at the lead are the built-in null. The lead-to-site distance
(`<sample>.pair_sites.tsv.gz`) records exposure to phase switches between the lead and a
fragment's heterozygous site.

CAVEATS THE REPORT CARRIES
--------------------------
Cell composition cancels within a donor only if allelic effects do not differ by cell
type; some reference mapping bias survives WASP.
"""
from __future__ import annotations

import argparse
import bisect
import json
import subprocess
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.ase_switch_direction import (
    MIN_DONORS_INFORMATIVE, MIN_READS_PER_DONOR, REGIONS, _merge)

# A separate cache from the `_GTF_CACHE` other modules share: that file is written by
# `isograph.explain.structure.annotate_switch_pairs` in its own schema, and a lean exon
# table written to the same path could be read back by a consumer expecting more.
_EXON_CACHE = stage_out("tmp", "gencode.v47.primary_assembly.exons.parquet")
_N_TRANSCRIPTS_V47 = 387_944  # see inputs/build_transcript_annotation.py

MIN_BASE_QUALITY = 10
MIN_ANCHOR = 8
QVAL_GATE_FAMILY = 0.05
# Pre-specified: fewer passing gate-family pairs than this drops the allele-specific arm.
MIN_PAIRS_TO_PROCEED = 30

# BAM flags that disqualify an alignment outright.
_BAD_FLAGS = 0x4 | 0x100 | 0x200 | 0x800


def quest_config(region: str) -> dict:
    cfg = load_yaml("configs/data_sources.yaml").get("quest", {}).get("ase_junction", {})
    if region not in cfg.get("regions", {}):
        raise SystemExit(f"configs/data_sources.yaml has no quest.ase_junction.regions.{region}")
    out = dict(cfg["regions"][region])
    out["genotype_vcf_template"] = cfg["genotype_vcf_template"]
    out["gtf"] = cfg["gtf"]
    return out


def out_dir(region: str) -> Path:
    return ensure_dir(stage_out("mechanism", "ase_junction_switch", region))


def _store(region: str) -> Path:
    from isograph_benchmark.paths import region_store
    return region_store("brainseq", region)


# --------------------------------------------------------------------------- #
# Annotation: exons -> introns -> isoform-specific junctions
# --------------------------------------------------------------------------- #
def load_exons(gtf: Path) -> pd.DataFrame:
    if _EXON_CACHE.exists():
        return pd.read_parquet(_EXON_CACHE)
    rows = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "exon":
                continue
            attr = f[8]
            i = attr.find('transcript_id "')
            j = attr.find('gene_id "')
            if i < 0 or j < 0:
                continue
            tx = attr[i + 15: attr.find('"', i + 15)]
            gene = attr[j + 9: attr.find('"', j + 9)]
            rows.append((tx, gene, f[0], f[6], "exon", int(f[3]), int(f[4])))
    ex = pd.DataFrame(rows, columns=["transcript_id", "gene_id", "chrom", "strand",
                                     "feature", "start", "end"])
    n = ex["transcript_id"].nunique()
    if n != _N_TRANSCRIPTS_V47:
        raise SystemExit(f"{gtf} carries {n:,} transcripts, expected {_N_TRANSCRIPTS_V47:,} "
                         f"(GENCODE v47 primary assembly); wrong annotation file?")
    ensure_dir(_EXON_CACHE.parent)
    ex.to_parquet(_EXON_CACHE, index=False)
    return ex


def transcript_introns(exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """1-based closed introns between consecutive exons -- CIGAR `N` coordinates."""
    ex = _merge(sorted(exons))
    return [(a[1] + 1, b[0] - 1) for a, b in zip(ex, ex[1:])]


def specific_junctions(i1: list[tuple[int, int]], i2: list[tuple[int, int]]
                       ) -> tuple[set, set]:
    s1, s2 = set(i1), set(i2)
    return s1 - s2, s2 - s1


def build_junctions(pairs: pd.DataFrame, exons: pd.DataFrame
                    ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per pair: specific junctions of each side, and the full intron sets.

    Returns (pairs_out, junctions). `junctions` has one row per (pair, side, intron) for
    EVERY intron of either isoform, with `specific` marking the discriminating ones; the
    full set is what the compatibility check in `count` needs.
    """
    by_tx = {t: g for t, g in exons.groupby("transcript_id", sort=False)}
    prow, jrow = [], []
    for r in pairs.itertuples(index=False):
        t1, t2 = str(r.transcript_id_1), str(r.transcript_id_2)
        pair_id = f"{t1}|{t2}"
        if t1 not in by_tx or t2 not in by_tx:
            prow.append({"pair_id": pair_id, "gene_id": r.gene_id, "transcript_id_1": t1,
                         "transcript_id_2": t2, "in_gtf": False})
            continue
        e1, e2 = by_tx[t1], by_tx[t2]
        chrom, strand = str(e1["chrom"].iloc[0]), str(e1["strand"].iloc[0])
        i1 = transcript_introns(list(zip(e1["start"], e1["end"])))
        i2 = transcript_introns(list(zip(e2["start"], e2["end"])))
        sp1, sp2 = specific_junctions(i1, i2)
        span = (int(min(e1["start"].min(), e2["start"].min())),
                int(max(e1["end"].max(), e2["end"].max())))
        prow.append({"pair_id": pair_id, "gene_id": r.gene_id, "transcript_id_1": t1,
                     "transcript_id_2": t2, "in_gtf": True, "chrom": chrom,
                     "strand": strand, "span_start": span[0], "span_end": span[1],
                     "n_introns_1": len(i1), "n_introns_2": len(i2),
                     "n_specific_1": len(sp1), "n_specific_2": len(sp2),
                     "two_sided": bool(sp1) and bool(sp2)})
        for side, introns, spec in ((1, i1, sp1), (2, i2, sp2)):
            for s, e in introns:
                jrow.append({"pair_id": pair_id, "isoform": side, "chrom": chrom,
                             "start": s, "end": e, "specific": (s, e) in spec})
    return pd.DataFrame(prow), pd.DataFrame(jrow)


# --------------------------------------------------------------------------- #
# Stage: targets
# --------------------------------------------------------------------------- #
def _leads(region: str, arm: str) -> pd.DataFrame:
    f = stage_out("anchoring.brainseq_qtl", arm, region, "qtl", "cis_qtl_switch.parquet")
    if not f.exists():
        raise SystemExit(f"missing switch-QTL table: {f}")
    q = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se",
                                    "pval_nominal", "qval", "af"])
    return q.rename(columns={"phenotype_id": "gene_id"})


def _coloc_genes() -> tuple[set[str], set[str]]:
    b = stage_out("anchoring", "coloc_brainseq", "ea_only", "nominations.parquet")
    brainseq = set(pd.read_parquet(b)["gene"].astype(str)) if b.exists() else set()
    audit = set()
    for sub in ("abf", "susie"):
        f = stage_out("anchoring", "locus_event_audit", sub, "locus_event_audit.parquet")
        if f.exists():
            audit |= set(pd.read_parquet(f)["gene"].astype(str))
    return brainseq, audit


def _bcftools(args: list[str]) -> str:
    p = subprocess.run(["bcftools", *args], capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"bcftools {' '.join(args[:3])} ... failed: {p.stderr[-2000:]}")
    return p.stdout


def lead_genotypes(leads: pd.DataFrame, pairs: pd.DataFrame, template: str,
                   donors: list[str], dest: Path, window: int = 1_100_000) -> pd.DataFrame:
    """Phased genotype of every donor at every lead rsID, by one query per chromosome.

    The rsID is resolved inside the gene's cis window (tensorQTL's 1 Mb plus slack), so
    the query touches only the indexed blocks it needs instead of a whole chromosome.
    """
    span = pairs[pairs["in_gtf"]].groupby("gene_id").agg(
        chrom=("chrom", "first"), s=("span_start", "min"), e=("span_end", "max"))
    rows = []
    for chrom, g in leads.merge(span, left_on="gene_id", right_index=True).groupby("chrom"):
        vcf = template.format(chrom=chrom)
        if not Path(vcf).exists():
            print(f"  no genotype VCF for {chrom}; {len(g)} leads unresolved")
            continue
        bed = dest / f"_lead_windows.{chrom}.bed"
        ids = dest / f"_lead_ids.{chrom}.txt"
        win = [(max(0, int(r.s) - window), int(r.e) + window) for r in g.itertuples()]
        bed.write_text("".join(f"{chrom}\t{a}\t{b}\n" for a, b in _merge(win)))
        ids.write_text("\n".join(sorted(set(g["variant_id"]))) + "\n")
        smp = dest / "_donors.txt"
        smp.write_text("\n".join(donors) + "\n")
        out = _bcftools(["query", "-R", str(bed), "-S", str(smp), "--force-samples",
                         "-i", f"ID=@{ids}",
                         "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GT]\n", vcf])
        for line in out.splitlines():
            f = line.split("\t")
            if len(f[3]) != 1 or len(f[4]) != 1:
                continue  # biallelic SNV leads only; indels cannot be read off a base
            for cell in f[5:]:
                donor, gt = cell.rsplit("=", 1)
                rows.append((f[2], f[0], int(f[1]), f[3], f[4], donor, gt))
        for tmp in (bed, ids):
            tmp.unlink(missing_ok=True)
    (dest / "_donors.txt").unlink(missing_ok=True)
    gt = pd.DataFrame(rows, columns=["variant_id", "chrom", "pos", "ref", "alt",
                                     "genotype_sample_id", "gt"])
    return gt.drop_duplicates(["variant_id", "genotype_sample_id"])


def _quickcheck(bams: list[str]) -> list[str]:
    """BAMs failing `samtools quickcheck` (missing EOF block, unreadable header)."""
    p = subprocess.run(["samtools", "quickcheck", "-v", *bams], capture_output=True, text=True)
    return [line for line in p.stdout.splitlines() if line.strip()]


def run_targets(region: str) -> Path:
    cfg = quest_config(region)
    dest = out_dir(region)
    sp = _store(region) / "isograph_vae" / "module_interpret" / "structure_switch_pairs.parquet"
    if not sp.exists():
        raise SystemExit(f"missing switch pairs: {sp}")
    pairs_in = pd.read_parquet(sp)
    exons = load_exons(Path(cfg["gtf"]))
    pairs, junctions = build_junctions(pairs_in, exons)
    print(f"  {len(pairs):,} pairs; {int(pairs['in_gtf'].sum()):,} in the GTF; "
          f"{int(pairs['two_sided'].fillna(False).sum()):,} two-sided")

    # leads: all_samples is primary (within-donor contrasts cancel stratification), ea_only
    # is recorded alongside as a sensitivity check.
    lead = _leads(region, "all_samples").add_suffix("_all").rename(
        columns={"gene_id_all": "gene_id"})
    ea = _leads(region, "ea_only")[["gene_id", "variant_id", "qval"]].rename(
        columns={"variant_id": "variant_id_ea", "qval": "qval_ea"})
    pairs = pairs.merge(lead, on="gene_id", how="left").merge(ea, on="gene_id", how="left")
    pairs["gate_family"] = pairs["qval_all"] < QVAL_GATE_FAMILY
    pairs["lead_differs_ea"] = (pairs["variant_id_ea"].notna()
                                & (pairs["variant_id_ea"] != pairs["variant_id_all"]))
    brainseq, audit = _coloc_genes()
    bare = pairs["gene_id"].astype(str).str.split(".").str[0]
    pairs["coloc_brainseq"] = bare.isin(brainseq)
    pairs["coloc_audit"] = bare.isin(audit)

    man = pd.read_csv(cfg["manifest"], sep="\t",
                      usecols=["sample_id", "donor_id", "genotype_sample_id"])
    bam_dir = Path(cfg["bam_dir"])
    man["bam"] = [str(bam_dir / cfg["bam_pattern"].format(sample=s)) for s in man["sample_id"]]
    man["phaser_vcf"] = [str(Path(cfg["phaser_vcf_dir"]) / f"{s}.vcf.gz")
                         for s in man["sample_id"]]
    ok = (man["bam"].map(lambda p: Path(p).exists())
          & man["bam"].map(lambda p: Path(p + ".bai").exists())
          & man["phaser_vcf"].map(lambda p: Path(p).exists()))
    print(f"  samples: {int(ok.sum())} of {len(man)} have a BAM, its index and a phASER VCF "
          f"(missing: {', '.join(man.loc[~ok, 'sample_id']) or 'none'})")
    man = man[ok].reset_index(drop=True)
    # BAMs are staged onto Quest one region at a time; a truncated copy would fail a whole
    # array task hours in, so exclude anything samtools cannot read to its EOF marker.
    bad = _quickcheck(list(man["bam"]))
    if bad:
        print(f"  EXCLUDED {len(bad)} BAM(s) failing samtools quickcheck: "
              f"{', '.join(Path(b).name for b in bad)}")
        man = man[~man["bam"].isin(bad)].reset_index(drop=True)
    man.to_csv(dest / "samples.tsv", sep="\t", index=False)

    leads = pd.concat([
        pairs[["gene_id", "variant_id_all"]].rename(columns={"variant_id_all": "variant_id"}),
        pairs[["gene_id", "variant_id_ea"]].rename(columns={"variant_id_ea": "variant_id"}),
    ]).dropna().drop_duplicates()
    leads = leads[leads["gene_id"].isin(pairs.loc[pairs["in_gtf"], "gene_id"])]
    gt = lead_genotypes(leads, pairs, cfg["genotype_vcf_template"],
                        sorted(man["genotype_sample_id"].astype(str).unique()), dest)
    gt = gt.merge(man[["sample_id", "donor_id", "genotype_sample_id"]],
                  on="genotype_sample_id", how="inner")
    gt.to_parquet(dest / "lead_genotypes.parquet", index=False)
    res = set(gt["variant_id"])
    need = set(pairs.loc[pairs["gate_family"] & pairs["in_gtf"], "variant_id_all"])
    print(f"  gate-family leads resolved in the phased VCF: {len(need & res)} of {len(need)}")
    pairs["lead_resolved"] = pairs["variant_id_all"].isin(res)

    pairs.to_parquet(dest / "pairs.parquet", index=False)
    junctions.to_parquet(dest / "junctions.parquet", index=False)
    reg = pairs[pairs["in_gtf"]]
    bed = []
    for chrom, g in reg.groupby("chrom"):
        for a, b in _merge(list(zip(g["span_start"].astype(int), g["span_end"].astype(int)))):
            bed.append(f"{chrom}\t{a - 1}\t{b}\n")
    (dest / "regions.bed").write_text("".join(bed))
    print(f"  {len(bed):,} fetch regions; wrote {dest}")
    return dest


# --------------------------------------------------------------------------- #
# Stage: count -- read-level logic (pure functions, unit-tested)
# --------------------------------------------------------------------------- #
def read_junctions(cigartuples, ref_start: int, min_anchor: int = MIN_ANCHOR
                   ) -> list[tuple[int, int]]:
    """1-based introns spanned by `N` operations, each anchored on both sides.

    The anchor is the aligned reference length of the contiguous block (M/=/X/D, broken
    only by `N`) on each side of the intron. Short anchors are where spurious STAR
    junctions live, and one misplaced overhang would assign a fragment to the wrong
    isoform.
    """
    blocks: list[int] = []      # aligned length of each N-delimited block
    introns: list[tuple[int, int]] = []
    pos, cur = ref_start, 0
    for op, n in cigartuples:
        if op in (0, 2, 7, 8):
            pos += n
            cur += n
        elif op == 3:
            introns.append((pos + 1, pos + n))
            blocks.append(cur)
            cur = 0
            pos += n
    blocks.append(cur)
    return [iv for k, iv in enumerate(introns)
            if blocks[k] >= min_anchor and blocks[k + 1] >= min_anchor]


@dataclass
class HetSites:
    """Phased heterozygous SNVs on one chromosome: sorted positions + the two haplotype
    alleles at each (1-based)."""
    pos: list[int] = field(default_factory=list)
    hap_alleles: dict[int, tuple[str, str]] = field(default_factory=dict)
    ref_alt: dict[int, tuple[str, str]] = field(default_factory=dict)


def read_het_obs(read, sites: HetSites, min_bq: int = MIN_BASE_QUALITY
                 ) -> dict[int, int]:
    """{site: hap} for every phased het site this alignment covers with an aligned base.

    hap is 1 or 2 when the base matches one haplotype's allele, 0 when it matches neither
    (sequencing error or a third allele); sites under a low-quality base are skipped.
    """
    if not sites.pos:
        return {}
    hits = []
    for a, b in read.get_blocks():          # 0-based half-open aligned blocks
        i = bisect.bisect_left(sites.pos, a + 1)
        while i < len(sites.pos) and sites.pos[i] <= b:
            hits.append(sites.pos[i])
            i += 1
    if not hits:
        return {}
    want = {p - 1 for p in hits}
    qmap = {r: q for q, r in read.get_aligned_pairs(matches_only=True) if r in want}
    seq, qual = read.query_sequence, read.query_qualities
    out = {}
    for p in hits:
        q = qmap.get(p - 1)
        if q is None or seq is None:
            continue
        if qual is not None and qual[q] < min_bq:
            continue
        h1, h2 = sites.hap_alleles[p]
        base = seq[q]
        out[p] = 1 if base == h1 else 2 if base == h2 else 0
    return out


@dataclass
class ReadRec:
    junctions: tuple
    het: dict            # site -> hap (1/2/0)
    vw: int | None       # WASP tag value, None when the read carries none


def merge_mates(mates: list[ReadRec]) -> tuple[str, frozenset, dict]:
    """Combine mates into one fragment: (status, junctions, site -> hap).

    status is `wasp_fail` when any mate carries a failing WASP tag, `wasp_missing` when a
    mate covers a heterozygous site without any WASP tag (STAR tags every read overlapping
    a variant it was given, so this should not happen), else `ok`. A site covered by both
    mates with disagreeing haplotypes, or by a base matching neither allele, is dropped.
    """
    if any(m.vw is not None and m.vw != 1 for m in mates):
        return "wasp_fail", frozenset(), {}
    if any(m.het and m.vw is None for m in mates):
        return "wasp_missing", frozenset(), {}
    junc = frozenset(j for m in mates for j in m.junctions)
    seen: dict[int, set] = defaultdict(set)
    for m in mates:
        for p, h in m.het.items():
            seen[p].add(h)
    sites = {p: next(iter(h)) for p, h in seen.items() if len(h) == 1 and 0 not in h}
    return "ok", junc, sites


def fragment_hap(sites: dict[int, int]) -> int | str:
    """0 = no phased site; 1/2 = haplotype; 'conflict' = sites disagree."""
    haps = set(sites.values())
    if not haps:
        return 0
    return haps.pop() if len(haps) == 1 else "conflict"


@dataclass
class PairModel:
    pair_id: str
    span: tuple[int, int]
    introns: tuple[frozenset, frozenset]
    specific: tuple[frozenset, frozenset]


def assign_isoform(junc: frozenset, pm: PairModel) -> int | str | None:
    """1 or 2 for an informative fragment, None when it hits no specific junction,
    'conflict' when it hits both sides, 'incompatible' when it also carries a junction
    (inside the pair's span) that the assigned isoform does not splice."""
    h1 = bool(junc & pm.specific[0])
    h2 = bool(junc & pm.specific[1])
    if not (h1 or h2):
        return None
    if h1 and h2:
        return "conflict"
    iso = 1 if h1 else 2
    inside = {j for j in junc if j[0] >= pm.span[0] and j[1] <= pm.span[1]}
    return iso if inside <= pm.introns[iso - 1] else "incompatible"


# --------------------------------------------------------------------------- #
# Stage: count -- driver
# --------------------------------------------------------------------------- #
def load_het_sites(vcf_path: str, regions: list[tuple[str, int, int]]
                   ) -> tuple[dict[str, HetSites], Counter]:
    """Phased heterozygous SNVs in the regions, from the sample's phASER VCF.

    Phase comes from `PW` (phASER's genome-wide genotype, i.e. the gw-phased call) and
    falls back to `GT` only when PW is absent; the fallback count is reported.
    """
    import pysam
    qc: Counter = Counter()
    out: dict[str, HetSites] = defaultdict(HetSites)
    with pysam.VariantFile(vcf_path) as vf:
        smp = list(vf.header.samples)[0]
        has_pw = "PW" in vf.header.formats
        for chrom, a, b in regions:
            for rec in vf.fetch(chrom, a, b):
                if len(rec.ref) != 1 or len(rec.alts or ()) != 1 or len(rec.alts[0]) != 1:
                    continue
                s = rec.samples[smp]
                pw = s.get("PW") if has_pw else None
                if pw in ("0|1", "1|0"):
                    idx = (int(pw[0]), int(pw[2]))
                elif pw in (None, ".", "") and s.phased and tuple(s["GT"]) in ((0, 1), (1, 0)):
                    idx = tuple(s["GT"])
                    qc["het_gt_fallback"] += 1
                else:
                    continue
                alle = (rec.ref, rec.alts[0])
                h = out[chrom]
                if rec.pos in h.hap_alleles:
                    continue
                h.pos.append(rec.pos)
                h.hap_alleles[rec.pos] = (alle[idx[0]], alle[idx[1]])
                h.ref_alt[rec.pos] = alle
                qc["het_sites"] += 1
    for h in out.values():
        h.pos.sort()
    return out, qc


def _pair_models(pairs: pd.DataFrame, junctions: pd.DataFrame
                 ) -> tuple[list[PairModel], dict[tuple, list[int]]]:
    models, index = [], defaultdict(list)
    jg = {k: g for k, g in junctions.groupby("pair_id", sort=False)}
    for r in pairs.itertuples(index=False):
        g = jg.get(r.pair_id)
        if g is None:
            continue
        intr, spec = [], []
        for side in (1, 2):
            s = g[g["isoform"] == side]
            intr.append(frozenset(zip(s["start"].astype(int), s["end"].astype(int))))
            ss = s[s["specific"]]
            spec.append(frozenset(zip(ss["start"].astype(int), ss["end"].astype(int))))
        k = len(models)
        models.append(PairModel(r.pair_id, (int(r.span_start), int(r.span_end)),
                                (intr[0], intr[1]), (spec[0], spec[1])))
        for j in spec[0] | spec[1]:
            index[(r.chrom, j[0], j[1])].append(k)
    return models, index


def count_bam(bam_path: str, regions: list[tuple[str, int, int]], het: dict[str, HetSites],
              models: list[PairModel], index: dict[tuple, list[int]],
              min_anchor: int = MIN_ANCHOR, min_bq: int = MIN_BASE_QUALITY):
    """Walk the BAM over the regions once; returns (pair_counts, pair_sites, site_counts, qc).

    pair_counts[(pair_id, isoform, hap)] = fragments; pair_sites adds the het position (a
    fragment with two sites counts once per site there); site_counts[(chrom, pos, hap)]
    counts every passing fragment at each phased site, isoform-agnostic, for the
    concordance check against the phASER release.
    """
    import pysam
    pair_counts: Counter = Counter()
    pair_sites: Counter = Counter()
    site_counts: Counter = Counter()
    qc: Counter = Counter()

    def finalize(chrom: str, mates: list[ReadRec]) -> None:
        qc["fragments"] += 1
        status, junc, sites = merge_mates(mates)
        if status != "ok":
            qc[status] += 1
            return
        hap = fragment_hap(sites)
        if hap == "conflict":
            qc["hap_conflict"] += 1
        else:
            for p, h in sites.items():
                site_counts[(chrom, p, h)] += 1
        cand = {k for j in junc for k in index.get((chrom, j[0], j[1]), ())}
        if not cand:
            return
        qc["fragments_on_specific_junction"] += 1
        if hap == "conflict":
            return
        for k in cand:
            iso = assign_isoform(junc, models[k])
            if iso is None:
                continue
            if isinstance(iso, str):
                qc[f"iso_{iso}"] += 1
                continue
            pair_counts[(models[k].pair_id, iso, hap)] += 1
            for p in sites:
                pair_sites[(models[k].pair_id, iso, hap, p)] += 1

    with pysam.AlignmentFile(bam_path) as bam:
        prev: tuple[str, int, int] | None = None
        for chrom, a, b in regions:
            sites = het.get(chrom, HetSites())
            pending: dict[str, tuple[ReadRec, int]] = {}
            for read in bam.fetch(chrom, a, b):
                if (prev is not None and prev[0] == chrom and read.reference_start < a
                        and read.reference_start < prev[2]):
                    continue  # overlaps the previous region too; already counted there
                if read.flag & _BAD_FLAGS:
                    continue
                if read.has_tag("NH") and read.get_tag("NH") != 1:
                    qc["multimapped"] += 1
                    continue
                if read.is_duplicate:
                    qc["duplicate_flagged"] += 1
                rec = ReadRec(tuple(read_junctions(read.cigartuples, read.reference_start,
                                                   min_anchor)),
                              read_het_obs(read, sites, min_bq),
                              read.get_tag("vW") if read.has_tag("vW") else None)
                if (not read.is_paired or read.mate_is_unmapped
                        or read.next_reference_name != chrom):
                    finalize(chrom, [rec])
                    continue
                mate = pending.pop(read.query_name, None)
                if mate is not None:
                    finalize(chrom, [mate[0], rec])
                else:
                    pending[read.query_name] = (rec, read.next_reference_start)
            # A mate that never arrived lies outside this region. If it starts AFTER the
            # region, it has not been seen anywhere, so the fragment is counted here from
            # this mate alone; if it starts BEFORE, the fragment was already counted from
            # that mate in an earlier region, and counting it again would double it.
            for rec, mate_start in pending.values():
                if mate_start < a:
                    qc["mate_counted_elsewhere"] += 1
                    continue
                qc["mate_outside_region"] += 1
                finalize(chrom, [rec])
            prev = (chrom, a, b)
    return pair_counts, pair_sites, site_counts, qc


def _read_regions(bed: Path) -> list[tuple[str, int, int]]:
    out = []
    for line in bed.read_text().splitlines():
        c, a, b = line.split("\t")[:3]
        out.append((c, int(a), int(b)))
    return out


def run_count(region: str, sample: str | None, task: int | None,
              genes: set[str] | None, dest_override: Path | None,
              min_anchor: int, min_bq: int) -> Path:
    src = out_dir(region)
    samples = pd.read_csv(src / "samples.tsv", sep="\t")
    if sample is None:
        if task is None:
            raise SystemExit("pass --sample or --task (1-based row of samples.tsv)")
        if task > len(samples):
            # the array is sized from the manifest before `targets` has run, so it can
            # overshoot the samples that survived; those tasks have nothing to do
            print(f"  task {task} > {len(samples)} samples in {region}; nothing to count")
            return src
        sample = str(samples.iloc[task - 1]["sample_id"])
    row = samples[samples["sample_id"] == sample]
    if row.empty:
        raise SystemExit(f"{sample} is not in {src / 'samples.tsv'}")
    row = row.iloc[0]
    pairs = pd.read_parquet(src / "pairs.parquet")
    pairs = pairs[pairs["in_gtf"]]
    if genes:
        pairs = pairs[pairs["gene_id"].astype(str).str.split(".").str[0].isin(genes)]
        regions = []
        for chrom, g in pairs.groupby("chrom"):
            for a, b in _merge(list(zip(g["span_start"].astype(int),
                                        g["span_end"].astype(int)))):
                regions.append((chrom, a - 1, b))
    else:
        regions = _read_regions(src / "regions.bed")
    junctions = pd.read_parquet(src / "junctions.parquet")
    junctions = junctions[junctions["pair_id"].isin(pairs["pair_id"])]
    dest = ensure_dir(dest_override or (src / "counts"))

    t0 = time.time()
    models, index = _pair_models(pairs, junctions)
    het, hqc = load_het_sites(str(row["phaser_vcf"]), regions)
    pc, ps, sc, qc = count_bam(str(row["bam"]), regions, het, models, index,
                               min_anchor, min_bq)
    qc.update(hqc)

    def _write(name: str, cols: list[str], counter: Counter) -> None:
        df = pd.DataFrame([(*k, v) for k, v in counter.items()], columns=[*cols, "n_frag"])
        df.insert(0, "donor_id", row["donor_id"])
        df.insert(0, "sample_id", sample)
        df.to_csv(dest / f"{sample}.{name}.tsv.gz", sep="\t", index=False)

    _write("pair_counts", ["pair_id", "isoform", "hap"], pc)
    _write("pair_sites", ["pair_id", "isoform", "hap", "het_pos"], ps)
    _write("site_counts", ["chrom", "pos", "hap"], sc)
    info = {"sample_id": sample, "donor_id": row["donor_id"], "bam": row["bam"],
            "n_regions": len(regions), "n_pairs": len(models),
            "seconds": round(time.time() - t0, 1), **{k: int(v) for k, v in qc.items()}}
    (dest / f"{sample}.qc.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))
    return dest


# --------------------------------------------------------------------------- #
# Stage: concordance -- recount vs the phASER release at the same sites
# --------------------------------------------------------------------------- #
def run_concordance(region: str, sample: str, counts_dir: Path | None,
                    min_depth: int = 20) -> dict:
    cfg = quest_config(region)
    src = out_dir(region)
    cdir = counts_dir or (src / "counts")
    sc = pd.read_csv(cdir / f"{sample}.site_counts.tsv.gz", sep="\t")
    samples = pd.read_csv(src / "samples.tsv", sep="\t")
    vcf = samples.loc[samples["sample_id"] == sample, "phaser_vcf"].iloc[0]
    regions = sorted({(c, int(p) - 1, int(p)) for c, p in zip(sc["chrom"], sc["pos"])})
    het, _ = load_het_sites(vcf, regions)
    wide = sc.pivot_table(index=["chrom", "pos"], columns="hap", values="n_frag",
                          aggfunc="sum", fill_value=0).reset_index()
    for h in (1, 2):
        if h not in wide.columns:
            wide[h] = 0
    ref_n, alt_n = [], []
    for c, p, n1, n2 in zip(wide["chrom"], wide["pos"], wide[1], wide[2]):
        h1 = het[c].hap_alleles[p][0]
        ref = het[c].ref_alt[p][0]
        ref_n.append(n1 if h1 == ref else n2)
        alt_n.append(n2 if h1 == ref else n1)
    wide["ref_recount"], wide["alt_recount"] = ref_n, alt_n
    ac = pd.read_csv(Path(cfg["per_sample_dir"]) / f"{sample}.allelic_counts.txt", sep="\t",
                     usecols=["contig", "position", "refCount", "altCount", "totalCount"])
    m = wide.merge(ac, left_on=["chrom", "pos"], right_on=["contig", "position"])
    tot_re = m["ref_recount"] + m["alt_recount"]
    deep = m[(tot_re >= min_depth) & (m["totalCount"] >= min_depth)]
    r = float(np.corrcoef(deep["ref_recount"] / (deep["ref_recount"] + deep["alt_recount"]),
                          deep["refCount"] / deep["totalCount"])[0, 1]) if len(deep) > 2 else float("nan")
    out = {"sample_id": sample, "sites_recount": int(len(wide)),
           "sites_in_release": int(len(m)), "sites_ge_min_depth": int(len(deep)),
           "min_depth": min_depth, "pearson_ref_fraction": r,
           "depth_ratio_median": float((tot_re / m["totalCount"].clip(lower=1)).median())
           if len(m) else float("nan"),
           "ref_fraction_recount": float(m["ref_recount"].sum() / max(tot_re.sum(), 1)),
           "ref_fraction_release": float(m["refCount"].sum() / max(m["totalCount"].sum(), 1))}
    (cdir / f"{sample}.concordance.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    return out


# --------------------------------------------------------------------------- #
# Stage: screen
# --------------------------------------------------------------------------- #
def gate(counts: pd.DataFrame, min_donors: int, min_reads: int) -> pd.DataFrame:
    """The exonic screen's gate, applied to junction-informative fragments.

    Informative = a specific junction AND a phased haplotype (hap 1/2). Identical rules
    to `ase_switch_direction.run_screen`: >= `min_donors` donors with >= `min_reads`
    informative fragments, and informative fragments on BOTH isoforms.
    """
    inf = counts[counts["hap"].isin([1, 2])]
    per_donor = inf.groupby(["pair_id", "donor_id"], as_index=False)["n_frag"].sum()
    out = (inf.groupby("pair_id", as_index=False)
           .agg(n_donors=("donor_id", "nunique"), total_frag=("n_frag", "sum"),
                n_isoforms_informative=("isoform", "nunique")))
    ok = (per_donor[per_donor["n_frag"] >= min_reads].groupby("pair_id")["donor_id"]
          .nunique().rename("n_donors_ge_min_reads"))
    med = per_donor.groupby("pair_id")["n_frag"].median().rename("median_depth")
    both = (inf.groupby(["pair_id", "donor_id"])["isoform"].nunique()
            .eq(2).groupby(level="pair_id").sum().rename("n_donors_both_isoforms"))
    out = out.join(ok, on="pair_id").join(med, on="pair_id").join(both, on="pair_id")
    out["n_donors_ge_min_reads"] = out["n_donors_ge_min_reads"].fillna(0).astype(int)
    out["testable"] = ((out["n_donors_ge_min_reads"] >= min_donors)
                       & (out["n_isoforms_informative"] >= 2))
    return out


def run_screen(region: str, min_donors: int, min_reads: int, allow_partial: bool) -> Path:
    dest = out_dir(region)
    cdir = dest / "counts"
    samples = pd.read_csv(dest / "samples.tsv", sep="\t")
    have = {p.name.split(".")[0] for p in cdir.glob("*.qc.json")}
    missing = sorted(set(samples["sample_id"]) - have)
    if missing and not allow_partial:
        raise SystemExit(f"{len(missing)} of {len(samples)} count shards missing "
                         f"(first: {missing[:5]}); rerun those tasks or pass --allow-partial")
    shards = [pd.read_csv(cdir / f"{s}.pair_counts.tsv.gz", sep="\t")
              for s in sorted(have)]
    counts = pd.concat(shards, ignore_index=True)
    counts.to_parquet(dest / "junction_allelic_counts.parquet", index=False)
    qc = pd.DataFrame([json.loads((cdir / f"{s}.qc.json").read_text()) for s in sorted(have)])
    qc.to_parquet(dest / "count_qc.parquet", index=False)

    pairs = pd.read_parquet(dest / "pairs.parquet")
    g = gate(counts, min_donors, min_reads)
    per_pair = pairs.merge(g, on="pair_id", how="left")
    for c in ("n_donors", "n_donors_ge_min_reads", "n_isoforms_informative",
              "n_donors_both_isoforms", "total_frag"):
        per_pair[c] = per_pair[c].fillna(0).astype(int)
    per_pair["testable"] = per_pair["testable"].fillna(False).astype(bool)

    # descriptive only, never a gate criterion: informative donors that are heterozygous at
    # the lead, i.e. the donors stage 4 could actually orient
    lg = pd.read_parquet(dest / "lead_genotypes.parquet")
    het_lead = lg[lg["gt"].isin(["0|1", "1|0"])][["variant_id", "donor_id"]]
    inf = counts[counts["hap"].isin([1, 2])]
    pdn = inf.groupby(["pair_id", "donor_id"], as_index=False)["n_frag"].sum()
    pdn = pdn[pdn["n_frag"] >= min_reads].merge(
        pairs[["pair_id", "variant_id_all"]], on="pair_id")
    pdn = pdn.merge(het_lead, left_on=["variant_id_all", "donor_id"],
                    right_on=["variant_id", "donor_id"])
    per_pair = per_pair.merge(pdn.groupby("pair_id")["donor_id"].nunique()
                              .rename("n_donors_ge_min_reads_het_lead"),
                              on="pair_id", how="left")
    per_pair["n_donors_ge_min_reads_het_lead"] = (
        per_pair["n_donors_ge_min_reads_het_lead"].fillna(0).astype(int))
    per_pair = per_pair.sort_values(["gate_family", "n_donors_ge_min_reads"],
                                    ascending=False)
    per_pair.to_parquet(dest / "pair_feasibility.parquet", index=False)

    def fam(mask: pd.Series) -> dict:
        s = per_pair[mask & per_pair["in_gtf"]]
        return {"n_pairs": int(len(s)), "n_genes": int(s["gene_id"].nunique()),
                "n_two_sided": int(s["two_sided"].fillna(False).sum()),
                "n_any_informative": int((s["n_donors"] > 0).sum()),
                "n_testable": int(s["testable"].sum()),
                "n_genes_testable": int(s.loc[s["testable"], "gene_id"].nunique())}

    fams = {"gate_family": fam(per_pair["gate_family"].fillna(False)),
            "coloc_nominated": fam(per_pair["coloc_brainseq"] | per_pair["coloc_audit"]),
            "all_pairs": fam(pd.Series(True, index=per_pair.index))}
    tot = qc.sum(numeric_only=True)
    passes = fams["gate_family"]["n_testable"] >= MIN_PAIRS_TO_PROCEED
    summary = {
        "region": region, "min_donors_informative": min_donors,
        "min_reads_per_donor": min_reads, "min_pairs_to_proceed": MIN_PAIRS_TO_PROCEED,
        "n_samples": int(len(qc)), "n_samples_missing": len(missing),
        "families": fams, "gate_passed": bool(passes),
        "hap_conflict_rate": float(tot.get("hap_conflict", 0) / max(tot.get("fragments", 1), 1)),
        "wasp_fail_rate": float(tot.get("wasp_fail", 0) / max(tot.get("fragments", 1), 1)),
        "wasp_missing": int(tot.get("wasp_missing", 0)),
        "het_gt_fallback": int(tot.get("het_gt_fallback", 0)),
    }
    (dest / "screen_summary.json").write_text(json.dumps(summary, indent=2))
    _write_report(dest, region, per_pair, summary)
    print(json.dumps(summary, indent=2))
    return dest


def _write_report(dest: Path, region: str, per_pair: pd.DataFrame, summary: dict) -> None:
    L: list[str] = []
    A = L.append
    gf = summary["families"]["gate_family"]
    A(f"# Allele-aware junction screen — {region}")
    A("")
    A("Can reads that cross an isoform-specific splice junction **and** carry a phased "
      "heterozygous site support a within-donor allelic test of an isoform switch? This "
      "is the two-sided rescue of the exonic phASER screen (`ase_switch_direction.py`), "
      "counted from the WASP-tagged STAR BAMs on Quest.")
    A("")
    A("## Pre-registered gate (unchanged from the exonic screen)")
    A("")
    A(f"- per pair: at least **{summary['min_donors_informative']}** donors, each with at "
      f"least **{summary['min_reads_per_donor']}** informative fragments, informative on "
      f"**both** isoforms")
    A(f"- to proceed: at least **{summary['min_pairs_to_proceed']}** pairs pass in the gate "
      f"family (pairs in genes with an all_samples switch-QTL, q < {QVAL_GATE_FAMILY})")
    A("")
    A("## Result")
    A("")
    A("| family | pairs | genes | two-sided | any informative | **pass** | genes passing |")
    A("|---|---|---|---|---|---|---|")
    for name, f in summary["families"].items():
        A(f"| {name} | {f['n_pairs']:,} | {f['n_genes']:,} | {f['n_two_sided']:,} | "
          f"{f['n_any_informative']:,} | **{f['n_testable']:,}** | {f['n_genes_testable']:,} |")
    A("")
    if summary["gate_passed"]:
        A(f"**The gate passes: {gf['n_testable']} gate-family pairs clear it.** Stage 4 "
          "(allele orientation at the lead, then the beta-binomial contrast) may be built; "
          "its design is in the module docstring.")
    else:
        A(f"**The gate fails: {gf['n_testable']} gate-family pairs clear it, fewer than "
          f"{summary['min_pairs_to_proceed']}.** As pre-specified, the allele-specific arm "
          "is dropped and reported as a negative result.")
    A("")
    A("## Quality")
    A("")
    A(f"- samples counted: {summary['n_samples']} (missing: {summary['n_samples_missing']})")
    A(f"- haplotype-conflict rate among fragments: {summary['hap_conflict_rate']:.4f}")
    A(f"- WASP-fail rate among fragments: {summary['wasp_fail_rate']:.4f}; fragments with a "
      f"het site but no WASP tag: {summary['wasp_missing']}")
    A(f"- het sites phased from GT because PW was absent: {summary['het_gt_fallback']}")
    A("")
    A("## Caveats")
    A("")
    A("- Cell composition cancels within a donor only if allelic effects do not differ by "
      "cell type.")
    A("- Some reference mapping bias survives WASP.")
    A("- The per-sample WASP VCFs carry `0|0`/`1|1` rows as well as hets (the `-g het` "
      "filter was applied before sample subsetting); heterozygous sites here are taken "
      "from phASER's `PW` genotype, so this does not reach the counts.")
    A("")
    A("## Top gate-family pairs")
    A("")
    A("| gene | T1 | T2 | donors (>= min) | both-isoform donors | het at lead | "
      "isoforms | pass |")
    A("|---|---|---|---|---|---|---|---|")
    top = per_pair[per_pair["gate_family"].fillna(False)].head(25)
    for r in top.itertuples(index=False):
        A(f"| {r.gene_id} | {r.transcript_id_1} | {r.transcript_id_2} | "
          f"{r.n_donors_ge_min_reads} | {r.n_donors_both_isoforms} | "
          f"{r.n_donors_ge_min_reads_het_lead} | {r.n_isoforms_informative} | "
          f"{'**yes**' if r.testable else 'no'} |")
    A("")
    (dest / "ASE_JUNCTION_SCREEN.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True,
                    choices=("targets", "count", "concordance", "screen"))
    ap.add_argument("--region", choices=REGIONS, default="dlpfc")
    ap.add_argument("--sample", default=None, help="count/concordance: R##### sample id")
    ap.add_argument("--task", type=int, default=None,
                    help="count: 1-based row of samples.tsv (SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--genes", default=None,
                    help="count: comma-separated bare ENSG ids (smoke tests)")
    ap.add_argument("--out", type=Path, default=None,
                    help="count/concordance: shard directory (default <out>/counts)")
    ap.add_argument("--min-anchor", type=int, default=MIN_ANCHOR)
    ap.add_argument("--min-base-quality", type=int, default=MIN_BASE_QUALITY)
    ap.add_argument("--min-donors", type=int, default=MIN_DONORS_INFORMATIVE)
    ap.add_argument("--min-reads", type=int, default=MIN_READS_PER_DONOR)
    ap.add_argument("--allow-partial", action="store_true",
                    help="screen: proceed with missing shards (reported in the summary)")
    args = ap.parse_args(argv)

    if args.stage == "targets":
        run_targets(args.region)
    elif args.stage == "count":
        genes = {g.strip() for g in args.genes.split(",")} if args.genes else None
        run_count(args.region, args.sample, args.task, genes, args.out,
                  args.min_anchor, args.min_base_quality)
    elif args.stage == "concordance":
        if not args.sample:
            raise SystemExit("--stage concordance needs --sample")
        run_concordance(args.region, args.sample, args.out)
    else:
        run_screen(args.region, args.min_donors, args.min_reads, args.allow_partial)


if __name__ == "__main__":
    main()
