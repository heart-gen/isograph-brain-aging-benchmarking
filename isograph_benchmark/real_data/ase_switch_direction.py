"""Allele-specific direction of an isoform switch, from the BrainSEQ phASER release.

WHAT THIS IS FOR
----------------
Colocalization says a risk variant and a splicing signal share a causal variant. It does
not say the risk ALLELE drives the switch. The within-donor allelic test does: in a
heterozygote, the two haplotypes sit in the same nucleus, same cell-type mixture, same
environment, so a cis effect shows up as a difference between them.

WHY THIS MODULE LEADS WITH A FEASIBILITY SCREEN
-----------------------------------------------
The obvious implementation is wrong, and expensively so. `brainseq.gene_ae.tsv.gz` is
keyed on ENST and looks like allele-specific isoform quantification. It is not.
`phaser_gene_ae` aggregates haplotypic reads over arbitrary BED features; the pipeline
was handed transcript intervals, so an ENST there names a genomic INTERVAL. Reads and
variants in exons shared by overlapping transcripts are counted into several
"transcript" features at once. Nothing in that file distinguishes isoforms by itself.

An allelic read is only isoform-informative when it carries BOTH a heterozygous site
(which says which haplotype) AND sequence unique to one isoform of the pair (which says
which isoform). So the first question is not "what is the effect" but "does that read
exist, in how many donors, at what depth". The QC tables say to expect very little:
median `gene_ae` depth is 1-2 reads per sample and median haplotype depth is 2.

This stage answers that question and nothing else. A negative here is a real result and
is cheaper to publish than a strained positive.

  --stage screen   Per switch pair: find isoform-DISCRIMINATING exonic sequence, intersect
                   it with the heterozygous sites actually observed in the phASER allelic
                   counts, and report informative donors and depth. Gates everything else.

Stages 2 and 3 (allele orientation, and the transcript x haplotype beta-binomial
contrast) are deliberately NOT implemented until the screen says there is signal to
model. Their design is in the plan; building them first would be writing a test for data
that may not exist.

WHERE THE JUNCTION-LEVEL RESCUE HAS TO RUN -- NOT ON BRIDGES-2
--------------------------------------------------------------
The screen's result (see the report) is that depth is fine and the binding constraint is
structural: for most switch pairs only ONE isoform carries unique exonic sequence with a
heterozygous site, so there is no second side to contrast. The fix is to count reads that
cross an isoform-DISCRIMINATING JUNCTION while carrying a phased heterozygous site --
two-sided by construction, because each isoform splices a different junction.

That cannot be done from anything in this repo or on this cluster:

  * It needs read-level alignments. Every phASER product shipped here is already
    aggregated to sites, haplotype blocks or feature intervals -- `allelic_counts` is
    per-site, `haplotypic_counts` is per-phase-block. None of them retains which read
    crossed which junction, and `junction_tables/*_SJ.out.tab` are STAR junction counts
    with no allele information at all.
  * It needs the WASP-tagged, ASE-grade BAMs specifically, not any alignment. Reference
    mapping bias at heterozygous sites systematically favours the reference allele --
    alt-carrying reads mismap -- and that bias is the dominant confound in exactly the
    within-donor allelic contrast this analysis makes. `ASE_GENERATION.md` is explicit
    that a standard alignment without WASP-aware handling is not ASE-final. The existing
    phASER numbers already have that treatment; new counting from scratch must inherit it.
  * Those BAMs are NOT on Bridges-2. The validated manifest's `bam_path` entries are
    relative (`../../_m/R#####_Aligned.sortedByCoord.out.bam`) to the original processing
    tree on Northwestern's Quest under `/gpfs/projects/b1042/HEART-GeN-Lab/ase-processing/`,
    which is where the `source_file` column also points. A search of both PSC project
    shares turns up no BAM or CRAM at all.

So the junction-level arm is a **Quest (b1042) job**, where the alignments and the storage
live. CRAM would be fine there in principle -- samtools, pysam and GATK all read CRAM, and
the GRCh38 reference it needs is available -- but the format is not the constraint. The
constraint is that the WASP tags and the reads themselves only exist on Quest.

WHAT STAGE 2 WILL HAVE TO DO, RECORDED HERE SO IT IS NOT FORGOTTEN
------------------------------------------------------------------
`gw_phased = 1` makes haplotype A internally consistent WITHIN a donor. It does not make
haplotype A the same biological allele ACROSS donors. Pooling ~450 donors without first
orienting each one to the risk allele would average a real cis effect toward zero. That
orientation, not the model, is the hard part.
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel, stage_out

REGIONS: tuple[str, ...] = ("caudate", "dlpfc", "hippocampus")
_ASE_ROOT = rel("inputs", "raw", "brainseq", "ase")
_GTF_CACHE = stage_out("tmp", "gencode.v47.primary_assembly.annotation.gtf_cache.parquet")

# Bucket width for the coarse pre-filter over the allelic-count stream. The files are
# 2.5-3 GB gzipped per region and carry no index, so they are streamed once; bucketing
# the discriminating intervals into a hash lets awk reject the vast majority of rows
# without an interval search, and the exact assignment is done afterwards in pandas.
_BUCKET = 10_000

# Pre-registered before the screen was run. A pair is only worth modelling if a real
# number of donors carry a read that is simultaneously haplotype- and isoform-informative.
MIN_DONORS_INFORMATIVE = 30
MIN_READS_PER_DONOR = 5


def out_dir(region: str) -> Path:
    return ensure_dir(stage_out("mechanism", "ase_switch_direction", region))


def _store(region: str) -> Path:
    from isograph_benchmark.paths import region_store
    return region_store("brainseq", region)


# --------------------------------------------------------------------------- #
# Isoform-discriminating sequence
# --------------------------------------------------------------------------- #
def _merge(iv: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(a, b) for a, b in out]


def _subtract(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Interval difference a \\ b, both assumed merged and sorted."""
    out = []
    for s, e in a:
        cur = [(s, e)]
        for bs, be in b:
            nxt = []
            for cs, ce in cur:
                if be < cs or bs > ce:
                    nxt.append((cs, ce))
                    continue
                if bs > cs:
                    nxt.append((cs, bs - 1))
                if be < ce:
                    nxt.append((be + 1, ce))
            cur = nxt
            if not cur:
                break
        out.extend(cur)
    return _merge(out)


def discriminating_segments(region: str, genes: set[str] | None = None
                            ) -> pd.DataFrame:
    """Exonic sequence unique to ONE member of each switch pair.

    A heterozygous site inside such a segment is isoform-informative: a read covering it
    can only have come from that isoform. Sequence shared by both members carries the
    haplotype but not the isoform, which is exactly the ambiguity that makes `gene_ae`
    unusable for this question.
    """
    sp = _store(region) / "isograph_vae" / "module_interpret" / "structure_switch_pairs.parquet"
    if not sp.exists():
        raise SystemExit(f"missing switch pairs: {sp}")
    pairs = pd.read_parquet(sp)
    pairs["gene_bare"] = pairs["gene_id"].astype(str).str.split(".").str[0]
    if genes:
        pairs = pairs[pairs["gene_bare"].isin(genes)]
    if pairs.empty:
        raise SystemExit("no switch pairs for the requested genes")

    gtf = pd.read_parquet(_GTF_CACHE, columns=["transcript_id", "gene_id", "chrom",
                                               "strand", "feature", "start", "end"])
    gtf = gtf[gtf["feature"] == "exon"]
    ex = {t: g for t, g in gtf.groupby("transcript_id", sort=False)}

    rows = []
    for r in pairs.itertuples(index=False):
        t1, t2 = str(r.transcript_id_1), str(r.transcript_id_2)
        if t1 not in ex or t2 not in ex:
            continue
        e1, e2 = ex[t1], ex[t2]
        chrom = str(e1["chrom"].iloc[0])
        i1 = _merge(list(zip(e1["start"].astype(int), e1["end"].astype(int))))
        i2 = _merge(list(zip(e2["start"].astype(int), e2["end"].astype(int))))
        for tx, uniq in ((t1, _subtract(i1, i2)), (t2, _subtract(i2, i1))):
            for s, e in uniq:
                rows.append({"gene": r.gene_bare, "transcript_id_1": t1,
                             "transcript_id_2": t2, "unique_to": tx,
                             "chrom": chrom, "start": s, "end": e,
                             "length": e - s + 1})
    seg = pd.DataFrame(rows)
    if seg.empty:
        raise SystemExit("no isoform-discriminating exonic sequence found")
    return seg


# --------------------------------------------------------------------------- #
# Heterozygous sites actually observed
# --------------------------------------------------------------------------- #
def _stream_allelic_counts(region: str, seg: pd.DataFrame) -> pd.DataFrame:
    """One pass over the region's allelic counts, keeping sites in discriminating sequence.

    phASER's `allelic_counts` is a per-donor, per-site table with `refCount`/`altCount`,
    so a site with both alleles observed in a donor is a heterozygous site that carries
    allelic information. Only sites inside discriminating sequence are kept.
    """
    f = _ASE_ROOT / region / "brainseq.allelic_counts.tsv.gz"
    if not f.exists():
        raise SystemExit(f"missing phASER allelic counts: {f}")

    buckets = set()
    for r in seg.itertuples(index=False):
        for b in range(r.start // _BUCKET, r.end // _BUCKET + 1):
            buckets.add(f"{r.chrom}:{b}")
    keyfile = out_dir(region) / "_buckets.txt"
    keyfile.write_text("\n".join(sorted(buckets)) + "\n")
    print(f"  coarse filter: {len(buckets):,} {_BUCKET//1000}kb buckets over "
          f"{len(seg):,} discriminating segments "
          f"({seg['length'].sum()/1e3:,.0f} kb of unique exonic sequence)")

    awk = (
        'BEGIN{FS=OFS="\\t"; while((getline k < KF) > 0) keep[k]=1}'
        'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print "sample_id","donor_id","contig",'
        '"position","refCount","altCount","totalCount"; next}'
        '{b=int($c["position"]/BW); if(($c["contig"]":"b) in keep)'
        ' print $c["sample_id"],$c["donor_id"],$c["contig"],$c["position"],'
        '$c["refCount"],$c["altCount"],$c["totalCount"]}'
    )
    cmd = f'zcat {f} | awk -v KF={keyfile} -v BW={_BUCKET} \'{awk}\''
    print("  streaming allelic counts (single pass over a multi-GB gzip)...")
    p = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"allelic-count stream failed: {p.stderr[-2000:]}")
    from io import StringIO
    d = pd.read_csv(StringIO(p.stdout), sep="\t")
    print(f"  kept {len(d):,} (donor, site) rows in candidate buckets")
    return d


def _assign_exact(counts: pd.DataFrame, seg: pd.DataFrame) -> pd.DataFrame:
    """Exact interval assignment of sites to discriminating segments, per chromosome."""
    out = []
    for chrom, s in seg.groupby("chrom"):
        c = counts[counts["contig"] == chrom]
        if c.empty:
            continue
        s = s.sort_values("start")
        starts = s["start"].to_numpy()
        ends = s["end"].to_numpy()
        pos = c["position"].to_numpy()
        idx = np.searchsorted(starts, pos, side="right") - 1
        # searchsorted finds the last segment starting at or before the site; segments are
        # merged per (pair, transcript) but may still overlap ACROSS pairs, so this keeps
        # one assignment per site and the per-pair totals are computed by re-joining.
        ok = (idx >= 0) & (pos <= np.where(idx >= 0, ends[idx], -1))
        if not ok.any():
            continue
        hit = c[ok].copy()
        take = s.iloc[idx[ok]]
        for col in ("gene", "transcript_id_1", "transcript_id_2", "unique_to"):
            hit[col] = take[col].to_numpy()
        out.append(hit)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# --------------------------------------------------------------------------- #
# Stage: screen
# --------------------------------------------------------------------------- #
def run_screen(region: str, genes: set[str] | None = None,
               min_donors: int = MIN_DONORS_INFORMATIVE,
               min_reads: int = MIN_READS_PER_DONOR, reuse: bool = False) -> Path:
    dest = out_dir(region)
    cached = dest / "informative_sites.parquet"
    if reuse and cached.exists():
        # The stream is a single pass over a 2.5 GB gzip and takes ~15 minutes. Rewording
        # the report should not cost that, so the assigned sites are reusable. The
        # thresholds are applied downstream of this cache, so `--min-donors` /
        # `--min-reads` can be swept without re-reading the file.
        print(f"  reusing cached informative sites: {cached}")
        seg = pd.read_parquet(dest / "discriminating_segments.parquet")
        hits = pd.read_parquet(cached)
    else:
        seg = discriminating_segments(region, genes)
        seg.to_parquet(dest / "discriminating_segments.parquet", index=False)
        print(f"  {seg['gene'].nunique():,} genes, "
              f"{seg.groupby(['transcript_id_1','transcript_id_2']).ngroups:,} switch "
              f"pairs carry unique exonic sequence")
        counts = _stream_allelic_counts(region, seg)
        hits = _assign_exact(counts, seg)
        if hits.empty:
            raise SystemExit("no allelic sites fell inside discriminating sequence")
        # A site is allele-informative in a donor only when BOTH alleles are seen.
        hits["het_informative"] = (hits["refCount"] > 0) & (hits["altCount"] > 0)
        hits.to_parquet(cached, index=False)

    per_pair = (hits[hits["het_informative"]]
                .groupby(["gene", "transcript_id_1", "transcript_id_2"], as_index=False)
                .agg(n_donors=("donor_id", "nunique"),
                     n_sites=("position", "nunique"),
                     total_reads=("totalCount", "sum"),
                     median_depth=("totalCount", "median")))
    # Donors clearing the per-donor read floor, which is the number that actually powers
    # a within-donor test.
    per_donor = (hits[hits["het_informative"]]
                 .groupby(["gene", "transcript_id_1", "transcript_id_2", "donor_id"],
                          as_index=False)["totalCount"].sum())
    ok = (per_donor[per_donor["totalCount"] >= min_reads]
          .groupby(["gene", "transcript_id_1", "transcript_id_2"], as_index=False)
          .agg(n_donors_ge_min_reads=("donor_id", "nunique")))
    per_pair = per_pair.merge(ok, on=["gene", "transcript_id_1", "transcript_id_2"],
                              how="left")
    per_pair["n_donors_ge_min_reads"] = per_pair["n_donors_ge_min_reads"].fillna(0).astype(int)

    # BOTH isoforms must be informative: the statistic is a CONTRAST between them, so a
    # pair informative on only one side cannot be tested however deep that side is.
    both = (hits[hits["het_informative"]]
            .groupby(["gene", "transcript_id_1", "transcript_id_2"])["unique_to"]
            .nunique().rename("n_isoforms_informative").reset_index())
    per_pair = per_pair.merge(both, on=["gene", "transcript_id_1", "transcript_id_2"],
                              how="left")
    per_pair["testable"] = ((per_pair["n_donors_ge_min_reads"] >= min_donors)
                            & (per_pair["n_isoforms_informative"] >= 2))
    per_pair = per_pair.sort_values("n_donors_ge_min_reads", ascending=False)
    per_pair.to_parquet(dest / "pair_feasibility.parquet", index=False)

    enough_donors = per_pair["n_donors_ge_min_reads"] >= min_donors
    both_iso = per_pair["n_isoforms_informative"] >= 2
    summary = {
        "region": region,
        "n_fail_only_both_isoform_rule": int((enough_donors & ~both_iso).sum()),
        "n_fail_only_donor_count": int((~enough_donors & both_iso).sum()),
        "median_depth_of_donor_rich_failures": float(
            per_pair.loc[enough_donors & ~both_iso, "median_depth"].median())
            if int((enough_donors & ~both_iso).sum()) else float("nan"),
        "n_pairs_with_unique_sequence": int(
            seg.groupby(["transcript_id_1", "transcript_id_2"]).ngroups),
        "n_pairs_with_any_informative_site": int(len(per_pair)),
        "n_pairs_testable": int(per_pair["testable"].sum()),
        "n_genes_testable": int(per_pair.loc[per_pair["testable"], "gene"].nunique()),
        "min_donors_informative": min_donors,
        "min_reads_per_donor": min_reads,
        "median_donors_per_pair": float(per_pair["n_donors_ge_min_reads"].median()),
        "max_donors_per_pair": int(per_pair["n_donors_ge_min_reads"].max()),
    }
    (dest / "screen_summary.json").write_text(json.dumps(summary, indent=2))
    _write_report(dest, region, seg, per_pair, summary)
    print(f"\n  TESTABLE PAIRS: {summary['n_pairs_testable']:,} of "
          f"{summary['n_pairs_with_unique_sequence']:,} "
          f"(>= {min_donors} donors with >= {min_reads} informative reads on BOTH isoforms)")
    print(f"  wrote {dest}")
    return dest


def _write_report(dest: Path, region: str, seg: pd.DataFrame, per_pair: pd.DataFrame,
                  summary: dict) -> None:
    L: list[str] = []
    A = L.append
    A(f"# phASER feasibility screen — {region}")
    A("")
    A("Can the existing BrainSEQ phASER release support an allele-specific test of an "
      "isoform switch? A read is only isoform-informative if it carries a heterozygous "
      "site **and** sequence unique to one isoform of the pair. `gene_ae` cannot answer "
      "this: its ENST labels name genomic intervals, and reads in shared exons are "
      "counted into several transcript features at once.")
    A("")
    A("## Pre-registered thresholds")
    A("")
    A(f"- at least **{summary['min_donors_informative']}** donors")
    A(f"- each with at least **{summary['min_reads_per_donor']}** informative reads")
    A("- informative on **both** isoforms of the pair (the statistic is a contrast)")
    A("")
    A("## Result")
    A("")
    A(f"- switch pairs with isoform-unique exonic sequence: "
      f"**{summary['n_pairs_with_unique_sequence']:,}**")
    A(f"- pairs with any informative site: **{summary['n_pairs_with_any_informative_site']:,}**")
    A(f"- **pairs passing the screen: {summary['n_pairs_testable']:,}** "
      f"({summary['n_genes_testable']:,} genes)")
    A(f"- median donors per pair: {summary['median_donors_per_pair']:.0f}; "
      f"best pair: {summary['max_donors_per_pair']:,}")
    A("")
    A("## Why pairs fail — and it is not what the QC tables suggested")
    A("")
    A(f"- pairs failing **only** the both-isoforms rule, with donors and depth to spare: "
      f"**{summary['n_fail_only_both_isoform_rule']}**")
    A(f"- pairs failing **only** on donor count: {summary['n_fail_only_donor_count']}")
    A(f"- median depth among those donor-rich failures: "
      f"**{summary['median_depth_of_donor_rich_failures']:.0f} reads**")
    A("")
    A("The phASER QC tables report a median `gene_ae` depth of 1-2 reads, which reads "
      "like a depth problem. It is not. That figure is the whole-transcript-interval "
      "aggregation; once attention is restricted to isoform-DISCRIMINATING exonic "
      "sequence, the donors that have any signal have plenty of it. **The binding "
      "constraint is structural: for most switch pairs only ONE of the two isoforms "
      "carries unique exonic sequence containing a heterozygous site.** A pair where one "
      "isoform's sequence is a subset of the other's — a truncation, an alternative last "
      "exon — has no unique sequence on that side at all, so there is nothing to contrast "
      "against however deeply it is sequenced.")
    A("")
    if summary["n_pairs_testable"] == 0:
        A("**No pair clears the screen**, and the reason above says how it might be "
          "rescued rather than abandoned.")
    else:
        A(f"**{summary['n_pairs_testable']} pairs clear the screen** and could proceed to "
          "allele orientation. That is too few to carry a claim on its own.")
    A("")
    A("### The rescue this points at")
    A("")
    A("Exonic-unique segments are a ONE-SIDED discriminator. A discriminating **junction** "
      "is two-sided by construction: isoform A splices one junction where isoform B "
      "splices another, so a junction-spanning read is informative for whichever isoform "
      "it came from. Counting reads that simultaneously cross a discriminating junction "
      "and carry a phased heterozygous site would test many of the 25 pairs that fail "
      "here for want of a second side. That needs allele-aware counting from the BAMs; "
      "`junction_tables/*_SJ.out.tab` are STAR junction counts and are **not** "
      "allele-aware, so they cannot substitute.")
    A("")
    A("**This cannot run on Bridges-2.** It needs read-level alignments, and every "
      "phASER product shipped here is already aggregated to sites, haplotype blocks or "
      "feature intervals — none retains which read crossed which junction. It also needs "
      "the **WASP-tagged, ASE-grade BAMs** specifically: reference mapping bias at "
      "heterozygous sites favours the reference allele, which is the dominant confound in "
      "exactly this within-donor contrast. Those BAMs live on Northwestern's **Quest** "
      "under `/gpfs/projects/b1042/HEART-GeN-Lab/ase-processing/` — the manifest's "
      "`bam_path` and `source_file` both point there, and neither PSC project share holds "
      "any BAM or CRAM. CRAM would be workable on Quest (samtools/pysam/GATK all read it, "
      "and the GRCh38 reference is available), but format is not the constraint: the WASP "
      "tags and the reads only exist there. **Run the junction-level arm on Quest (b1042), "
      "where the alignments and the storage are.**")
    A("")
    A("If stage 2 is attempted, note that orientation — not the model — is the hard part: "
      "`gw_phased = 1` makes haplotype A consistent WITHIN a donor but not the same "
      "biological allele ACROSS donors, so each donor must be oriented to the risk allele "
      "before pooling or a real cis effect averages toward zero.")
    A("")
    A("## Top pairs by informative donors")
    A("")
    A("| gene | T1 | T2 | donors (>= min reads) | sites | isoforms informative | testable |")
    A("|---|---|---|---|---|---|---|")
    for r in per_pair.head(25).itertuples(index=False):
        A(f"| {r.gene} | {r.transcript_id_1} | {r.transcript_id_2} | "
          f"{r.n_donors_ge_min_reads} | {r.n_sites} | {r.n_isoforms_informative} | "
          f"{'**yes**' if r.testable else 'no'} |")
    A("")
    (dest / "ASE_SCREEN.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("screen",), default="screen")
    ap.add_argument("--region", choices=REGIONS, default="caudate")
    ap.add_argument("--genes", type=str, default=None,
                    help="comma-separated bare ENSG ids; default = the audited "
                         "colocalization nominations")
    ap.add_argument("--min-donors", type=int, default=MIN_DONORS_INFORMATIVE)
    ap.add_argument("--min-reads", type=int, default=MIN_READS_PER_DONOR)
    ap.add_argument("--reuse", action="store_true",
                    help="reuse the cached assigned sites instead of re-streaming the "
                         "multi-GB allelic-count file (thresholds still re-applied)")
    args = ap.parse_args(argv)

    if args.genes:
        genes = {g.strip() for g in args.genes.split(",") if g.strip()}
    else:
        f = stage_out("anchoring", "locus_event_audit", "abf", "locus_event_audit.parquet")
        if not f.exists():
            raise SystemExit(f"no gene list given and {f} is absent; "
                             f"run locus_event_audit first or pass --genes")
        genes = set(pd.read_parquet(f)["gene"].astype(str))
        print(f"  gene list: {len(genes)} audited colocalization nominations")
    run_screen(args.region, genes=genes, min_donors=args.min_donors,
               min_reads=args.min_reads, reuse=args.reuse)


if __name__ == "__main__":
    main()
