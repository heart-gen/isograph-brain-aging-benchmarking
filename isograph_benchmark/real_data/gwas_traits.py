"""GWAS trait registry for the genetic-anchoring capstone (S-LDSC + coloc).

Single source of truth for the summary statistics used beyond MAGMA: one `TraitSpec`
per trait carrying the file path, column names, effect encoding, sample size, genome
build, and the LD-confounded regions to drop in fine-mapping. Both the S-LDSC munge
step (`05_genetic_anchoring/_h/03d.ldsc_munge_h2.sh`) and the coloc pipeline read this, so a trait is described once.

Two disease "cases" share the machinery:
  * disease   — schizophrenia (PGC3 wave 3, European), anchored in the SCZD switch layer.
  * aging     — neurodegeneration (AD, PD, LBD, ALS), anchored in the pooled
                aging switch layer. These are the aging-case counterpart to SCZ.

Build note: PGC3 SCZ is GRCh37 (matches the g1000_eur LD panel). The neurodegeneration
sumstats are GRCh38 harmonised GWAS-catalog files; S-LDSC is build-agnostic (munge
matches HapMap3 by rsID), and the coloc pipeline re-derives hg19 positions from the LD
panel bim by rsID, so `build` is informational for both.
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import pandas as pd

from isograph_benchmark.paths import ensure_dir

_GWAS = Path("/ocean/projects/bio250020p/shared/resources/gwas")
_GWAS_2600 = Path("/ocean/projects/bio260021p/shared/resources/gwas")

# LD-confounded regions to exclude from per-locus fine-mapping, in hg19 (the coloc
# panel build). MHC is the long-range-LD sink; APOE is the AD/LBD long-range signal.
MHC_HG19 = (6, 25_000_000, 34_000_000)
APOE_HG19 = (19, 44_400_000, 46_500_000)


@dataclass(frozen=True)
class TraitSpec:
    key: str
    label: str
    case: str                 # "disease" | "aging"
    path: Path
    build: str                # "hg19" | "hg38"
    snp_col: str              # rsID column
    chr_col: str
    pos_col: str
    a1_col: str               # effect allele
    a2_col: str               # other allele
    effect_col: str           # signed effect (beta on log-odds scale)
    effect_null: float        # null value for --signed-sumstats (0 for beta)
    se_col: str
    p_col: str
    n_cas_col: str | None = None
    n_con_col: str | None = None
    fixed_n: int | None = None
    comment: str | None = None            # leading meta-line prefix to strip (PGC VCF)
    exclude_hg19: tuple[tuple[int, int, int], ...] = field(default_factory=lambda: (MHC_HG19,))


TRAITS: dict[str, TraitSpec] = {
    # --- disease case ---------------------------------------------------------
    "scz": TraitSpec(
        key="scz", label="SCZ", case="disease",
        path=_GWAS / "PGC/SCZ/PGC3/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
        build="hg19", comment="##",
        snp_col="ID", chr_col="CHROM", pos_col="POS", a1_col="A1", a2_col="A2",
        effect_col="BETA", effect_null=0.0, se_col="SE", p_col="PVAL",
        n_cas_col="NCAS", n_con_col="NCON",
        exclude_hg19=(MHC_HG19,),
    ),
    # --- aging / neurodegeneration case --------------------------------------
    # AD: Bellenguez et al. 2022 (GCST90027158). Per-SNP case/control N; APOE excluded.
    "ad": TraitSpec(
        key="ad", label="AD", case="aging",
        path=_GWAS / "alz/bellenguez2022/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
        build="hg38",
        snp_col="hm_rsid", chr_col="hm_chrom", pos_col="hm_pos",
        a1_col="hm_effect_allele", a2_col="hm_other_allele",
        effect_col="hm_beta", effect_null=0.0, se_col="standard_error", p_col="p_value",
        n_cas_col="n_cas", n_con_col="n_con",
        exclude_hg19=(MHC_HG19, APOE_HG19),
    ),
    # PD: Nalls et al. 2019, no-23andMe public (GCST009325). Per-SNP case/control N.
    "pd": TraitSpec(
        key="pd", label="PD", case="aging",
        path=_GWAS / "PD/data/GCST009325.h.tsv.gz",
        build="hg38",
        snp_col="rsid", chr_col="chromosome", pos_col="base_pair_location",
        a1_col="effect_allele", a2_col="other_allele",
        effect_col="beta", effect_null=0.0, se_col="standard_error", p_col="p_value",
        n_cas_col="N_cases", n_con_col="N_controls",
        exclude_hg19=(MHC_HG19,),
    ),
    # LBD: Chia et al. 2021 (GCST90001390). No per-SNP N -> fixed total (2591 cases +
    # 4027 controls). APOE excluded (genuine but long-range LD signal).
    "lbd": TraitSpec(
        key="lbd", label="LBD", case="aging",
        path=_GWAS / "LBD/data/harmonised/33589841-GCST90001390-EFO_0006792.h.tsv.gz",
        build="hg38",
        snp_col="hm_rsid", chr_col="hm_chrom", pos_col="hm_pos",
        a1_col="hm_effect_allele", a2_col="hm_other_allele",
        effect_col="hm_beta", effect_null=0.0, se_col="standard_error", p_col="p_value",
        fixed_n=6618,
        exclude_hg19=(MHC_HG19, APOE_HG19),
    ),
    # ALS: van Rheenen et al. 2021 EUR (GCST90027164). No per-SNP N -> fixed total.
    "als": TraitSpec(
        key="als", label="ALS", case="aging",
        path=_GWAS_2600 / "als/34873335-GCST90027164-MONDO_0004976.h.tsv.gz",
        build="hg38",
        snp_col="hm_rsid", chr_col="hm_chrom", pos_col="hm_pos",
        a1_col="hm_effect_allele", a2_col="hm_other_allele",
        effect_col="hm_beta", effect_null=0.0, se_col="standard_error", p_col="p_value",
        fixed_n=138086,
        exclude_hg19=(MHC_HG19,),
    ),
    # FTD intentionally omitted: the only available sumstats are IFGC (Ferrari et al.
    # 2014); the intended cohort is a 2024 FTD GWAS not yet provided. harmonize_ftd.py
    # remains for when those sumstats arrive, but no FTD TraitSpec is registered so it
    # never enters the anchoring runs with a 2014/2024 provenance mismatch.
}

AGING_TRAITS = [k for k, s in TRAITS.items() if s.case == "aging"]
DISEASE_TRAITS = [k for k, s in TRAITS.items() if s.case == "disease"]


def get(trait: str) -> TraitSpec:
    if trait not in TRAITS:
        raise SystemExit(f"Unknown trait {trait!r}; known: {sorted(TRAITS)}")
    spec = TRAITS[trait]
    if not spec.path.exists():
        raise SystemExit(f"Sumstats for {trait!r} not found: {spec.path}")
    return spec


def write_munge_input(spec: TraitSpec, clean_dir: Path) -> Path:
    """Stream the raw sumstats to a canonical `SNP A1 A2 BETA SE P N` tsv.gz.

    Harmonised GWAS-catalog files carry BOTH `hm_*` and plain allele columns, which
    makes munge_sumstats' column auto-detection ambiguous ("2 different A2 columns").
    Resolving the columns we want by name and re-emitting a slim, canonically-named
    file removes that ambiguity and makes the munge flags identical for every trait.
    Also strips the PGC VCF `##` meta lines. Cached by output path.
    """
    ensure_dir(clean_dir)
    clean = clean_dir / f"{spec.key}.munge.tsv.gz"
    if clean.exists():
        return clean
    idx = _resolve_indices(spec)
    body = ('{ if (RS ~ /^rs/ && (SE+0)>0 && (P+0)>0) '
            'print RS"\\t"A1"\\t"A2"\\t"B"\\t"SE"\\t"P"\\t"N }')
    prog = _awk_prog(spec, idx, body,
                     begin='BEGIN{print "SNP\\tA1\\tA2\\tBETA\\tSE\\tP\\tN"}')
    tmp = clean.with_suffix(".tmp.gz")
    zcat = subprocess.Popen(["zcat", str(spec.path)], stdout=subprocess.PIPE)
    awk = subprocess.Popen(["awk", "-F", "\t", prog], stdin=zcat.stdout,
                           stdout=subprocess.PIPE)
    zcat.stdout.close()
    with open(tmp, "wb") as fh:
        gz = subprocess.Popen(["gzip", "-c"], stdin=awk.stdout, stdout=fh)
        awk.stdout.close()
        gz.wait()
    awk.wait(); zcat.wait()
    if gz.returncode or awk.returncode:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"munge-input stream failed for {spec.key}")
    tmp.rename(clean)
    return clean


def munge_args(spec: TraitSpec, clean_dir: Path) -> list[str]:
    """munge_sumstats.py flags for this trait (uniform canonical columns)."""
    sumstats = write_munge_input(spec, clean_dir)
    return [
        "--sumstats", str(sumstats),
        "--snp", "SNP", "--a1", "A1", "--a2", "A2",
        "--signed-sumstats", "BETA,0", "--p", "P", "--N-col", "N",
    ]


# --------------------------------------------------------------------------------
# Build-agnostic streaming for coloc: pull normalised (rsid, a1, a2, beta, se, p, n)
# rows straight from the raw sumstats by resolving column *names* to awk indices, so
# hg38 harmonised files and the hg19 PGC VCF are read the same way. The coloc pipeline
# re-derives hg19 positions from the LD-panel bim by rsID, so no GWAS position is used
# here and no liftover is needed.
# --------------------------------------------------------------------------------
_STREAM_COLS = ["rsid", "a1", "a2", "beta", "se", "p", "n"]


def _resolve_indices(spec: TraitSpec) -> dict[str, int]:
    """1-based column indices for the fields we stream (header resolved by name)."""
    with gzip.open(spec.path, "rt") as fh:
        header = None
        for line in fh:
            if spec.comment and line.startswith(spec.comment):
                continue
            header = line.rstrip("\n").split("\t")
            break
    if header is None:
        raise SystemExit(f"{spec.key}: no header line found")
    want = [spec.snp_col, spec.a1_col, spec.a2_col, spec.effect_col, spec.se_col, spec.p_col]
    if spec.n_cas_col and spec.n_con_col:
        want += [spec.n_cas_col, spec.n_con_col]
    idx = {}
    for c in want:
        if c not in header:
            raise SystemExit(f"{spec.key}: column {c!r} not in header {header[:8]}...")
        idx[c] = header.index(c) + 1
    return idx


def _n_expr(spec: TraitSpec, idx: dict[str, int]) -> str:
    if spec.n_cas_col and spec.n_con_col:
        return f"(${idx[spec.n_cas_col]}+${idx[spec.n_con_col]})"
    return str(spec.fixed_n)


def _awk_prog(spec: TraitSpec, idx: dict[str, int], body: str,
              two_file: bool = False, begin: str = "") -> str:
    """Assemble an awk program that skips meta/header lines then runs `body`.

    `body` may reference the resolved fields via the placeholders RS,A1,A2,B,SE,P,N
    (expanded to `$<index>`), and W for an rsID lookup set (built from a first file when
    `two_file`). Meta/header skipping is gated to the *second* input (NR>FNR) so a
    leading rsID-set file passes through untouched.
    """
    guard = "NR>FNR" if two_file else "1"
    comment_skip = ""
    if spec.comment:
        comment_skip = (f'{guard} && substr($0,1,{len(spec.comment)})=="{spec.comment}"'
                        f'{{next}}\n')
    # header: first non-comment line of the (second) input
    hdr = f'{guard} && hdr==0 {{hdr=1; next}}\n'
    repl = {"RS": f"${idx[spec.snp_col]}", "A1": f"${idx[spec.a1_col]}",
            "A2": f"${idx[spec.a2_col]}", "B": f"${idx[spec.effect_col]}",
            "SE": f"${idx[spec.se_col]}", "P": f"${idx[spec.p_col]}",
            "N": _n_expr(spec, idx)}
    for k, v in repl.items():
        body = body.replace(k, v)
    setload = 'NR==FNR{W[$1]=1; next}\n' if two_file else ""
    begin_block = (begin + "\n") if begin else ""
    return begin_block + setload + comment_skip + hdr + guard + " " + body


def _run_awk(spec: TraitSpec, prog: str, set_file: Path | None) -> pd.DataFrame:
    zcat = subprocess.Popen(["zcat", str(spec.path)], stdout=subprocess.PIPE)
    cmd = ["awk", "-F", "\t", prog]
    cmd += ([str(set_file), "-"] if set_file else ["-"])
    res = subprocess.run(cmd, stdin=zcat.stdout, capture_output=True, text=True)
    zcat.stdout.close(); zcat.wait()
    if res.returncode != 0:
        raise RuntimeError(f"awk failed for {spec.key}: {res.stderr[:400]}")
    if not res.stdout.strip():
        return pd.DataFrame(columns=_STREAM_COLS)
    d = pd.read_csv(StringIO(res.stdout), sep="\t", header=None, names=_STREAM_COLS,
                    dtype={"rsid": str, "a1": str, "a2": str})
    for c in ["beta", "se", "p", "n"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["beta", "se", "p"])


def stream_sig_snps(spec: TraitSpec, p_thresh: float) -> pd.DataFrame:
    """rsID-keyed genome-wide SNPs with P < p_thresh (rsid, a1, a2, beta, se, p, n)."""
    idx = _resolve_indices(spec)
    body = (f'{{ pv=P+0; if (pv>0 && pv<{p_thresh:g} && RS ~ /^rs/) '
            f'print RS"\\t"A1"\\t"A2"\\t"B"\\t"SE"\\t"pv"\\t"N }}')
    return _run_awk(spec, _awk_prog(spec, idx, body), None)


def stream_snps_in_set(spec: TraitSpec, rsids: set[str], tmp_dir: Path) -> pd.DataFrame:
    """rsID-keyed SNP rows for rsid in `rsids` (for per-locus GWAS fine-mapping)."""
    idx = _resolve_indices(spec)
    ensure_dir(tmp_dir)
    set_file = tmp_dir / f"_{spec.key}_rsids.txt"
    set_file.write_text("\n".join(sorted(rsids)) + "\n")
    body = ('{ if ((RS in W) && (SE+0)>0) '
            'print RS"\\t"A1"\\t"A2"\\t"B"\\t"SE"\\t"P"\\t"N }')
    out = _run_awk(spec, _awk_prog(spec, idx, body, two_file=True), set_file)
    set_file.unlink(missing_ok=True)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="GWAS trait registry helpers.")
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("munge-args", help="print munge_sumstats.py flags, one per line")
    m.add_argument("--trait", required=True)
    m.add_argument("--clean-dir", required=True, help="dir for a cleaned copy if needed")

    sub.add_parser("list", help="list registered traits")

    args = p.parse_args()
    if args.cmd == "list":
        for k, s in TRAITS.items():
            print(f"{k}\t{s.label}\t{s.case}\t{s.build}\t{'OK' if s.path.exists() else 'MISSING'}")
    elif args.cmd == "munge-args":
        for a in munge_args(get(args.trait), Path(args.clean_dir)):
            print(a)


if __name__ == "__main__":
    main()
