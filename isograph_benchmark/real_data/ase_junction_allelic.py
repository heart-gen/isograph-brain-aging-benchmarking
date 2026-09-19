"""Allele-specific switch test on the junction recount (PI item 10a, step 4).

WHAT THIS IS FOR
----------------
Colocalization says a risk variant and a switch share a causal variant; it does not say the
ALLELE drives the switch. In a donor heterozygous at the switch-QTL lead, the two
haplotypes share a nucleus, a cell-type mixture and an environment, so a cis effect of the
lead on isoform choice shows up as a difference between them. `ase_junction_switch.py`
counted the fragments that say both things at once -- an isoform-specific junction (which
isoform) and a phased heterozygous site (which haplotype) -- and its pre-registered gate
decided whether this step may run. It runs only where that gate passed.

ORIENTATION
-----------
`hap` in the counts is phASER's genome-wide haplotype `PW` (hap 1 = the left allele). `PW`
is anchored to the same population phasing the genotype VCF carries (DLPFC: `PW` equals the
input `GT` at 9,688,924 of 9,688,947 records of a checked sample), so the lead's phased
`GT` is in the same frame: `1|0` puts ALT on hap 1, `0|1` on hap 2. Donors homozygous at
the lead have no ALT haplotype; they keep hap 1 as an arbitrary pseudo-ALT and are the
built-in null.

THE MODEL
---------
Units are donor x haplotype; y = T1 fragments, n = T1 + T2 fragments, x = 1 on ALT.

    logit p = mu + u_donor + beta * x,   u_donor ~ N(0, sigma^2),   y ~ BetaBinomial(n, p, rho)

beta is the within-donor log odds of T1 on the ALT haplotype against the REF haplotype,
tested by likelihood ratio. The donor intercept is RANDOM on purpose: with two units per
donor a fixed intercept is an incidental parameter and the unconditional MLE of beta is
biased away from zero. x is balanced within every donor, so the donor spread (which
carries each donor's cell composition and baseline isoform ratio) is separated from beta by
design. rho absorbs overdispersion within a haplotype.

Alongside it, never instead of it:
  * a model-free stratified score test (per donor U = T1_alt - n_alt * T1 / n, Z = sum U /
    sqrt(sum U^2)), valid under any within-donor overdispersion;
  * the same GLMM on donors homozygous at the lead, which should be null;
  * the between-donor Spearman of lead dosage vs each donor's T1 fraction, from the same
    reads (hap-0 fragments included), the population-level counterpart of beta.

ISOGRAPH MODULES
----------------
Every switch pair comes from a gene of an age-selected IsoGraph co-switching module. Each
pair carries its module, the gene's module role, the module's age association and its
module polarity r(T1) - r(T2) (the transcripts' correlations with the module score). The
module-direction check asks, orientation-free, whether beta tracks that polarity across a
gene's pairs -- whether the lead moves the isoforms along the module's own switch axis --
against a within-gene permutation null. The test itself says a module gene's switch is
cis-regulated; it does not explain why genes co-switch.

  --stage test   Orient, fit, write `allelic_test.parquet`, `allelic_donor_counts.parquet`,
                 `allelic_summary.json` and `ASE_JUNCTION_ALLELIC.md` beside the screen.

Pairs tested: those passing the screen's gate whose lead resolved in the phased VCF, in the
gate family or coloc-nominated. BH is computed over the fitted gate-family pairs only.
beta is oriented to the QTL ALT allele; `--risk-alleles` (rsid, risk_allele; built on
Bridges-2 with `coloc_direction._gwas_risk`) re-orients it to the GWAS risk allele.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.real_data.ase_junction_switch import _store, out_dir
from isograph_benchmark.real_data.ase_switch_direction import MIN_READS_PER_DONOR, REGIONS

# Pre-specified before the test was run: a pair is fitted only with at least this many
# lead-heterozygous donors informative on BOTH haplotypes; below it there are too few
# within-donor contrasts to separate beta from the donor spread.
MIN_PAIRED_HET_DONORS = 10
QVAL_ALLELIC = 0.05

_HET = {"0|1": 2, "1|0": 1}   # phased lead GT -> the haplotype carrying ALT (hap 1 = left)
_GENOTYPES = ("0|0", "0|1", "1|0", "1|1")

# Gauss-Hermite nodes for the donor random intercept.
_GH_T, _GH_W = np.polynomial.hermite.hermgauss(24)
_GH_LOGW = np.log(_GH_W / np.sqrt(np.pi))


# --------------------------------------------------------------------------- #
# Orientation
# --------------------------------------------------------------------------- #
def orient(counts: pd.DataFrame, lead_gt: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    """Per pair x donor: T1/T2 fragments on the ALT and REF haplotype at the pair's lead.

    Homozygous donors get `alt_hap` = 0 and hap 1 as the pseudo-ALT.
    """
    inf = counts[counts["hap"].isin([1, 2])]
    wide = inf.pivot_table(index=["pair_id", "donor_id"], columns=["isoform", "hap"],
                           values="n_frag", aggfunc="sum", fill_value=0)
    wide = wide.reindex(columns=pd.MultiIndex.from_product([[1, 2], [1, 2]]), fill_value=0)
    wide.columns = [f"t{i}_h{h}" for i, h in wide.columns]
    wide = wide.reset_index().merge(pairs[["pair_id", "gene_id", "variant_id_all"]],
                                    on="pair_id")
    g = lead_gt[["variant_id", "donor_id", "gt"]].drop_duplicates(["variant_id", "donor_id"])
    wide = wide.merge(g, left_on=["variant_id_all", "donor_id"],
                      right_on=["variant_id", "donor_id"]).drop(columns="variant_id")
    wide = wide[wide["gt"].isin(_GENOTYPES)].copy()
    wide["alt_hap"] = wide["gt"].map(_HET).fillna(0).astype(int)
    alt1 = (wide["alt_hap"] != 2).to_numpy()
    for iso in (1, 2):
        a, b = wide[f"t{iso}_h1"].to_numpy(), wide[f"t{iso}_h2"].to_numpy()
        wide[f"t{iso}_alt"] = np.where(alt1, a, b)
        wide[f"t{iso}_ref"] = np.where(alt1, b, a)
    return wide.drop(columns=[f"t{i}_h{h}" for i in (1, 2) for h in (1, 2)])


# --------------------------------------------------------------------------- #
# The model
# --------------------------------------------------------------------------- #
def _bb_ll(y, n, p, rho):
    """Beta-binomial log-likelihood without the binomial coefficient (constant in theta)."""
    from scipy.special import betaln
    s = (1.0 - rho) / rho
    a, b = p * s, (1.0 - p) * s
    return betaln(y + a, n - y + b) - betaln(a, b)


def _glmm_nll(theta, y, n, x, starts):
    from scipy.special import expit, logsumexp
    mu, log_sigma, beta, lrho = theta
    u = np.sqrt(2.0) * np.exp(log_sigma) * _GH_T
    p = np.clip(expit(mu + beta * x[:, None] + u[None, :]), 1e-10, 1 - 1e-10)
    ll = _bb_ll(y[:, None], n[:, None], p, expit(lrho))
    per_donor = np.add.reduceat(ll, starts, axis=0)
    return -float(logsumexp(per_donor + _GH_LOGW[None, :], axis=1).sum())


def fit_bb_glmm(y, n, x, donor) -> dict:
    """Beta-binomial GLMM with a donor random intercept; LRT for beta = 0."""
    from scipy.optimize import minimize
    from scipy.stats import chi2
    order = np.argsort(np.asarray(donor), kind="stable")
    y, n, x, donor = (np.asarray(v)[order] for v in (y, n, x, donor))
    y, n, x = y.astype(float), n.astype(float), x.astype(float)
    starts = np.flatnonzero(np.r_[True, donor[1:] != donor[:-1]])
    pbar = np.clip(y.sum() / n.sum(), 0.02, 0.98)
    x0 = np.array([np.log(pbar / (1 - pbar)), np.log(0.5), 0.0, np.log(0.05 / 0.95)])
    bounds = [(-15, 15), (-8, 3), (-10, 10), (-12, 4)]
    opts = {"maxiter": 500}

    def nll0(t):
        return _glmm_nll(np.array([t[0], t[1], 0.0, t[2]]), y, n, x, starts)

    null = minimize(nll0, x0[[0, 1, 3]], method="L-BFGS-B",
                    bounds=[bounds[i] for i in (0, 1, 3)], options=opts)
    full = minimize(lambda t: _glmm_nll(t, y, n, x, starts),
                    np.array([null.x[0], null.x[1], 0.0, null.x[2]]),
                    method="L-BFGS-B", bounds=bounds, options=opts)
    beta = float(full.x[2])
    lrt = max(0.0, 2.0 * (null.fun - full.fun))
    return {"beta": beta,
            # signed-root SE: the Wald SE that reproduces the LRT
            "se": abs(beta) / np.sqrt(lrt) if lrt > 0 else float("nan"),
            "lrt": lrt, "pval": float(chi2.sf(lrt, 1)),
            "sigma_donor": float(np.exp(full.x[1])),
            "rho": float(1.0 / (1.0 + np.exp(-full.x[3]))),
            "converged": bool(full.success and null.success),
            "at_bound": bool(abs(beta) >= bounds[2][1] - 1e-6)}


def robust_score(tab: pd.DataFrame) -> dict:
    """Model-free stratified score test; mean-zero under the null at any overdispersion."""
    from scipy.stats import norm
    n_alt = tab["t1_alt"] + tab["t2_alt"]
    n_tot = n_alt + tab["t1_ref"] + tab["t2_ref"]
    u = tab["t1_alt"] - n_alt * (tab["t1_alt"] + tab["t1_ref"]) / n_tot
    ss = float((u ** 2).sum())
    z = float(u.sum() / np.sqrt(ss)) if ss > 0 else float("nan")
    return {"z": z, "pval": float(2 * norm.sf(abs(z))) if np.isfinite(z) else float("nan")}


def _units(tab: pd.DataFrame) -> tuple[np.ndarray, ...]:
    """Donor x haplotype units with at least one fragment."""
    y = np.r_[tab["t1_alt"].to_numpy(), tab["t1_ref"].to_numpy()]
    n = y + np.r_[tab["t2_alt"].to_numpy(), tab["t2_ref"].to_numpy()]
    x = np.r_[np.ones(len(tab)), np.zeros(len(tab))]
    d = np.r_[np.arange(len(tab)), np.arange(len(tab))]
    keep = n > 0
    return y[keep], n[keep], x[keep], d[keep]


def _is_paired(tab: pd.DataFrame) -> pd.Series:
    return ((tab["t1_alt"] + tab["t2_alt"]) > 0) & ((tab["t1_ref"] + tab["t2_ref"]) > 0)


def allelic_contrast(tab: pd.DataFrame, min_paired: int = MIN_PAIRED_HET_DONORS) -> dict:
    """One pair: the lead-heterozygous test, then the homozygous-at-lead null."""
    out: dict = {}
    for arm, sub in (("het", tab[tab["alt_hap"] > 0]), ("hom", tab[tab["alt_hap"] == 0])):
        paired = _is_paired(sub)
        out[f"n_donors_{arm}"] = int(len(sub))
        out[f"n_paired_{arm}"] = int(paired.sum())
        y, n, x, d = _units(sub)
        ok = bool(paired.sum() >= min_paired and 0 < y.sum() < n.sum())
        if arm == "het":
            tot = sub[["t1_alt", "t2_alt", "t1_ref", "t2_ref"]].sum()
            out.update({f"frag_{k}": int(v) for k, v in tot.items()})
            out["status"] = "fitted" if ok else (
                "too_few_paired_het_donors" if paired.sum() < min_paired
                else "one_isoform_in_het_donors")
            rs = (robust_score(sub[paired]) if paired.any()
                  else {"z": float("nan"), "pval": float("nan")})
            out["score_z"], out["score_pval"] = rs["z"], rs["pval"]
        if not ok:
            continue
        f = fit_bb_glmm(y, n, x, d)
        if arm == "het":
            out.update({k: f[k] for k in ("beta", "se", "lrt", "pval", "sigma_donor", "rho",
                                          "converged", "at_bound")})
        else:
            out["beta_hom_null"], out["pval_hom_null"] = f["beta"], f["pval"]
    return out


def between_donor(iso: pd.DataFrame, dosage: pd.Series, min_reads: int) -> dict:
    """Spearman of lead ALT dosage vs the donor's T1 fraction over all its junction reads."""
    from scipy.stats import spearmanr
    t = iso.set_index("donor_id")
    n = t["t1_all"] + t["t2_all"]
    ds = dosage.reindex(t.index)
    keep = (n >= min_reads) & ds.notna()
    if keep.sum() < 10 or ds[keep].nunique() < 2:
        return {"rho_between": float("nan"), "pval_between": float("nan")}
    r = spearmanr(ds[keep], (t["t1_all"] / n)[keep])
    return {"rho_between": float(r.statistic), "pval_between": float(r.pvalue)}


def _work(job):
    pair_id, tab, iso, dosage, min_paired, min_reads = job
    out = {"pair_id": pair_id, **allelic_contrast(tab, min_paired)}
    out.update(between_donor(iso, dosage, min_reads) if iso is not None
               else {"rho_between": float("nan"), "pval_between": float("nan")})
    return out


def _bh(p: pd.Series) -> pd.Series:
    from scipy.stats import false_discovery_control
    q = pd.Series(np.nan, index=p.index)
    v = p.dropna()
    if len(v):
        q.loc[v.index] = false_discovery_control(v.to_numpy())
    return q


def _lambda_gc(p: pd.Series) -> float | None:
    from scipy.stats import chi2
    p = p.dropna()
    return float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)) if len(p) else None


def lead_site_distance(counts_dir: Path, pairs: pd.DataFrame, lead_gt: pd.DataFrame,
                       tested: set[str]) -> pd.DataFrame:
    """Fragment-weighted median |het site - lead| in lead-het donors (phase-switch exposure)."""
    lead_pos = lead_gt.drop_duplicates("variant_id").set_index("variant_id")["pos"]
    het = lead_gt[lead_gt["gt"].isin(list(_HET))][["variant_id", "donor_id"]]
    rows = []
    for f in sorted(counts_dir.glob("*.pair_sites.tsv.gz")):
        s = pd.read_csv(f, sep="\t")
        s = s[s["pair_id"].isin(tested) & s["hap"].isin([1, 2])]
        if len(s):
            rows.append(s)
    if not rows:
        return pd.DataFrame(columns=["pair_id", "median_lead_site_kb"])
    s = pd.concat(rows, ignore_index=True).merge(
        pairs[["pair_id", "variant_id_all"]], on="pair_id")
    s = s.merge(het, left_on=["variant_id_all", "donor_id"],
                right_on=["variant_id", "donor_id"])
    s["dist_kb"] = (s["het_pos"] - s["variant_id_all"].map(lead_pos)).abs() / 1e3

    def wmed(g):
        g = g.sort_values("dist_kb")
        c = g["n_frag"].cumsum().to_numpy()
        return float(g["dist_kb"].to_numpy()[np.searchsorted(c, c[-1] / 2)])

    return (s.groupby("pair_id")[["dist_kb", "n_frag"]].apply(wmed)
            .rename("median_lead_site_kb").reset_index())


# --------------------------------------------------------------------------- #
# IsoGraph modules: which co-switching program a pair belongs to, and which way it moves
# --------------------------------------------------------------------------- #
def module_annotation(region: str) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Gene -> module/role, (module, transcript) -> polarity r, module -> age association.

    Polarity is `transcript_polarity_table.r`, a transcript's correlation with its module
    score, written by module interpretation for the age-selected modules only (the modules
    the switch pairs were drawn from).
    """
    root = _store(region) / "isograph_vae"
    roles = pd.read_parquet(root / "module_gene_roles.parquet",
                            columns=["gene_id", "module_id", "module_role"])
    pol = []
    for f in sorted((root / "module_interpret").glob("M*/transcript_polarity_table.parquet")):
        p = pd.read_parquet(f, columns=["transcript_id", "r"])
        pol.append(p.assign(module_id=f.parent.name))
    polarity = (pd.concat(pol).set_index(["module_id", "transcript_id"])["r"] if pol
                else pd.Series(dtype=float))
    age = pd.read_parquet(root / "module_interpret" / "module_selection.parquet",
                          columns=["module_id", "trait", "effect"]).rename(
        columns={"trait": "module_age_trait", "effect": "module_age_effect"})
    return roles, polarity, age


def add_module_columns(out: pd.DataFrame, roles: pd.DataFrame, polarity: pd.Series,
                       age: pd.DataFrame) -> pd.DataFrame:
    """module_polarity = r(T1) - r(T2): > 0 when T1 rises over T2 as the module score rises.

    beta_along_module = beta * sign(module_polarity) is the ALT allele's effect along the
    module's switch direction. ALT is arbitrary ACROSS genes, so its sign only means
    something within a gene (see `module_direction_check`) or once oriented to a risk
    allele (`risk_along_module`).
    """
    out = out.merge(roles, on="gene_id", how="left").merge(age, on="module_id", how="left")

    def r(mod, tx):
        return polarity.get((mod, tx), np.nan) if isinstance(mod, str) else np.nan

    out["module_polarity"] = [r(m, a) - r(m, b) for m, a, b in
                              zip(out["module_id"], out["transcript_id_1"],
                                  out["transcript_id_2"])]
    s = np.sign(out["module_polarity"]).replace(0, np.nan)
    out["beta_along_module"] = out["beta"] * s
    out["risk_along_module"] = out["beta_risk"] * s
    return out


def module_direction_check(fitted: pd.DataFrame, qval: float = QVAL_ALLELIC,
                           min_pairs: int = 3, n_perm: int = 1000, seed: int = 0) -> dict:
    """Does the lead move a gene's isoforms along its module's switch axis?

    Orientation-free: within a gene, ALT is one fixed allele, so the question is whether
    beta tracks module polarity across that gene's pairs. Statistic: the median over genes
    (>= `min_pairs` fitted pairs, polarity not constant) of |Spearman(beta, polarity)|;
    null: polarity permuted within each gene. Pairs of a gene share transcripts, so the
    permutation (not the pair count) is the reference. Also: among genes with >= 2
    significant pairs, how many have sign(beta * polarity) the same for all of them.
    """
    from scipy.stats import rankdata
    f = fitted.dropna(subset=["beta", "module_polarity"])
    genes = []
    for _, g in f.groupby("gene_id"):
        if len(g) >= min_pairs and g["module_polarity"].nunique() > 1 and g["beta"].nunique() > 1:
            genes.append((rankdata(g["beta"]), rankdata(g["module_polarity"])))
    if not genes:
        return {"n_genes": 0}

    def abs_rho(a, b):   # rows of b are permutations; Pearson on ranks = Spearman
        a = (a - a.mean()) / a.std()
        b = (b - b.mean(axis=-1, keepdims=True)) / b.std(axis=-1, keepdims=True)
        return np.abs((a * b).mean(axis=-1))

    obs = float(np.median([abs_rho(a, b[None, :])[0] for a, b in genes]))
    rng = np.random.default_rng(seed)
    null = np.median(np.stack([abs_rho(a, rng.permuted(np.tile(b, (n_perm, 1)), axis=1))
                               for a, b in genes]), axis=0)
    sig = f[f["qval"] < qval]
    consistent, multi = 0, 0
    for _, g in sig.groupby("gene_id"):
        prod = np.sign(g["beta"] * g["module_polarity"])
        prod = prod[prod != 0]
        if len(prod) >= 2:
            multi += 1
            consistent += int(prod.nunique() == 1)
    return {"n_genes": len(genes), "min_pairs": min_pairs,
            "median_abs_rho": obs, "null_median_abs_rho": float(np.mean(null)),
            "null_q95": float(np.quantile(null, 0.95)),
            "perm_p": float((1 + (null >= obs).sum()) / (1 + n_perm)),
            "n_genes_ge2_sig_pairs": multi, "n_genes_all_sig_pairs_aligned": consistent}


def module_table(out: pd.DataFrame) -> list[dict]:
    """Per module: gate-family pairs/genes fitted and significant."""
    f = out[(out["status"] == "fitted") & out["gate_family"]]
    rows = []
    for mod, g in f.groupby("module_id"):
        sig = g[g["qval"] < QVAL_ALLELIC]
        rows.append({"module_id": mod, "module_age_trait": g["module_age_trait"].iloc[0],
                     "module_age_effect": float(g["module_age_effect"].iloc[0]),
                     "pairs_fitted": int(len(g)), "genes_fitted": int(g["gene_id"].nunique()),
                     "pairs_q05": int(len(sig)), "genes_q05": int(sig["gene_id"].nunique())})
    return rows


def _coloc_traits() -> dict[str, str]:
    """Bare gene id -> comma-joined GWAS traits it colocalizes with (BrainSEQ S_g coloc
    nominations and the GTEx locus-event audit, the same sources as `coloc_nominated`)."""
    from isograph_benchmark.paths import stage_out
    frames = []
    for f in (stage_out("anchoring", "coloc_brainseq", "ea_only", "nominations.parquet"),
              stage_out("anchoring", "locus_event_audit", "abf", "locus_event_audit.parquet"),
              stage_out("anchoring", "locus_event_audit", "susie", "locus_event_audit.parquet")):
        if f.exists():
            frames.append(pd.read_parquet(f, columns=["gene", "trait"]))
    if not frames:
        return {}
    t = pd.concat(frames).dropna()
    t["gene"] = t["gene"].astype(str).str.split(".").str[0]
    return t.groupby("gene")["trait"].agg(lambda s: ",".join(sorted(set(s)))).to_dict()


# --------------------------------------------------------------------------- #
# Stage: test
# --------------------------------------------------------------------------- #
def run_test(region: str, min_paired: int, min_reads: int, risk_alleles: Path | None,
             threads: int) -> Path:
    from concurrent.futures import ProcessPoolExecutor
    dest = out_dir(region)
    summ = json.loads((dest / "screen_summary.json").read_text())
    if not summ["gate_passed"]:
        # the pre-specified negative is a result, not an error: end the SLURM chain cleanly
        print(f"  {region}: the screen gate failed; the allelic test is not run (pre-specified)")
        return dest
    counts = pd.read_parquet(dest / "junction_allelic_counts.parquet")
    feas = pd.read_parquet(dest / "pair_feasibility.parquet")
    lead_gt = pd.read_parquet(dest / "lead_genotypes.parquet")
    feas["gate_family"] = feas["gate_family"].fillna(False).astype(bool)
    feas["coloc_nominated"] = feas["coloc_brainseq"] | feas["coloc_audit"]
    sel = feas[feas["testable"] & feas["lead_resolved"]
               & (feas["gate_family"] | feas["coloc_nominated"])].copy()
    counts = counts[counts["pair_id"].isin(sel["pair_id"])]

    oriented = orient(counts, lead_gt, sel)
    oriented["dosage"] = oriented["gt"].str.count("1")
    oriented.to_parquet(dest / "allelic_donor_counts.parquet", index=False)

    # between-donor check: every junction-informative fragment (hap 0 too), every donor
    # genotyped at the lead
    iso = (counts.groupby(["pair_id", "donor_id", "isoform"])["n_frag"].sum()
           .unstack("isoform").reindex(columns=[1, 2]).fillna(0))
    iso.columns = ["t1_all", "t2_all"]
    iso = dict(tuple(iso.reset_index().groupby("pair_id")))
    g = lead_gt[lead_gt["gt"].isin(_GENOTYPES)].drop_duplicates(["variant_id", "donor_id"])
    dos = {v: d.set_index("donor_id")["gt"].str.count("1")
           for v, d in g.groupby("variant_id")}
    lead_of = sel.set_index("pair_id")["variant_id_all"]

    jobs = [(pid, tab, iso.get(pid), dos.get(lead_of[pid], pd.Series(dtype=float)),
             min_paired, min_reads) for pid, tab in oriented.groupby("pair_id")]
    print(f"  {len(jobs):,} pairs with oriented donors ({region}); {threads} worker(s)")
    t0 = time.time()
    if threads > 1:
        with ProcessPoolExecutor(threads) as ex:
            res = list(ex.map(_work, jobs, chunksize=4))
    else:
        res = [_work(j) for j in jobs]
    print(f"  fitted in {time.time() - t0:.0f}s")

    out = sel[["pair_id", "gene_id", "transcript_id_1", "transcript_id_2", "variant_id_all",
               "slope_all", "qval_all", "gate_family", "coloc_nominated", "coloc_brainseq",
               "coloc_audit", "two_sided", "n_donors_ge_min_reads"]].merge(
        pd.DataFrame(res), on="pair_id", how="left")
    out["status"] = out["status"].fillna("no_oriented_donors")
    # everything the Bridges-2 side needs to orient beta to a GWAS risk allele from this
    # table alone: the lead's alleles, whether they are strand-ambiguous, and the traits
    alleles = lead_gt.drop_duplicates("variant_id").set_index("variant_id")
    out["lead_ref"] = out["variant_id_all"].map(alleles["ref"])
    out["lead_alt"] = out["variant_id_all"].map(alleles["alt"])
    out["palindromic"] = [{a, b} in ({"A", "T"}, {"C", "G"})
                          for a, b in zip(out["lead_ref"], out["lead_alt"])]
    out["coloc_traits"] = out["gene_id"].astype(str).str.split(".").str[0].map(_coloc_traits())
    fit = out["status"] == "fitted"
    out["qval"] = np.nan
    gf = fit & out["gate_family"]
    out.loc[gf, "qval"] = _bh(out.loc[gf, "pval"])
    both = fit & out["rho_between"].notna()
    out["sign_agrees_between"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    out.loc[both, "sign_agrees_between"] = (np.sign(out.loc[both, "beta"])
                                            == np.sign(out.loc[both, "rho_between"]))
    out = out.merge(lead_site_distance(dest / "counts", sel, lead_gt,
                                       set(out.loc[fit, "pair_id"])), on="pair_id", how="left")

    out["risk_allele"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out["beta_risk"] = np.nan
    if risk_alleles is not None:
        ra = pd.read_csv(risk_alleles, sep="\t", usecols=["rsid", "risk_allele"])
        ra = ra.drop_duplicates("rsid").set_index("rsid")["risk_allele"].str.upper()
        out["risk_allele"] = out["variant_id_all"].map(ra).astype("string")
        on_alt = out["risk_allele"] == out["lead_alt"]
        on_ref = out["risk_allele"] == out["lead_ref"]
        out["beta_risk"] = out["beta"] * np.where(on_alt.fillna(False), 1.0,
                                                  np.where(on_ref.fillna(False), -1.0, np.nan))

    out = add_module_columns(out, *module_annotation(region))
    out = out.sort_values(["gate_family", "pval"], ascending=[False, True])
    out.to_parquet(dest / "allelic_test.parquet", index=False)

    fitted = out[out["status"] == "fitted"]
    summary = {
        "region": region, "min_paired_het_donors": min_paired,
        "min_reads_between": min_reads, "qval_threshold": QVAL_ALLELIC,
        "families": {"gate_family": _family(out[out["gate_family"]]),
                     "coloc_nominated": _family(out[out["coloc_nominated"]])},
        "status_counts": {k: int(v) for k, v in out["status"].value_counts().items()},
        "hom_null": {"n_fitted": int(fitted["pval_hom_null"].notna().sum()),
                     "frac_p05": _frac(fitted["pval_hom_null"] < 0.05,
                                       fitted["pval_hom_null"].notna()),
                     "lambda_gc": _lambda_gc(fitted["pval_hom_null"])},
        "lambda_gc_het": _lambda_gc(fitted["pval"]),
        "score_vs_glmm_spearman": (float(fitted[["score_z", "beta"]].corr("spearman").iloc[0, 1])
                                   if len(fitted) > 2 else None),
        "median_lead_site_kb": (float(fitted["median_lead_site_kb"].median())
                                if fitted["median_lead_site_kb"].notna().any() else None),
        "risk_alleles": str(risk_alleles) if risk_alleles else None,
        "modules": module_table(out),
        "module_direction": module_direction_check(fitted[fitted["gate_family"]]),
    }
    (dest / "allelic_summary.json").write_text(json.dumps(summary, indent=2))
    _write_report(dest, region, out, summary)
    print(json.dumps(summary, indent=2))
    return dest


def _frac(hit: pd.Series, among: pd.Series) -> float | None:
    return float(hit[among].mean()) if among.any() else None


def _family(s: pd.DataFrame) -> dict:
    f = s[s["status"] == "fitted"]
    sig = f[f["qval"] < QVAL_ALLELIC]
    agree = f["sign_agrees_between"].dropna().astype(bool)
    sig_agree = sig["sign_agrees_between"].dropna().astype(bool)
    return {"n_pairs": int(len(s)), "n_fitted": int(len(f)),
            "n_genes_fitted": int(f["gene_id"].nunique()),
            "n_nominal_p05": int((f["pval"] < 0.05).sum()),
            "n_q05": int(len(sig)), "n_genes_q05": int(sig["gene_id"].nunique()),
            "sign_agreement_between": float(agree.mean()) if len(agree) else None,
            "sign_agreement_between_q05": float(sig_agree.mean()) if len(sig_agree) else None,
            "n_not_converged": int((~f["converged"].astype(bool)).sum()),
            "n_at_bound": int(f["at_bound"].astype(bool).sum()),
            "n_at_bound_q05": int(sig["at_bound"].astype(bool).sum())}


def _fmt(v, spec: str) -> str:
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, spec)


def _write_report(dest: Path, region: str, out: pd.DataFrame, s: dict) -> None:
    L: list[str] = []
    A = L.append
    A(f"# Allele-specific switch test — {region}")
    A("")
    A("Within a donor heterozygous at the switch-QTL lead, the two haplotypes share a nucleus, "
      "a cell-type mixture and an environment, so a cis effect of the lead on isoform choice "
      "shows up as a difference between them. Each fragment crosses an isoform-specific "
      "junction (which isoform) and carries a phased heterozygous site (which haplotype); "
      "the haplotype carrying the lead's ALT allele is read off the donor's phased genotype "
      "at the lead, in the same phase frame as phASER's `PW`.")
    A("")
    A("## Model")
    A("")
    A("- units: donor x haplotype; `y` = T1 fragments, `n` = T1 + T2 fragments")
    A("- `logit p = mu + u_donor + beta * ALT`, `u_donor ~ N(0, sigma^2)`, "
      "`y ~ BetaBinomial(n, p, rho)`; likelihood-ratio test of `beta = 0`")
    A("- `beta` = within-donor log odds of T1 on the ALT haplotype vs the REF haplotype; the "
      "random donor intercept absorbs each donor's baseline isoform ratio, and with it the "
      "donor's cell composition")
    A(f"- fitted when >= {s['min_paired_het_donors']} lead-heterozygous donors have fragments "
      "on both haplotypes; BH over the fitted gate-family pairs")
    A("- built-in null: the same model on donors HOMOZYGOUS at the lead (hap 1 as a "
      "pseudo-ALT)")
    A("- sensitivity: a model-free stratified score test with empirical variance, and the "
      "between-donor Spearman of lead dosage vs the donor's T1 fraction from the same reads")
    A("")
    A("## Result")
    A("")
    A("| family | pairs | fitted | genes fitted | p < 0.05 | **q < 0.05** | genes q < 0.05 | "
      "sign = between-donor (all / q < 0.05) |")
    A("|---|---|---|---|---|---|---|---|")
    for name, f in s["families"].items():
        A(f"| {name} | {f['n_pairs']:,} | {f['n_fitted']:,} | {f['n_genes_fitted']:,} | "
          f"{f['n_nominal_p05']:,} | **{f['n_q05']:,}** | {f['n_genes_q05']:,} | "
          f"{_fmt(f['sign_agreement_between'], '.2f')} / "
          f"{_fmt(f['sign_agreement_between_q05'], '.2f')} |")
    A("")
    A("q-values are computed within the gate family only; a coloc-nominated pair outside the "
      "gate family is fitted but carries no q. Pairs of one gene share donors and fragments, "
      "so pair-level counts are not independent; the gene counts are the conservative read.")
    A("")
    A("## Calibration")
    A("")
    h = s["hom_null"]
    A(f"- homozygous-at-lead null: {h['n_fitted']:,} fits, fraction p < 0.05 = "
      f"{_fmt(h['frac_p05'], '.3f')}, lambda_GC = {_fmt(h['lambda_gc'], '.2f')}")
    A(f"- lead-heterozygous test: lambda_GC = {_fmt(s['lambda_gc_het'], '.2f')} (inflation is "
      "expected here where the leads are real switch-QTLs)")
    A(f"- GLMM beta vs score-test Z, Spearman: {_fmt(s['score_vs_glmm_spearman'], '.2f')}")
    A(f"- median lead-to-het-site distance: {_fmt(s['median_lead_site_kb'], '.0f')} kb")
    gf = s["families"]["gate_family"]
    A(f"- gate family: {gf['n_not_converged']} fits did not converge; {gf['n_at_bound']} hit "
      f"the |beta| = 10 bound ({gf['n_at_bound_q05']} of them at q < 0.05). A bound fit is "
      "quasi-separation (one haplotype carries one isoform only): its sign is informative, "
      "its magnitude is not.")
    A("")
    A("## IsoGraph modules")
    A("")
    A("Every switch pair belongs to a gene of an age-selected IsoGraph co-switching module. "
      "The within-donor test asks whether that gene's switch is under cis-genetic control "
      "that cell composition cannot produce; it does not explain why genes co-switch.")
    A("")
    A("| module | age trait | age effect | pairs fitted | genes fitted | pairs q < 0.05 | "
      "genes q < 0.05 |")
    A("|---|---|---|---|---|---|---|")
    for m in s["modules"]:
        A(f"| {m['module_id']} | {m['module_age_trait']} | {m['module_age_effect']:+.2f} | "
          f"{m['pairs_fitted']:,} | {m['genes_fitted']:,} | {m['pairs_q05']:,} | "
          f"{m['genes_q05']:,} |")
    A("")
    md = s["module_direction"]
    if md.get("n_genes"):
        A(f"**Does the lead move the isoforms along the module's switch axis?** Within a gene, "
          f"ALT is one fixed allele, so beta should track each pair's module polarity "
          f"(r(T1) - r(T2) against the module score) if the variant and the module move the "
          f"same axis. Over {md['n_genes']} gate-family genes with >= {md['min_pairs']} "
          f"fitted pairs, median |Spearman(beta, polarity)| = {md['median_abs_rho']:.2f}, "
          f"against {md['null_median_abs_rho']:.2f} with polarity permuted within gene "
          f"(95th percentile {md['null_q95']:.2f}; permutation p = {md['perm_p']:.3g}). "
          f"Among genes with >= 2 significant pairs, {md['n_genes_all_sig_pairs_aligned']} "
          f"of {md['n_genes_ge2_sig_pairs']} have every significant pair on the same side "
          "of the module axis.")
        A("")
    A("`beta_along_module` in `allelic_test.parquet` is beta signed by the pair's module "
      "polarity; across genes its sign is arbitrary until oriented to a risk allele "
      "(`risk_along_module`: > 0 when the risk allele shifts the isoforms the way the "
      "module score rises).")
    A("")
    A("## Caveats")
    A("")
    A("- Cell composition cancels within a donor only if allelic effects do not differ by "
      "cell type.")
    A("- Some reference mapping bias survives WASP. It shifts ALT vs REF fragment totals; it "
      "moves the T1 fraction within a haplotype only if it differs between the two isoforms' "
      "junction reads.")
    A("- A statistical phase switch between the lead and a fragment's heterozygous site swaps "
      "ALT and REF for that fragment and pulls beta toward zero; the lead-to-site distance is "
      "recorded per pair.")
    A("- The lead is the all_samples switch-QTL lead, not necessarily the causal variant; a "
      "causal variant in imperfect LD with it attenuates beta.")
    if s["risk_alleles"] is None:
        A("- `beta` is oriented to the QTL ALT allele. Orienting it to the GWAS risk allele "
          "needs the rsID -> risk-allele table (`--risk-alleles`, built on Bridges-2 with "
          "`coloc_direction._gwas_risk`).")
    A("")
    A("## Top gate-family pairs")
    A("")
    A("| gene | T1 | T2 | lead | het donors (paired) | beta | p | q | score p | "
      "between rho | hom-null p |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    top = out[(out["status"] == "fitted") & out["gate_family"]].head(25)
    for r in top.itertuples(index=False):
        A(f"| {r.gene_id} | {r.transcript_id_1} | {r.transcript_id_2} | {r.variant_id_all} | "
          f"{r.n_donors_het} ({r.n_paired_het}) | {r.beta:+.2f}{' (bound)' if r.at_bound else ''} | "
          f"{r.pval:.1e} | "
          f"{r.qval:.1e} | {r.score_pval:.1e} | {_fmt(r.rho_between, '+.2f')} | "
          f"{_fmt(r.pval_hom_null, '.2f')} |")
    A("")
    (dest / "ASE_JUNCTION_ALLELIC.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", required=True, choices=("test",))
    ap.add_argument("--region", choices=REGIONS, default="dlpfc")
    ap.add_argument("--min-paired-het-donors", type=int, default=MIN_PAIRED_HET_DONORS)
    ap.add_argument("--min-reads", type=int, default=MIN_READS_PER_DONOR,
                    help="between-donor check: fragments a donor needs to enter it")
    ap.add_argument("--risk-alleles", type=Path, default=None,
                    help="TSV with rsid, risk_allele (the GWAS risk allele per lead rsID)")
    ap.add_argument("--threads", type=int, default=1, help="worker processes")
    args = ap.parse_args(argv)
    run_test(args.region, args.min_paired_het_donors, args.min_reads, args.risk_alleles,
             args.threads)


if __name__ == "__main__":
    main()
