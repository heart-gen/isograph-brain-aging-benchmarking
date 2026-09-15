"""Do SCZ-risk loci converge on age-sensitive isoform-switch programs disrupted in disease?

This tests the accelerated-aging / genetic-convergence hypothesis for the IsoGraph switch layer,
projecting age-sensitive co-switching MODULES (defined out-of-cohort in the independent aging
fits) onto an independent schizophrenia case/control cohort (BrainSeq caudate_sczd, n~390) and
asking whether disease disrupts them in the direction genetics and aging predict. Structure is
learned independently; behaviour is measured where power lives (the disease cohort has numeric
age, 152 SCZD / 238 Control, and TOPMed genotypes).

Five layers (A is the headline):

  A. GWAS-directed switch concordance (genotype-anchored). For each SCZ-colocalized locus we
     regress the gene's switch axis on (i) the risk-allele dosage and (ii) diagnosis, *in the
     same disease cohort on the same switch score* — so the arbitrary switch-axis orientation
     cancels and sign(cis effect)==sign(Dx effect) is a clean test that the risk allele and the
     disease state push the isoform the SAME way. Aggregated by a binomial sign test; the GTEx
     coloc risk_qtl_effect is carried as an external cross-cohort direction QC. Abundance-
     conditioned (DTU-without-DGE).
  A3. Genotype x diagnosis. risk_dosage, Dx and risk_dosage:Dx on the driver switch and the
     module eigengene.
  B. Aging<->disease direction concordance. Per module gene, sign of the control-only switch-vs-
     age slope vs sign of the switch-vs-Dx effect: does SCZ recapitulate the age program gene by
     gene?
  C. Module preservation. Is the aging module's co-switch structure coherent in the disease
     switch network (median intramodular |r| vs a size-matched permutation null)?
  D. Age-state deviation. Fit each module eigengene's age trajectory in CONTROLS; test whether
     cases deviate from the control-expected state (accelerated-aging residual ~ Dx).
  Convergence. Do multiple SCZ-risk loci land in the same age-sensitive modules, and are their
     driver switches directionally disrupted the same way? **The answer is no**, once the
     background is the coloc-TESTED gene pool rather than all module genes: anchored modules are
     defined by MAGMA SCZ enrichment and so enter the coloc test preferentially, and conditioning
     on that removes the effect (48% vs a 45% tested-pool background, P=0.40; the superseded
     all-genes background gave 25% and P=0.004). See module_coloc_convergence.py, which repeats
     this for all five traits against a size-matched permutation null and finds no concentration
     anywhere. Layers B/C/D (disruption of the age-sensitive modules in disease) are unaffected.

Outputs under 08_integration/_m/scz_age_projection/ (parquet + Manubot SCZ_AGE_PROJECTION.md).
Genotype dosages are produced upstream by plink2 (eqtl env) in the SLURM wrapper.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import region_store, rel, stage_out
from isograph.io.artifacts import load_dataset_bundle
from isograph_benchmark.real_data.run_models import BRAINSEQ_COVARIATES
from isograph_benchmark.real_data.qtl_anchoring import _bare

SEED = 13
CASE, CTRL = "SCZD", "Control"
DX_COVS = ["Age", "Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
           "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"]
AGE_COVS = ["Sex", "RIN", "mapping_rate", "mito_rate",
            "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"]
_MIN_MODULE = 3
_OUT = stage_out("integration", "scz_age_projection")

# aging module sources (independent of the disease cohort's discovery), caudate-matched
_AGING_SOURCES = {
    "gtex_caudate_bg": {
        "modules": region_store("gtex", "caudate_basal_ganglia", "isograph_vae", "modules.parquet"),
        "magma_prefix": "gtex__caudate_basal_ganglia__",
        "rbp_region": "caudate_basal_ganglia",
    },
    "brainseq_caudate": {
        "modules": region_store("brainseq", "caudate", "isograph_vae", "modules.parquet"),
        "magma_prefix": "brainseq__caudate__",
        "rbp_region": "caudate",
    },
}
# candidate trans-regulators: per-module RBP regulons (rbp_regulon.py --scope combined)
_RBP_REGULON = stage_out("regulation", "rbp", "rbp_regulon_combined.parquet")
_DISEASE_FS = region_store("brainseq", "caudate_sczd", "isograph_vae", "feature_scores.parquet")
_DISEASE_BUNDLE = rel("inputs", "bundles", "brainseq_sczd", "caudate")
_MAGMA = stage_out("anchoring.gwas", "magma_results_combined.parquet")
_COLOC = stage_out("anchoring.coloc", "coloc_isoform_events_combined.parquet")


# --------------------------------------------------------------------------- #
# regression helpers
# --------------------------------------------------------------------------- #
def _design(frame: pd.DataFrame, term: str, covs: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Return (X, y) for `y ~ term + covs`; term is column 1 of X. Categorical covs dummied."""
    use = [term] + [c for c in covs if c in frame.columns]
    d = frame[["_y"] + use].replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 15:
        return np.empty((0, 0)), np.empty(0)
    y = d["_y"].to_numpy(float)
    t = pd.to_numeric(d[term], errors="coerce").to_numpy(float)[:, None] \
        if np.issubdtype(pd.to_numeric(d[term], errors="coerce").dtype, np.number) or d[term].dtype != object \
        else pd.get_dummies(d[term], drop_first=True).to_numpy(float)
    covcols = [c for c in covs if c in d.columns]
    cov = pd.get_dummies(d[covcols], drop_first=True).to_numpy(float) if covcols else np.empty((len(d), 0))
    X = np.hstack([np.ones((len(d), 1)), t, cov])
    return X, y


def _ols(X: np.ndarray, y: np.ndarray, idx: int = 1) -> tuple[float, float, int]:
    """(beta_idx, two-sided p, n). NaN on rank-deficiency."""
    if X.size == 0 or np.linalg.matrix_rank(X) <= idx or np.nanstd(y) == 0:
        return np.nan, np.nan, len(y)
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    dfres = max(len(y) - rank, 1)
    resid = y - X @ beta
    s2 = float(resid @ resid) / dfres
    se = np.sqrt(max(s2 * np.linalg.pinv(X.T @ X)[idx, idx], 1e-300))
    t = beta[idx] / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(t), dfres)) if np.isfinite(t) else np.nan
    return float(beta[idx]), p, len(y)


def _fit_term(y: np.ndarray, samp_frame: pd.DataFrame, term: str, term_vals: np.ndarray,
              covs: list[str]) -> tuple[float, float, int]:
    f = samp_frame.copy()
    f["_y"] = y
    f[term] = term_vals
    X, yy = _design(f, term, covs)
    return _ols(X, yy, idx=1)


# --------------------------------------------------------------------------- #
# data loading
# --------------------------------------------------------------------------- #
def _channel(fs: pd.DataFrame, ftype: str, samp: list[str]) -> pd.DataFrame:
    return fs[fs["feature_type"] == ftype].set_index("gene_id")[samp]


def _eigengene(channel: pd.DataFrame, genes) -> tuple[np.ndarray | None, int]:
    present = channel.index.intersection(pd.Index(genes).unique())
    if len(present) == 0:
        return None, 0
    return channel.loc[present].to_numpy(float).mean(axis=0), len(present)


def load_disease():
    fs = pd.read_parquet(_DISEASE_FS)
    st = load_dataset_bundle(_DISEASE_BUNDLE).sample_table.copy()
    sid = "sample_id" if "sample_id" in st.columns else st.columns[0]
    st = st.rename(columns={sid: "sample_id"})
    samp = [c for c in fs.columns if c in set(st["sample_id"])]
    st = st.set_index("sample_id").loc[samp].reset_index()
    SW = _channel(fs, "switch", samp)
    AB = _channel(fs, "abundance", samp)
    st["is_case"] = (st["Dx"] == CASE).astype(float)
    return SW, AB, st, samp


def load_scz_modules():
    """SCZ-GWAS-anchored aging caudate modules (MAGMA SCZ P<0.05) -> {(source, mid): genes}."""
    magma = pd.read_parquet(_MAGMA)
    scz = magma[(magma.trait == "SCZ") & (magma.backend == "isograph_vae")]
    p_of = dict(zip(scz["FULL_NAME"].astype(str), scz["P"]))
    fdr_of = dict(zip(scz["FULL_NAME"].astype(str), scz["FDR"]))
    out = {}
    for src, cfg in _AGING_SOURCES.items():
        mods = pd.read_parquet(cfg["modules"])
        for mid, g in mods.groupby("module_id"):
            fn = cfg["magma_prefix"] + str(mid)
            genes = list(pd.unique(g["gene_id"]))
            out[(src, str(mid))] = {"genes": genes, "scz_p": p_of.get(fn, np.nan),
                                    "scz_fdr": fdr_of.get(fn, np.nan),
                                    "anchored": p_of.get(fn, 1.0) < 0.05}
    return out


def load_genotypes(dosage_raw: Path, st: pd.DataFrame) -> pd.DataFrame | None:
    """risk-allele dosage matrix (sample_id x locus). Aligns counted->risk allele."""
    if not dosage_raw.exists():
        return None
    raw = pd.read_csv(dosage_raw, sep="\t")
    coloc = pd.read_parquet(_COLOC)
    coloc = coloc[coloc.trait.astype(str).str.lower().eq("scz")]
    risk_of = dict(zip(coloc.best_rsid.astype(str), coloc.risk_allele.astype(str)))
    br2sid = dict(zip(st["BrNum"].astype(str), st["sample_id"]))
    raw["sample_id"] = raw["FID"].astype(str).map(br2sid)
    raw = raw.dropna(subset=["sample_id"]).set_index("sample_id")
    dos = {}
    for col in raw.columns:
        if not col.startswith("chr"):
            continue
        toks = col.split("_")
        rs = next((t for t in toks if t.startswith("rs")), None)
        counted = toks[-1]
        if rs is None or rs not in risk_of:
            continue
        d = pd.to_numeric(raw[col], errors="coerce")
        if counted != risk_of[rs]:      # flip to count the risk allele
            d = 2 - d
        dos[rs] = d
    if not dos:
        return None
    return pd.DataFrame(dos)


# --------------------------------------------------------------------------- #
# A + A3: genotype-anchored GWAS-directed concordance and GxD
# --------------------------------------------------------------------------- #
def layer_A(SW, AB, st, geno, scz_mods):
    coloc = pd.read_parquet(_COLOC)
    coloc = coloc[coloc.trait.astype(str).str.lower().eq("scz")].copy()
    coloc["gene"] = coloc["gene"].astype(str)
    # feature_scores/modules gene_ids are versioned; coloc gene is bare -> resolve
    bare2ver = {}
    for gid in SW.index.unique():
        bare2ver.setdefault(_bare(pd.Series([gid]))[0], gid)
    # map bare gene -> its aging module (first anchored hit) for the module-score GxD outcome
    gene2mod = {}
    for (src, mid), rec in scz_mods.items():
        if rec["anchored"]:
            for g in rec["genes"]:
                gene2mod.setdefault(_bare(pd.Series([g]))[0], (src, mid, rec["genes"]))
    rows, gxd = [], []
    samp = list(SW.columns)
    stf = st.set_index("sample_id").loc[samp]
    for r in coloc.itertuples():
        bare = _bare(pd.Series([r.gene]))[0]
        gene = bare2ver.get(bare)
        if gene is None or gene not in SW.index:
            continue
        sw = SW.loc[[gene]].to_numpy(float).mean(axis=0)          # gene switch axis (fixed orient)
        ab = AB.loc[[gene]].to_numpy(float).mean(axis=0) if gene in AB.index else None
        covs_dx = DX_COVS + (["_abund"] if ab is not None else [])
        base = stf.copy()
        if ab is not None:
            base["_abund"] = ab
        # Dx effect on switch (abundance-conditioned)
        b_dx, p_dx, _ = _fit_term(sw, base, "is_case", stf["is_case"].to_numpy(float), covs_dx)
        rec = {"gene": gene, "gene_name": getattr(r, "gene_name", gene), "kind": r.kind,
               "go_invisible": bool(r.go_invisible), "risk_allele": r.risk_allele,
               "gtex_risk_qtl_effect": r.risk_qtl_effect, "clpp": r.clpp,
               "dx_beta_switch": b_dx, "dx_p": p_dx, "has_geno": False}
        if geno is not None and getattr(r, "best_rsid", None) in geno.columns:
            g = geno[r.best_rsid].reindex(samp).to_numpy(float)
            b_g, p_g, n_g = _fit_term(sw, base, "_dose", g, DX_COVS)   # cis effect in disease
            rec.update({"has_geno": True, "cis_beta_switch": b_g, "cis_p": p_g, "n_geno": n_g,
                        "concordant_cis_dx": (np.sign(b_g) == np.sign(b_dx))
                        if np.isfinite(b_g) and np.isfinite(b_dx) else np.nan,
                        "gtex_dir_matches_cis": (np.sign(r.risk_qtl_effect) == np.sign(b_g))
                        if np.isfinite(b_g) else np.nan})
            # A3: GxD on switch and on module eigengene
            for outcome_name, yv in _gxd_outcomes(sw, bare, gene2mod, SW):
                f = base.copy()
                f["_y"] = yv
                f["_dose"] = g
                f["_gxd"] = g * stf["is_case"].to_numpy(float)
                # explicit design: intercept, dose, dx, dose:dx, covs
                d = f[["_y", "_dose", "is_case", "_gxd"] + [c for c in DX_COVS if c in f.columns]]
                d = d.replace([np.inf, -np.inf], np.nan).dropna()
                if len(d) < 20:
                    continue
                yy = d["_y"].to_numpy(float)
                cov = pd.get_dummies(d[[c for c in DX_COVS if c in d.columns]], drop_first=True).to_numpy(float)
                Xg = np.hstack([np.ones((len(d), 1)), d[["_dose", "is_case", "_gxd"]].to_numpy(float), cov])
                bd, pd_, _ = _ols(Xg, yy, 1)
                bx, px, _ = _ols(Xg, yy, 2)
                bi, pi, n = _ols(Xg, yy, 3)
                gxd.append({"gene": gene, "gene_name": getattr(r, "gene_name", gene),
                            "outcome": outcome_name, "beta_dose": bd, "p_dose": pd_,
                            "beta_dx": bx, "p_dx": px, "beta_gxd": bi, "p_gxd": pi, "n": n})
        rows.append(rec)
    A = pd.DataFrame(rows)
    GXD = pd.DataFrame(gxd)
    return A, GXD


def _gxd_outcomes(sw, gene, gene2mod, SW):
    out = [("driver_switch", sw)]
    if gene in gene2mod:
        _, _, genes = gene2mod[gene]
        eig, n = _eigengene(SW, genes)
        if eig is not None:
            out.append(("module_eigengene", eig))
    return out


# --------------------------------------------------------------------------- #
# B: aging<->disease per-gene direction concordance
# --------------------------------------------------------------------------- #
def layer_B(SW, st, scz_mods):
    samp = list(SW.columns)
    stf = st.set_index("sample_id").loc[samp]
    ctrl = stf["Dx"] == CTRL
    age = pd.to_numeric(stf["Age"], errors="coerce").to_numpy(float)
    dx = stf["is_case"].to_numpy(float)
    rows = []
    for (src, mid), rec in scz_mods.items():
        if not rec["anchored"]:
            continue
        signs = []
        for gene in rec["genes"]:
            if gene not in SW.index:
                continue
            y = SW.loc[[gene]].to_numpy(float).mean(axis=0)
            # control-only age slope
            fa = stf.loc[ctrl].copy(); fa["_y"] = y[ctrl.to_numpy()]
            Xa, ya = _design(fa.assign(Age=age[ctrl.to_numpy()]), "Age", AGE_COVS)
            b_age, p_age, _ = _ols(Xa, ya, 1)
            # disease Dx effect
            b_dx, p_dx, _ = _fit_term(y, stf, "is_case", dx, DX_COVS)
            if np.isfinite(b_age) and np.isfinite(b_dx):
                signs.append(np.sign(b_age) == np.sign(b_dx))
        if len(signs) >= _MIN_MODULE:
            k = int(np.sum(signs)); n = len(signs)
            rows.append({"source": src, "module_id": mid, "n_genes": n, "n_concordant": k,
                         "concordance_rate": k / n, "scz_p": rec["scz_p"],
                         "binom_p": float(stats.binomtest(k, n, 0.5, "greater").pvalue)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# C: module preservation in disease switch network
# --------------------------------------------------------------------------- #
def layer_C(SW, scz_mods, n_perm=2000):
    rng = np.random.default_rng(SEED)
    pool = list(SW.index.unique())
    rows = []
    for (src, mid), rec in scz_mods.items():
        if not rec["anchored"]:
            continue
        present = SW.index.intersection(pd.Index(pd.unique(rec["genes"])))
        if len(present) < _MIN_MODULE:
            continue
        M = SW.loc[present].to_numpy(float)
        obs = _mean_abs_cor(M)
        null = np.array([_mean_abs_cor(SW.loc[rng.choice(pool, len(present), replace=False)].to_numpy(float))
                         for _ in range(n_perm)])
        null = null[np.isfinite(null)]
        z = (obs - null.mean()) / (null.std() + 1e-12) if len(null) else np.nan
        emp = float((null >= obs).mean()) if len(null) else np.nan
        rows.append({"source": src, "module_id": mid, "n_present": int(len(present)),
                     "mean_abs_cor": obs, "null_mean": float(null.mean()) if len(null) else np.nan,
                     "preservation_z": float(z), "emp_p": emp, "scz_p": rec["scz_p"]})
    return pd.DataFrame(rows)


def _mean_abs_cor(M):
    if M.shape[0] < 2:
        return np.nan
    C = np.corrcoef(M)
    iu = np.triu_indices_from(C, k=1)
    v = C[iu]
    v = v[np.isfinite(v)]
    return float(np.abs(v).mean()) if len(v) else np.nan


# --------------------------------------------------------------------------- #
# D: age-state deviation (accelerated aging)
# --------------------------------------------------------------------------- #
def layer_D(SW, st, scz_mods):
    samp = list(SW.columns)
    stf = st.set_index("sample_id").loc[samp]
    ctrl = (stf["Dx"] == CTRL).to_numpy()
    age = pd.to_numeric(stf["Age"], errors="coerce").to_numpy(float)
    dx = stf["is_case"].to_numpy(float)
    rows = []
    for (src, mid), rec in scz_mods.items():
        if not rec["anchored"]:
            continue
        eig, n = _eigengene(SW, rec["genes"])
        if eig is None or n < _MIN_MODULE:
            continue
        # fit eigengene ~ age (linear) in CONTROLS, predict expected, residual for all
        fc = pd.DataFrame({"_y": eig[ctrl], "Age": age[ctrl]}).dropna()
        if len(fc) < 15 or fc["Age"].std() == 0:
            continue
        b = np.polyfit(fc["Age"], fc["_y"], 1)
        resid = eig - np.polyval(b, age)
        # do cases deviate from the control-expected state? residual ~ Dx (+covs)
        b_dx, p_dx, nn = _fit_term(resid, stf, "is_case", dx, DX_COVS)
        rows.append({"source": src, "module_id": mid, "n_present": int(n),
                     "age_slope_control": float(b[0]), "case_deviation_beta": b_dx,
                     "case_deviation_p": p_dx, "n": nn, "scz_p": rec["scz_p"]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# convergence
# --------------------------------------------------------------------------- #
def layer_convergence(A, scz_mods):
    """How many SCZ-risk loci land in each anchored age-sensitive module, are their driver
    switches disrupted the same way, and are coloc loci enriched in the anchored (age-sensitive)
    modules vs all aging caudate modules (hypergeometric)?"""
    rows = []
    hits = A.dropna(subset=["dx_beta_switch"])
    coloc_genes = set(hits["gene"])
    # Hypergeometric: are coloc genes concentrated in anchored (age-sensitive) modules?
    #
    # THE DENOMINATOR MATTERS AND AN EARLIER VERSION OF THIS GOT IT WRONG. Using ALL
    # module genes as the background gave 15/31 = 48% vs 25%, P=0.004 — the published
    # convergence headline. But a gene can only appear in `coloc_genes` if it was
    # COLOC-TESTED, i.e. if it sat under a SCZ GWAS peak with a QTL credible set. Anchored
    # modules are defined by MAGMA SCZ P<0.05, so their genes preferentially sit under SCZ
    # peaks and preferentially enter the tested pool: the tested pool is already ~45%
    # anchored, and against it 48% is null (P=0.40). MAGMA anchoring and coloc testing
    # select on the same GWAS signal, so the all-genes background measures that shared
    # ascertainment, not convergence.
    #
    # The tested pool is therefore the background, and the all-genes figure is retained
    # only so the inflation stays visible. See module_coloc_convergence.py, which runs
    # this comparison for all five traits with a size-matched permutation null.
    anch_genes, all_genes = set(), set()
    for (src, mid), rec in scz_mods.items():
        all_genes |= set(rec["genes"])
        if rec["anchored"]:
            anch_genes |= set(rec["genes"])
    tested = set(A["gene"]) & all_genes          # every gene that could have been a hit
    tested_anch = tested & anch_genes
    N, K = len(tested), len(tested_anch)
    drawn = coloc_genes & tested
    x = len(coloc_genes & tested_anch)
    enrich_p = float(stats.hypergeom.sf(x - 1, N, K, len(drawn))) if (len(drawn) and N) else np.nan
    N_all, K_all = len(all_genes), len(anch_genes)
    enrich = {"n_coloc_in_pool": len(drawn), "n_coloc_in_anchored": x,
              "n_tested_pool": N, "frac_anchored_pool": K / N if N else np.nan,
              "frac_coloc_anchored": x / len(drawn) if drawn else np.nan,
              "hyperg_p": enrich_p,
              # superseded all-module-genes background, kept to show the ascertainment
              "frac_anchored_allgenes": K_all / N_all if N_all else np.nan,
              "hyperg_p_allgenes_denom": float(
                  stats.hypergeom.sf(x - 1, N_all, K_all, len(drawn))) if len(drawn) else np.nan}
    for (src, mid), rec in scz_mods.items():
        if not rec["anchored"]:
            continue
        genes = set(rec["genes"])
        sub = hits[hits["gene"].isin(genes)]
        if sub.empty:
            continue
        signs = np.sign(sub["dx_beta_switch"].to_numpy())
        signs = signs[np.isfinite(signs)]
        rows.append({"source": src, "module_id": mid, "n_coloc_loci": int(sub["gene"].nunique()),
                     "n_go_invisible": int(sub["go_invisible"].sum()),
                     "frac_same_dx_direction": float(max((signs > 0).mean(), (signs < 0).mean()))
                     if len(signs) else np.nan, "scz_p": rec["scz_p"]})
    conv = pd.DataFrame(rows).sort_values("n_coloc_loci", ascending=False)
    conv.attrs["enrich"] = enrich
    return conv


# --------------------------------------------------------------------------- #
def layer_mechanism(scz_mods, conv, rbp_q=0.05):
    """Attach candidate trans-regulators to each anchored age-sensitive module by joining the
    per-module RBP regulons (rbp_regulon --scope combined). Turns 'convergence exists' into a
    named-splicing-factor hypothesis per module. Returns (per-hit table, {(src,mid): 'RBP(q); ...'})."""
    if not _RBP_REGULON.exists():
        return pd.DataFrame(), {}
    reg = pd.read_parquet(_RBP_REGULON)
    region_of = {src: cfg["rbp_region"] for src, cfg in _AGING_SOURCES.items()}
    rows = []
    for (src, mid), rec in scz_mods.items():
        if not rec["anchored"]:
            continue
        sub = reg[(reg.region == region_of.get(src)) & (reg.module_id.astype(str) == str(mid))
                  & (reg.q < rbp_q)].sort_values("q")
        for r in sub.itertuples():
            rows.append({"source": src, "module_id": mid, "region": region_of.get(src),
                         "rbp": r.rbp, "enrichment": r.enrichment, "q": r.q,
                         "go_invisible": bool(r.go_invisible),
                         "pool_source": getattr(r, "pool_source", "switch_genes"),
                         "scz_p": rec["scz_p"]})
    mech = pd.DataFrame(rows)
    top = {}
    if not mech.empty:
        for (src, mid), g in mech.groupby(["source", "module_id"]):
            g = g.sort_values("q").head(5)
            top[(src, mid)] = "; ".join(f"{r.rbp}(q={r.q:.2g})" for r in g.itertuples())
    return mech, top


def _fdr(df, pcol, into):
    if not df.empty and df[pcol].notna().any():
        m = df[pcol].notna()
        df.loc[m, into] = stats.false_discovery_control(df.loc[m, pcol], method="bh")
    return df


def _write_report(A, GXD, B, C, D, conv, mech=None):
    L = ["# Age-sensitive isoform-switch programs in schizophrenia",
         "",
         "Age-sensitive co-switching modules are defined out-of-cohort in the independent aging "
         "caudate fits (GTEx caudate basal ganglia + BrainSeq caudate) and restricted to those "
         "enriched for schizophrenia GWAS (MAGMA SCZ P<0.05). Their **behaviour** is tested in an "
         "independent SCZ case/control cohort (BrainSeq caudate_sczd) with numeric age and TOPMed "
         "genotypes. Effects are abundance-conditioned (DTU-without-DGE).", ""]
    def _boolsum(s):
        v = pd.to_numeric(s.map({True: 1, False: 0}), errors="coerce").dropna()
        return int(v.sum()), int(len(v))
    nB = int((B.binom_p < 0.05).sum()) if B is not None and not B.empty else 0
    nD = int((D.case_deviation_p < 0.05).sum()) if D is not None and not D.empty else 0
    nC = int((C.emp_p < 0.05).sum()) if C is not None and not C.empty else 0

    # ---- HEADLINE: genetic convergence + module-level disruption (best-powered result) --------
    if conv is not None and not conv.empty:
        en = conv.attrs.get("enrich", {})
        L += ["## SCZ-risk loci and age-sensitive switch programs", ""]
        if en and np.isfinite(en.get("hyperg_p", np.nan)):
            sig = en["hyperg_p"] < 0.05
            L += [
                ("Schizophrenia-colocalized switch genes are **concentrated in the "
                 if sig else
                 "Schizophrenia-colocalized switch genes are **not concentrated in the ") +
                f"age-sensitive (SCZ-GWAS-enriched) modules**: "
                f"{en['n_coloc_in_anchored']}/{en['n_coloc_in_pool']} "
                f"({en['frac_coloc_anchored']:.0%}) of coloc genes fall in anchored modules, against a "
                f"{en['frac_anchored_pool']:.0%} background among the genes that were coloc-TESTED "
                f"(hypergeometric P={en['hyperg_p']:.3g}).",
                "",
                "> **The background is the whole result, and an earlier version of this report used the "
                "wrong one.** Against *all* module genes the background is only "
                f"{en.get('frac_anchored_allgenes', float('nan')):.0%} and the same counts give "
                f"P={en.get('hyperg_p_allgenes_denom', float('nan')):.3g} — the previously reported "
                "convergence headline. That comparison is confounded by ascertainment: a gene can only "
                "colocalize if it sat under a SCZ GWAS peak with a QTL credible set, and anchored modules "
                "are *defined* by MAGMA SCZ enrichment, so their genes enter the tested pool "
                "preferentially. Conditioning on what could have been a hit removes the effect. "
                "`module_coloc_convergence.py` repeats this for all five traits with a size-matched "
                "permutation null and finds no concentration anywhere (P = 0.19–1.00).",
                "",
                "The modules carrying the most colocalized loci are listed below; with these counts the "
                "per-module numbers are descriptive, not evidence of convergence.", ""]
        L += ["| source | module | # coloc genes | # GO-invisible | max same-dir frac | SCZ MAGMA P |",
              "|--------|--------|-----------|----------------|-------------------|-------------|"]
        for r in conv.head(6).itertuples():
            L.append(f"| {r.source} | {r.module_id} | {r.n_coloc_loci} | {r.n_go_invisible} | "
                     f"{r.frac_same_dx_direction:.2f} | {r.scz_p:.2g} |")
        L += ["",
              "Independently of that null, the age-sensitive modules are directionally disrupted in disease: they recapitulate "
              f"the aging switch direction gene-by-gene in **{nB}/{len(B) if B is not None else 0}** "
              f"modules (B), show case deviation from the control age trajectory in **{nD}/"
              f"{len(D) if D is not None else 0}** modules (D), and stay co-switch-coherent in disease "
              f"in **{nC}/{len(C) if C is not None else 0}** modules (C).", ""]
        # candidate trans-regulators — the mechanism hypothesis per convergent module
        if mech is not None and not mech.empty:
            L += ["### Candidate trans-regulators (mechanism)", "",
                  "Each listed module's members are tested for shared RBP binding-site switching "
                  "(rbp_regulon --scope combined; mature+intronic motif scan). Significant RBPs "
                  "(q<0.05) are candidate trans regulators coordinating the co-switch program — a "
                  "named, testable hypothesis. These modules are NOT established as points of "
                  "SCZ-risk convergence (see the background caveat above); the regulators are "
                  "candidates for the modules' own co-switching, not for a convergence effect:", "",
                  "| source | module | # SCZ loci | top candidate RBP regulators (q) |",
                  "|--------|--------|-----------|----------------------------------|"]
            convtop = conv.head(6) if conv is not None else pd.DataFrame()
            for r in convtop.itertuples():
                reg = getattr(r, "candidate_rbp_regulators", "") or "—"
                L.append(f"| {r.source} | {r.module_id} | {r.n_coloc_loci} | {reg} |")
            L += ["", "_Motif-based candidate regulation (predicted binding-site gain/loss between "
                  "switch isoforms), not experimental validation._", ""]

    # ---- B / D / C detail ---------------------------------------------------------------------
    if B is not None and not B.empty:
        top = B.sort_values("binom_p").head(5)
        L += ["## B. Aging↔disease direction concordance (accelerated-aging recapitulation)", "",
              f"Per module gene, sign of the control-only switch-vs-age slope vs the switch-vs-Dx "
              f"effect. **{nB}/{len(B)}** anchored modules show above-chance concordance (binomial "
              "P<0.05) — SCZ recapitulates the age switch program gene by gene. Top:", "",
              "| source | module | concordant/n | rate | binom P |",
              "|--------|--------|--------------|------|---------|"]
        for r in top.itertuples():
            L.append(f"| {r.source} | {r.module_id} | {r.n_concordant}/{r.n_genes} | "
                     f"{r.concordance_rate:.2f} | {r.binom_p:.3g} |")
        L.append("")
    if D is not None and not D.empty:
        L += ["## D. Age-state deviation", "",
              f"Each module eigengene's age trajectory is fit in controls; **{nD}/{len(D)}** anchored "
              "modules show cases deviating from the control-expected state (residual ~ Dx, P<0.05).", ""]
    if C is not None and not C.empty:
        L += ["## C. Module preservation", "",
              f"**{nC}/{len(C)}** anchored aging modules keep coherent co-switch structure in the disease "
              "switch network (median intramodular |r| above a size-matched permutation null) — a "
              "prerequisite the projections rely on.", ""]

    # ---- Supplementary: single-locus resolution (honest, underpowered) ------------------------
    L += ["---", "",
          "## Supplementary — single-locus genotype resolution", "",
          "Per-locus genotype tests are reported for completeness; single-locus QTL power in the "
          "disease cohort is far below the module-level analyses above, so these do not carry the "
          "narrative.", ""]
    if not A.empty:
        g = A[A.has_geno]
        k, n = _boolsum(g["concordant_cis_dx"])
        L += ["### S1. Genotype-anchored cis↔diagnosis switch concordance", ""]
        if n:
            bp = stats.binomtest(k, n, 0.5, "greater").pvalue
            gk, gn = _boolsum(g["gtex_dir_matches_cis"])
            L += [f"For **{n}** genotyped SCZ-colocalized loci, the risk-allele cis effect and the "
                  f"diagnosis effect on the same disease switch axis agree in sign in **{k}/{n}** "
                  f"({k/n:.0%}; binomial P={bp:.3g}) — **null at single-locus resolution**. The in-cohort "
                  f"cis effects are weak (median cis P={g['cis_p'].median():.2f}); QTL power in n~{n} is "
                  "far below GTEx discovery (GTEx cis-direction match "
                  f"{gk}/{gn}), which is why the directional signal is resolved at the module level "
                  "(headline convergence + B/D) rather than per locus.", ""]
        sig = A[(A.dx_p < 0.05)]
        L += [f"- coloc loci with a nominal disease switch shift (Dx P<0.05, abundance-conditioned): "
              f"**{len(sig)}/{len(A)}** ({int((A.go_invisible & (A.dx_p<0.05)).sum())} GO-invisible)", ""]
    if GXD is not None and not GXD.empty:
        gsig = GXD[(GXD.p_gxd < 0.05)]
        L += ["### S2. Genotype × diagnosis", "",
              f"- genotype×diagnosis interactions at P<0.05: **{len(gsig)}** of {len(GXD)} "
              "locus×outcome tests (exploratory; single-locus power caveat applies).", ""]
    (_OUT / "SCZ_AGE_PROJECTION.md").write_text("\n".join(L) + "\n")


def run(dosage_raw: Path, n_perm: int = 2000) -> None:
    _OUT.mkdir(parents=True, exist_ok=True)
    SW, AB, st, samp = load_disease()
    print(f"disease cohort: {len(samp)} samples "
          f"({int((st['Dx']==CASE).sum())} {CASE} / {int((st['Dx']==CTRL).sum())} {CTRL}); "
          f"switch genes {SW.shape[0]}")
    scz_mods = load_scz_modules()
    n_anch = sum(v["anchored"] for v in scz_mods.values())
    print(f"aging caudate modules: {len(scz_mods)} ({n_anch} SCZ-anchored, MAGMA P<0.05)")
    geno = load_genotypes(dosage_raw, st)
    print(f"genotypes: {'none' if geno is None else f'{geno.shape[1]} loci x {geno.shape[0]} samples'}")

    A, GXD = layer_A(SW, AB, st, geno, scz_mods)
    A = _fdr(A, "dx_p", "dx_fdr")
    if not GXD.empty:
        GXD = _fdr(GXD, "p_gxd", "fdr_gxd")
    B = layer_B(SW, st, scz_mods)
    B = _fdr(B, "binom_p", "fdr")
    C = layer_C(SW, scz_mods, n_perm=n_perm)
    C = _fdr(C, "emp_p", "fdr")
    D = layer_D(SW, st, scz_mods)
    D = _fdr(D, "case_deviation_p", "fdr")
    conv = layer_convergence(A, scz_mods)
    mech, top_rbp = layer_mechanism(scz_mods, conv)
    if not conv.empty:
        conv["candidate_rbp_regulators"] = conv.apply(
            lambda r: top_rbp.get((r["source"], r["module_id"]), ""), axis=1)
    print(f"mechanism: {0 if mech.empty else len(mech)} module x RBP regulon hits "
          f"across {0 if mech.empty else mech.groupby(['source','module_id']).ngroups} anchored modules")

    for name, df in [("A_gwas_directed_concordance", A), ("A3_gxd_models", GXD),
                     ("B_direction_concordance", B), ("C_module_preservation", C),
                     ("D_age_deviation", D), ("convergence", conv), ("mechanism_rbp", mech)]:
        df.to_parquet(_OUT / f"{name}.parquet", index=False)
        df.to_csv(_OUT / f"{name}.tsv", sep="\t", index=False)
    _write_report(A, GXD, B, C, D, conv, mech)
    print(f"wrote 7 tables + SCZ_AGE_PROJECTION.md to {_OUT}")


def main() -> None:
    p = argparse.ArgumentParser(description="SCZ-risk convergence on age-sensitive switch modules.")
    p.add_argument("--dosage", type=Path,
                   default=_OUT / "genotypes" / "scz_loci_dosage.raw",
                   help="plink2 --export A dosage table for SCZ coloc lead loci (from wrapper).")
    p.add_argument("--n-perm", type=int, default=2000)
    args = p.parse_args()
    run(args.dosage, args.n_perm)


if __name__ == "__main__":
    main()
