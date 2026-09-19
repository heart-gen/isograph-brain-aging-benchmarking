#!/usr/bin/env Rscript
## Locus-specific LD robustness audit for a single fine-mapped GWAS locus.
##
## WHY THIS EXISTS
## ---------------
## The uniform primary coloc grid uses MAX_SNPS = 12,000 in 03c.coloc_gwas_susie.R, which
## drops 61 of 579 loci before any fit. Raising the guard recovers them, but the recovered
## loci were EMPIRICALLY ENRICHED for greater GWAS-reference-LD inconsistency: in
## aging__ad all 11 recovered loci had s_rss >= 0.310 against a median of 0.255 over all
## 49 fitted, worst 0.739. That is a correlation between a compute threshold and
## statistical difficulty. It is NOT a demonstration of why the two travel together, and
## nothing here should be written as though large loci are ill-conditioned BECAUSE of
## long-range LD -- this audit does not test that and cannot support it.
##
## The consequence is a policy, not a global re-run: the recovered loci are a scoped
## SENSITIVITY arm, and any one of them that is going to carry a biological claim must
## clear a locus-specific audit first. This script is that audit. PICALM
## (aging__ad / locus60_chr11 / ENSG00000073921) is the locus it was written for, but
## nothing here is PICALM-specific -- it takes the locus as an argument so that any other
## recovered locus faces the identical bar.
##
## THE FIVE CHECKS
## ---------------
##   1 convergence   susie_rss converged, in how many iterations, with how many credible
##                   sets. A fit that hit maxit is not evidence of anything.
##   2 purity        per credible set: size, purity (min/mean/median |r| within the set),
##                   lead variant and PIP, and the max |r| BETWEEN credible-set leads.
##                   Two "independent" signals in near-perfect LD are the signature of a
##                   fit splitting one signal under LD mismatch, not of two signals.
##   3 boundaries    refit on windows trimmed by 100/250/500 kb at each edge. A real
##                   signal should not depend on where an arbitrary window was cut. The
##                   comparison is the credible set's lead variant and its membership
##                   Jaccard against the primary fit.
##   4 coloc         rerun coloc.susie against the gene's sQTL introns for every variant
##                   fit above. What matters for the narrative is not whether a credible
##                   set exists but whether the COLOCALIZATION conclusion moves.
##   5 ld_mismatch   the LD-mismatch-aware pass. `kriging_rss` gives per-variant logLR for
##                   z-scores inconsistent with the reference panel; the fit is repeated
##                   after dropping the flagged variants, and again with
##                   `estimate_residual_variance = TRUE`, which is susieR's own
##                   accommodation for exactly this failure. Three answers that agree are
##                   worth far more than one s_rss.
##
## Usage: Rscript 08d.locus_ld_robustness.R <analysis> <LOCUS_ID> <gene_bare> <tissue[,tissue...]>
## e.g.   Rscript 08d.locus_ld_robustness.R aging__ad locus60_chr11 ENSG00000073921 Brain_Cortex
## Output: 05_genetic_anchoring/_m/locus_ld_robustness/<analysis>__<LOCUS_ID>/
##         (per locus; every tissue lands in the same coloc_stability.parquet, keyed by `tissue`)
suppressPackageStartupMessages({
    library(data.table); library(arrow); library(susieR); library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("need <analysis> <LOCUS_ID> <gene_bare> <tissue[,tissue...]>")
ANALYSIS <- args[1]; LID <- args[2]; GENE <- args[3]
TISSUES  <- unique(trimws(strsplit(args[4], ",", fixed = TRUE)[[1]]))
TISSUES  <- TISSUES[nzchar(TISSUES)]
if (!length(TISSUES)) stop("no tissue given")

ROOT      <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
GTEX      <- "/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11"
PANEL_DIR <- "/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink"
BASE   <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc_signal_susie")
CDIR   <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc", ANALYSIS, "susie")
BRIDGE <- file.path(ROOT, "inputs", "raw", "gtex_v11", "variant_bridge")
OUTD   <- file.path(ROOT, "05_genetic_anchoring", "_m", "locus_ld_robustness",
                    sprintf("%s__%s", ANALYSIS, LID))
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

## Identical to stages A and B. The audit must not be more permissive than the pipeline
## it is auditing, or a locus could clear here and fail there.
L_MAX      <- 10L
MIN_SHARED <- 50L
## The same sweep stage B runs, not just the primary: check 4 asks whether the coloc
## CONCLUSION is stable, and a conclusion that holds only at one prior is not stable.
P12_SWEEP  <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
P12        <- 1e-5                       # the primary prior
TRIMS      <- c(0L, 100000L, 250000L, 500000L)
## susieR flags a variant when the observed z is far from the value the reference LD
## predicts from its neighbours. 2 is susieR's own suggested logLR cut in the LD-mismatch
## vignette; it is a screen, not a test, so the count of flagged variants is reported
## alongside every result that depends on it.
KRIGING_LOGLR <- 2

is_ambig <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

read_ld <- function(lid) {
    vf <- file.path(CDIR, paste0(lid, ".unphased.vcor1.bin.vars"))
    bf <- file.path(CDIR, paste0(lid, ".unphased.vcor1.bin"))
    if (!file.exists(vf) || !file.exists(bf)) return(NULL)
    v <- readLines(vf); p <- length(v)
    R <- matrix(readBin(bf, "double", n = p * p, size = 4), p, p)
    dimnames(R) <- list(v, v); R
}

## ------------------------------------------------------------------ GWAS side
loci <- fread(file.path(CDIR, "loci_testable.tsv"))
Lrow <- loci[LOCUS_ID == LID]
if (!nrow(Lrow)) stop("locus not in loci_testable.tsv: ", LID)
CH <- Lrow$chr[1]

R0 <- read_ld(LID); if (is.null(R0)) stop("no LD matrix for ", LID)
gwas <- fread(file.path(CDIR, paste0(LID, ".gwas.tsv")))
excl_f <- file.path(CDIR, "exclude_regions.tsv")
excl   <- if (file.exists(excl_f)) fread(excl_f) else
          data.table(chr = 6L, start = 25e6, stop = 34e6)
for (j in seq_len(nrow(excl)))
    if (CH == excl$chr[j]) gwas <- gwas[!(pos >= excl$start[j] & pos <= excl$stop[j])]

bim <- fread(file.path(PANEL_DIR, sprintf("1000G.EUR.QC.%d.bim", CH)), header = FALSE,
             col.names = c("bchr","rsid","cm","bpos","A1","A2"))
bim <- bim[rsid %in% rownames(R0), .(rsid, REF = A2, ALT = A1)]
dt  <- merge(gwas, bim, by = "rsid")
dt  <- dt[((a1==REF & a2==ALT) | (a1==ALT & a2==REF)) & !is_ambig(a1, a2)]
dt[, z_ref := z * fifelse(a1==REF, 1, -1)]
dt  <- dt[rsid %in% rownames(R0)][!duplicated(rsid)]
NN  <- max(1000, round(Lrow$min_neff[1]))
message(sprintf("[%s / %s] %d usable SNPs, n_eff=%d, span %s:%d-%d",
                ANALYSIS, LID, nrow(dt), NN, CH, min(dt$pos), max(dt$pos)))

fit_window <- function(keep_rsid) {
    d <- dt[rsid %in% keep_rsid]
    if (nrow(d) < 50) return(NULL)
    Rk <- R0[d$rsid, d$rsid, drop = FALSE]
    f <- tryCatch(susie_rss(z = d$z_ref, R = Rk, n = NN, L = L_MAX,
                            estimate_residual_variance = FALSE),
                  error = function(e) { message("  susie failed: ", conditionMessage(e)); NULL })
    if (is.null(f)) return(NULL)
    list(fit = coloc::annotate_susie(f, d$rsid, Rk), d = d, R = Rk)
}

## ---------------------------------------------------------- 1 + 2: primary fit
P <- fit_window(dt$rsid)
if (is.null(P)) stop("primary fit failed for ", LID)

cs_table <- function(f, Rk, label) {
    if (!length(f$sets$cs)) return(data.table(variant = label, cs = NA_integer_))
    rows <- lapply(seq_along(f$sets$cs), function(k) {
        ix <- f$sets$cs[[k]]
        nm <- names(ix); if (is.null(nm)) nm <- rownames(Rk)[ix]
        sub <- abs(Rk[nm, nm, drop = FALSE]); diag(sub) <- NA_real_
        lead <- nm[which.max(f$pip[nm])]
        ## A one-variant set has no off-diagonal pairs: min() of nothing is Inf and mean()
        ## NaN. Its purity is 1 by definition, which is also what susieR reports.
        single <- length(nm) == 1L
        data.table(variant = label, cs = f$sets$cs_index[k], cs_size = length(ix),
                   purity_min = if (single) 1 else min(sub, na.rm = TRUE),
                   purity_mean = if (single) 1 else mean(sub, na.rm = TRUE),
                   purity_median = if (single) 1 else median(sub, na.rm = TRUE),
                   lead_rsid = lead, lead_pip = unname(f$pip[lead]),
                   members = paste(nm, collapse = ";"))
    })
    rbindlist(rows, fill = TRUE)
}

## Max |r| between the lead variants of DIFFERENT credible sets. Near 1 means the fit has
## split one signal in two rather than found two, which is the classic LD-mismatch
## artefact and is invisible in a per-set purity number.
cross_cs_r <- function(f, Rk) {
    if (length(f$sets$cs) < 2) return(NA_real_)
    leads <- vapply(f$sets$cs, function(ix) {
        nm <- names(ix); if (is.null(nm)) nm <- rownames(Rk)[ix]
        nm[which.max(f$pip[nm])]
    }, character(1))
    m <- abs(Rk[leads, leads, drop = FALSE]); diag(m) <- NA_real_
    max(m, na.rm = TRUE)
}

variants <- list(primary = P)
conv <- list(data.table(variant = "primary", converged = isTRUE(P$fit$converged),
                        niter = P$fit$niter, n_snp = nrow(P$d),
                        n_cs = length(P$fit$sets$cs),
                        cross_cs_max_abs_r = cross_cs_r(P$fit, P$R)))

## ------------------------------------------------------------- 3: boundary trims
for (tr in TRIMS[TRIMS > 0]) {
    lo <- min(dt$pos) + tr; hi <- max(dt$pos) - tr
    keep <- dt[pos >= lo & pos <= hi, rsid]
    lab <- sprintf("trim%dkb", tr / 1000L)
    W <- fit_window(keep)
    if (is.null(W)) { message("  ", lab, ": fit failed / too few SNPs"); next }
    variants[[lab]] <- W
    conv[[length(conv)+1]] <- data.table(
        variant = lab, converged = isTRUE(W$fit$converged), niter = W$fit$niter,
        n_snp = nrow(W$d), n_cs = length(W$fit$sets$cs),
        cross_cs_max_abs_r = cross_cs_r(W$fit, W$R))
    message(sprintf("  %-12s %d SNPs, %d credible sets", lab, nrow(W$d),
                    length(W$fit$sets$cs)))
}

## ------------------------------------------------- 5: LD-mismatch-aware variants
s_rss <- tryCatch(susieR::estimate_s_rss(z = P$d$z_ref, R = P$R, n = NN),
                  error = function(e) NA_real_)
message(sprintf("  s_rss = %.4f", s_rss))

kr <- tryCatch(susieR::kriging_rss(z = P$d$z_ref, R = P$R, n = NN),
               error = function(e) { message("  kriging_rss failed: ",
                                             conditionMessage(e)); NULL })
n_flagged <- NA_integer_
if (!is.null(kr)) {
    cond <- as.data.table(kr$conditional_dist)
    cond[, rsid := P$d$rsid]
    flag <- cond[logLR > KRIGING_LOGLR, rsid]
    n_flagged <- length(flag)
    message(sprintf("  kriging_rss: %d / %d variants flagged at logLR > %g",
                    n_flagged, nrow(cond), KRIGING_LOGLR))
    fwrite(cond, file.path(OUTD, "kriging_conditional_dist.tsv.gz"), sep = "\t")
    if (n_flagged > 0 && n_flagged < nrow(cond) - 50) {
        W <- fit_window(setdiff(P$d$rsid, flag))
        if (!is.null(W)) {
            variants[["drop_kriging_outliers"]] <- W
            conv[[length(conv)+1]] <- data.table(
                variant = "drop_kriging_outliers", converged = isTRUE(W$fit$converged),
                niter = W$fit$niter, n_snp = nrow(W$d),
                n_cs = length(W$fit$sets$cs),
                cross_cs_max_abs_r = cross_cs_r(W$fit, W$R))
        }
    }
}

## susieR's own LD-mismatch accommodation: let the residual variance absorb the
## inconsistency instead of forcing it to 1. If a credible set survives BOTH this and the
## outlier drop, the signal is not an artefact of the reference panel.
##
## susieR prints a HINT here that `estimate_residual_variance = TRUE` assumes in-sample
## LD. That is not an oversight: we are deliberately using it off-label as a mismatch
## probe, because the residual variance is the parameter that absorbs z-vs-LD
## inconsistency. Treat it as a THIRD OPINION that should agree with the other two, never
## as a better-specified fit than the primary.
fe <- tryCatch(susie_rss(z = P$d$z_ref, R = P$R, n = NN, L = L_MAX,
                         estimate_residual_variance = TRUE),
               error = function(e) { message("  erv fit failed: ",
                                             conditionMessage(e)); NULL })
if (!is.null(fe)) {
    variants[["estimate_resid_var"]] <- list(
        fit = coloc::annotate_susie(fe, P$d$rsid, P$R), d = P$d, R = P$R)
    conv[[length(conv)+1]] <- data.table(
        variant = "estimate_resid_var", converged = isTRUE(fe$converged),
        niter = fe$niter, n_snp = nrow(P$d), n_cs = length(fe$sets$cs),
        cross_cs_max_abs_r = cross_cs_r(fe, P$R))
}

convergence <- rbindlist(conv, fill = TRUE)
convergence[, `:=`(analysis = ANALYSIS, LOCUS_ID = LID, s_rss = s_rss,
                   n_kriging_flagged = n_flagged)]

purity <- rbindlist(lapply(names(variants),
                           function(v) cs_table(variants[[v]]$fit, variants[[v]]$R, v)),
                    fill = TRUE)

## Membership agreement with the primary fit: does a perturbation find the SAME set?
prim_sets <- if (length(P$fit$sets$cs))
    lapply(P$fit$sets$cs, function(ix) names(ix)) else list()
purity[, jaccard_vs_primary := NA_real_]
if (length(prim_sets) && "members" %in% names(purity)) {
    for (i in seq_len(nrow(purity))) {
        m <- strsplit(purity$members[i], ";", fixed = TRUE)[[1]]
        if (!length(m) || is.na(purity$cs[i])) next
        purity$jaccard_vs_primary[i] <- max(vapply(prim_sets, function(ps)
            length(intersect(ps, m)) / length(union(ps, m)), numeric(1)))
    }
}

## ------------------------------------------------------------------- 4: coloc
## One or more tissues. Checks 1-3 and 5 concern the GWAS fit alone, so a locus that
## colocalizes in several tissues is audited in ONE run instead of re-fitting the GWAS side
## per tissue -- and, since the output dir is per locus, instead of a second tissue's run
## overwriting the first.
br   <- as.data.table(read_parquet(file.path(BRIDGE, sprintf("chr%d.parquet", CH))))
bim2 <- fread(file.path(PANEL_DIR, sprintf("1000G.EUR.QC.%d.bim", CH)),
              header = FALSE, col.names = c("bchr","rsid","cm","bpos","A1","A2"))
bim2 <- bim2[, .(rsid, REF = A2, ALT = A1)]

coloc_rows <- list()
for (TISSUE in TISSUES) {
    sf <- file.path(GTEX, "GTEx_Analysis_v11_sQTL_all_associations",
                    sprintf("%s.v11.cis_sqtl.allpairs.chr%d.parquet", TISSUE, CH))
    COVF <- file.path(GTEX, "GTEx_Analysis_v11_eQTL_covariates",
                      sprintf("%s.v11.covariates.txt", TISSUE))
    if (!file.exists(sf) || !file.exists(COVF)) {
        message(sprintf("  %s: sQTL all-pairs or covariates not on disk -- skipped", TISSUE))
        next
    }
    N_QTL <- length(strsplit(readLines(COVF, n = 1L), "\t")[[1]]) - 1L

    sd_ <- open_dataset(sf)
    pid <- as.data.table(sd_ |> dplyr::distinct(phenotype_id) |> dplyr::collect())
    pid[, g := sub("\\..*$", "", sub("^.*:", "", phenotype_id))]
    want <- pid[g == GENE, phenotype_id]
    message(sprintf("  %s: %d intron phenotypes for %s", TISSUE, length(want), GENE))
    if (!length(want)) next

    sq <- as.data.table(sd_ |> dplyr::filter(phenotype_id %in% want) |>
          dplyr::select(phenotype_id, variant_id, pval_nominal, slope, slope_se) |>
          dplyr::collect())
    sq <- merge(sq, br, by = "variant_id")
    parts <- tstrsplit(sq$variant_id, "_", fixed = TRUE)
    sq[, `:=`(g_ref = parts[[3]], g_alt = parts[[4]])]
    sq <- merge(sq, bim2, by = "rsid")
    sq <- sq[((g_ref == REF & g_alt == ALT) | (g_ref == ALT & g_alt == REF))
             & !is_ambig(g_ref, g_alt)]
    sq[, z_ref := (slope / slope_se) * fifelse(g_alt == REF, 1, -1)]
    sq <- sq[is.finite(z_ref)]

    for (vn in names(variants)) {
        V <- variants[[vn]]
        if (!length(V$fit$sets$cs)) next
        for (PH in unique(sq$phenotype_id)) {
            qp <- sq[phenotype_id == PH][!duplicated(rsid)]
            shared <- intersect(qp$rsid, V$d$rsid)
            if (length(shared) < MIN_SHARED) next
            qp <- qp[match(shared, rsid)]
            Rq <- V$R[shared, shared, drop = FALSE]
            fq <- tryCatch(susie_rss(z = qp$z_ref, R = Rq, n = N_QTL, L = L_MAX,
                                     estimate_residual_variance = FALSE),
                           error = function(e) NULL)
            if (is.null(fq) || !length(fq$sets$cs)) next
            fq <- coloc::annotate_susie(fq, shared, Rq)
            for (pp in P12_SWEEP) {
                cr <- tryCatch(
                    suppressWarnings(coloc::coloc.susie(V$fit, fq, p12 = pp)),
                    error = function(e) NULL)
                if (is.null(cr) || is.null(cr$summary) || !nrow(cr$summary)) next
                s <- as.data.table(cr$summary)
                if (!"PP.H4.abf" %in% names(s)) next
                s[, `:=`(variant = vn, phenotype_id = PH, tissue = TISSUE,
                         gene = GENE, n_shared = length(shared), p12 = pp)]
                coloc_rows[[length(coloc_rows)+1]] <- s
            }
        }
    }
}
coloc_dt <- if (length(coloc_rows)) rbindlist(coloc_rows, fill = TRUE) else
            data.table(variant = character())

write_parquet(convergence, file.path(OUTD, "convergence.parquet"))
write_parquet(purity,      file.path(OUTD, "credible_sets.parquet"))
if (nrow(coloc_dt)) write_parquet(coloc_dt, file.path(OUTD, "coloc_stability.parquet"))

message("\n==== convergence ====")
print(convergence[, .(variant, converged, niter, n_snp, n_cs, cross_cs_max_abs_r)])
message("\n==== credible sets ====")
if (nrow(purity) && "cs_size" %in% names(purity))
    print(purity[, .(variant, cs, cs_size, purity_min, lead_rsid, lead_pip,
                     jaccard_vs_primary)])
message("\n==== coloc.susie PP4 by tissue and variant (max over introns) ====")
if (nrow(coloc_dt)) {
    pr <- coloc_dt[p12 == P12]
    if (nrow(pr)) {
        best <- pr[, .(best_PP4 = max(PP.H4.abf),
                       best_pheno = phenotype_id[which.max(PP.H4.abf)],
                       n_pairs = .N), by = .(tissue, variant)]
        print(best)
    }
    message("\n---- PP4 of the winning intron across the p12 sweep ----")
    print(dcast(coloc_dt[, .(PP4 = max(PP.H4.abf)), by = .(tissue, variant, p12)],
                tissue + variant ~ p12, value.var = "PP4"))
} else message("  (no coloc rows produced)")

message(sprintf("\nwrote %s", OUTD))
cat("\nReproducibility\n"); print(Sys.time()); print(sessionInfo())
