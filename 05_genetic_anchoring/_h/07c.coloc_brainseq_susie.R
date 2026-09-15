#!/usr/bin/env Rscript
## BrainSEQ signal-level coloc: QTL SuSiE on in-sample EA LD + coloc.susie, and coloc.abf for
## every cell, for ONE (analysis, region) task.
##
## For every (locus, gene) target this fits the gene's cis QTL on BOTH BrainSEQ axes -- the
## switch coordinate S_g and the abundance channel A_g, mapped in the same EA donors -- and
## colocalizes SIGNAL AGAINST SIGNAL with the cached GWAS SuSiE from stage A
## (03c.coloc_gwas_susie.R), the identical fit the GTEx layer (23.*) uses.
##
## THE QTL-SIDE LD IS IN-SAMPLE, AND RESIDUALIZED LIKE THE Z-SCORES
##   The GTEx layer has to fine-map its QTL on the 1000G EUR panel and guard the result with an
##   agreement filter against GTEx's own credible sets. BrainSEQ's genotypes are here, so the LD
##   comes from exactly the donors the region was mapped in (regions.tsv / donors_<region>.txt,
##   written by `coloc_brainseq --stage prep`).
##
##   It is the correlation of ALT dosages RESIDUALIZED on that axis's covariate matrix
##   (covariates_used_<axis>.txt), not of raw genotypes. tensorQTL regresses the covariates out
##   of genotype and phenotype alike, so the residualized genotypes are the matrix the z-scores
##   actually came from. Raw-genotype LD is a mismatch that grows as donors shrink relative to
##   covariates (41-42 here): measured 2026-09-11 on DLPFC (169 donors), it doubled
##   estimate_s_rss (0.025 vs 0.013) and made susie_rss abort ("estimated prior variance is
##   unreasonably large") in 54 cells, all in DLPFC or hippocampus. On a cell that fit either
##   way, the two matrices gave the same credible set (s_rss 0.0108 vs 0.0103). The axes have
##   different hidden factors, so each axis gets its own matrix.
##
##   The GWAS side stays on the stage-A panel fit; coloc.susie matches the two fits by rsID, so
##   each fit only has to be consistent with its own LD.
##
## ALLELE ORIENTATION
##   tensorQTL's slope is per the ALT allele (the positive control confirmed it counts ALT), and
##   the dosages are oriented to ALT from the plink2 export's counted-allele suffix, so z and LD
##   share one coding. Alleles are still checked against the 1000G panel, and strand-ambiguous or
##   mismatched variants dropped, so an rsID only joins the GWAS fit when it is the same variant.
##
## Usage:  Rscript 05_genetic_anchoring/_h/07c.coloc_brainseq_susie.R <analysis> <region> [chr]
##         (a chr argument is a smoke test: it writes under <root>/smoke/, never into the grid)
## Env:    COLOC_BRAINSEQ_ARM  genotype arm (only "ea_only")
##         PLINK2              plink2 binary (default: the eqtl env's)
## Output: 05_genetic_anchoring/_m/coloc_brainseq/<arm>/{susie,abf,fits}/<analysis>__<region>.parquet
suppressPackageStartupMessages({
    library(data.table); library(arrow); library(susieR); library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: 07c.coloc_brainseq_susie.R <analysis> <region> [chr]")
ANALYSIS <- args[1]
REGION   <- args[2]
ONE_CHR  <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
ROOT     <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
ARM      <- Sys.getenv("COLOC_BRAINSEQ_ARM", unset = "ea_only")
if (ARM != "ea_only")
    stop("BrainSEQ coloc runs on the ea_only arm only: the GWAS and the GWAS-side LD are European")
PLINK2   <- Sys.getenv("PLINK2", unset = "/ocean/projects/bio260021p/shared/opt/envs/eqtl/bin/plink2")
if (!file.exists(PLINK2)) stop("plink2 not found at ", PLINK2)
THREADS  <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
MEM_MB   <- THREADS * 1000L   # half of the PSC allocation; R holds the LD matrices in the rest
PANEL_DIR <- "/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink"
GENO     <- file.path("/ocean/projects/bio260021p/shared/resources/processed-data/genotypes/qtl", ARM)

M     <- file.path(ROOT, "05_genetic_anchoring", "_m")
CROOT <- file.path(M, "coloc_brainseq", ARM)
QDIR  <- file.path(M, "brainseq_switch_qtl", ARM, REGION, "qtl")
GDIR  <- file.path(M, "coloc_signal_susie", "gwas_susie", ANALYSIS)
CDIR  <- file.path(M, "coloc", ANALYSIS)
for (d in c(CROOT, QDIR, GDIR, CDIR)) if (!dir.exists(d)) stop("missing ", d)
OUT_BASE <- if (is.na(ONE_CHR)) CROOT else file.path(CROOT, "smoke")

TMPD <- Sys.getenv("LOCAL", unset = "")
if (!nzchar(TMPD) || !dir.exists(TMPD)) TMPD <- tempdir()
TMPD <- file.path(TMPD, sprintf("coloc_brainseq_%s_%s_%d", ANALYSIS, REGION, Sys.getpid()))
dir.create(TMPD, recursive = TRUE, showWarnings = FALSE)

L_MAX       <- 10L
MIN_SHARED  <- 50L        # emit the row; the meta stage applies the real (100) floor
P12_SWEEP   <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
P12_PRIMARY <- 1e-5
SDY         <- 1          # phenotypes are inverse-normal transformed before mapping
AXES        <- c(S_g = "switch", A_g = "abundance")

regions <- fread(file.path(CROOT, "regions.tsv"))
rg <- regions[region == REGION]
if (nrow(rg) != 1L) stop("region ", REGION, " is not in regions.tsv")
N_QTL <- as.integer(rg$n_donors)
KEEP  <- file.path(CROOT, rg$donors_file)
message(sprintf("[%s / %s] arm=%s donors=%d threads=%d", ANALYSIS, REGION, ARM, N_QTL, THREADS))

## The covariate matrix each axis was mapped with: donors x covariates, as tensorQTL was given it.
COV <- lapply(AXES, function(ft) {
    f <- file.path(QDIR, sprintf("covariates_used_%s.txt", ft))
    x <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
    if (nrow(x) != N_QTL) stop(f, " has ", nrow(x), " donors, regions.tsv says ", N_QTL)
    x
})

targets <- as.data.table(read_parquet(file.path(CROOT, "targets.parquet")))
targets <- targets[analysis == ANALYSIS]
if (!nrow(targets)) { message("no targets for this analysis; nothing to do"); quit(status = 0) }
phen <- as.data.table(read_parquet(file.path(CROOT, "phenotypes.parquet")))
phen <- phen[region == REGION]

## Stage A writes an RDS only for a locus whose GWAS fine-mapped, so a missing RDS alone does
## not say why. Its status table does, and the fit outcomes below carry that reason.
gst <- fread(file.path(GDIR, "_status.tsv"))
GWAS_REASON <- setNames(
    ifelse(grepl("MAX_SNPS", gst$reason), "gwas_locus_over_max_snps",
    ifelse(grepl("credible set", gst$reason), "gwas_no_credible_set",
    ifelse(grepl("usable SNPs", gst$reason), "gwas_too_few_snps", "gwas_not_fit"))),
    gst$LOCUS_ID)

gm <- fread(file.path(CROOT, "gwas_meta.tsv"))
TRAIT_NS <- setNames(
    lapply(seq_len(nrow(gm)), function(i) list(s = gm$s[i], n_gwas = gm$n_gwas[i])), gm$trait)

is_ambig <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

## Per-locus GWAS with the coloc layer's long-range-LD exclusions (as 05a.coloc_modality_abf.R).
get_gwas <- function(lid, chr) {
    f <- file.path(CDIR, "susie", paste0(lid, ".gwas.tsv"))
    if (!file.exists(f)) return(data.table())
    g <- fread(f)
    ef <- file.path(CDIR, "exclude_regions.tsv")
    if (file.exists(ef)) {
        ex <- fread(ef)
        for (j in seq_len(nrow(ex)))
            if (chr == ex$chr[j]) g <- g[!(pos >= ex$start[j] & pos <= ex$stop[j])]
    }
    g <- g[!is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
    g[, .(rsid, gwas_beta = beta, gwas_varbeta = se^2, gwas_p = p)]
}

## REF/ALT per rsID from the arm's .pvar. An rsID that occurs twice (split multi-allelics) is
## dropped entirely: it cannot key a join or an LD row unambiguously.
read_pvar <- function(CH) {
    p <- fread(file.path(GENO, sprintf("chr%d.pvar", CH)), skip = "#CHROM",
               select = c("ID", "REF", "ALT"), colClasses = "character")
    dup <- p$ID[duplicated(p$ID)]
    setnames(p[!ID %in% dup], c("variant_id", "bs_ref", "bs_alt"))
}

read_qtl <- function(CH, AX, ids) {
    files <- Sys.glob(file.path(QDIR, sprintf("%s.chr%d.cis_qtl_pairs.*.parquet", AXES[[AX]], CH)))
    if (!length(files) || !length(ids)) return(data.table())
    as.data.table(open_dataset(files) |>
        dplyr::filter(phenotype_id %in% ids) |>
        dplyr::select(phenotype_id, variant_id, af, pval_nominal, slope, slope_se) |>
        dplyr::collect())
}

prep_qtl <- function(q, pv, bim) {
    if (!nrow(q)) return(q)
    q <- merge(q, pv, by = "variant_id")
    q <- merge(q, bim, by.x = "variant_id", by.y = "rsid")
    q <- q[((bs_ref == REF & bs_alt == ALT) | (bs_ref == ALT & bs_alt == REF))
           & !is_ambig(bs_ref, bs_alt)]
    q[, `:=`(rsid = variant_id, z = slope / slope_se,
             qtl_beta = slope, qtl_varbeta = slope_se^2)]
    q[is.finite(z) & qtl_varbeta > 0]
}

## ALT dosages (donors x SNPs) for `snps` in the region's mapped donors, from the arm's pgen.
read_dosage <- function(CH, snps, tag, alt_of) {
    pre <- file.path(TMPD, tag)
    writeLines(snps, paste0(pre, ".snps"))
    rc <- system2(PLINK2, c("--pfile", file.path(GENO, sprintf("chr%d", CH)),
                            "--keep", KEEP, "--extract", paste0(pre, ".snps"),
                            "--export", "A", "--threads", THREADS, "--memory", MEM_MB,
                            "--out", pre), stdout = FALSE, stderr = FALSE)
    rf <- paste0(pre, ".raw")
    if (rc != 0 || !file.exists(rf)) {
        lg <- paste0(pre, ".log")
        message("  plink2 export failed for ", tag, " (rc=", rc, ")",
                if (file.exists(lg)) paste0(": ", paste(tail(readLines(lg), 3), collapse = " | ")) else "")
        unlink(Sys.glob(paste0(pre, "*")))
        return(NULL)
    }
    raw <- fread(rf)
    unlink(Sys.glob(paste0(pre, "*")))
    X <- as.matrix(raw[, -(1:6)])
    rownames(X) <- as.character(raw$IID)
    cid <- colnames(X)
    rs <- sub("_[^_]+$", "", cid)
    counted <- sub("^.*_", "", cid)
    ## Orient every column to ALT dosage, whichever allele plink2 counted.
    X <- sweep(X, 2, ifelse(counted == alt_of[rs], 1, -1), "*")
    colnames(X) <- rs
    miss <- which(is.na(X), arr.ind = TRUE)
    if (nrow(miss)) X[miss] <- colMeans(X, na.rm = TRUE)[miss[, 2]]
    X
}

## LD of the genotypes an axis's z-scores came from: dosages with that axis's covariates
## (plus an intercept) regressed out. A SNP left with no variance has no defined correlation
## and is dropped.
residual_ld <- function(X, cov) {
    X <- X[rownames(cov), , drop = FALSE]
    Q <- qr.Q(qr(cbind(1, cov)))
    Xr <- X - Q %*% crossprod(Q, X)
    ok <- apply(Xr, 2, sd) > 1e-8
    R <- cor(Xr[, ok, drop = FALSE])
    R[!is.finite(R)] <- 0
    diag(R) <- 1
    R
}

run_abf <- function(qg, g, tn) {
    m <- merge(qg[, .(rsid, qtl_beta, qtl_varbeta)], g, by = "rsid")[!duplicated(rsid)]
    if (nrow(m) < MIN_SHARED) return(NULL)
    d_gwas <- list(beta = m$gwas_beta, varbeta = m$gwas_varbeta, snp = m$rsid,
                   type = "cc", s = tn$s, N = tn$n_gwas)
    d_qtl  <- list(beta = m$qtl_beta, varbeta = m$qtl_varbeta, snp = m$rsid,
                   type = "quant", sdY = SDY, N = N_QTL)
    out <- vector("list", length(P12_SWEEP))
    for (i in seq_along(P12_SWEEP)) {
        r <- tryCatch({
                 invisible(capture.output(
                     z <- suppressWarnings(suppressMessages(
                         coloc.abf(d_gwas, d_qtl, p12 = P12_SWEEP[i])))))
                 z
             }, error = function(e) NULL)
        if (is.null(r)) next
        s <- as.list(r$summary)
        top <- as.data.table(r$results)[which.max(SNP.PP.H4)]
        out[[i]] <- data.table(
            p12 = P12_SWEEP[i], nsnps = as.integer(s$nsnps),
            PP0 = s$PP.H0.abf, PP1 = s$PP.H1.abf, PP2 = s$PP.H2.abf,
            PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf,
            top_rsid = top$snp, top_pp4 = top$SNP.PP.H4,
            gwas_min_p = min(m$gwas_p, na.rm = TRUE))
    }
    rbindlist(out)
}

## One (gwas signal x qtl signal) posterior table for one cell, swept over p12.
coloc_cell <- function(fit_g, fit_q, meta) {
    if (!length(fit_g$sets$cs) || !length(fit_q$sets$cs)) return(NULL)
    out <- vector("list", length(P12_SWEEP))
    for (k in seq_along(P12_SWEEP)) {
        r <- tryCatch(suppressWarnings(coloc::coloc.susie(fit_g, fit_q, p12 = P12_SWEEP[k])),
                      error = function(e) NULL)
        if (is.null(r) || is.null(r$summary) || !nrow(r$summary)) next
        s <- as.data.table(r$summary)
        if (!"PP.H4.abf" %in% names(s)) next
        out[[k]] <- data.table(
            p12 = P12_SWEEP[k], idx1 = s$idx1, idx2 = s$idx2,
            nsnps = s$nsnps, hit1 = s$hit1, hit2 = s$hit2,
            PP0 = s$PP.H0.abf, PP1 = s$PP.H1.abf, PP2 = s$PP.H2.abf,
            PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf)
    }
    out <- rbindlist(out, fill = TRUE)
    if (!nrow(out)) return(NULL)
    cbind(as.data.table(meta), out)
}

susie_rows <- list(); abf_rows <- list(); fit_rows <- list(); susie_errors <- character()
chrs <- if (is.na(ONE_CHR)) sort(unique(targets$chr)) else ONE_CHR
for (CH in chrs) {
    tg <- targets[chr == CH]
    if (!nrow(tg)) next
    genes <- unique(tg$gene)
    pv  <- read_pvar(CH)
    bim <- fread(file.path(PANEL_DIR, sprintf("1000G.EUR.QC.%d.bim", CH)), header = FALSE,
                 col.names = c("bchr", "rsid", "cm", "bpos", "A1", "A2"))
    bim <- bim[!rsid %in% rsid[duplicated(rsid)], .(rsid, REF = A2, ALT = A1)]
    Q <- list()
    for (AX in names(AXES)) {
        ph <- phen[modality == AX & gene %in% genes]
        q  <- read_qtl(CH, AX, ph$phenotype_id)
        if (nrow(q)) q <- prep_qtl(merge(q, ph[, .(phenotype_id, gene)], by = "phenotype_id"), pv, bim)
        Q[[AX]] <- q
    }
    alt_of <- setNames(pv$bs_alt, pv$variant_id)
    rm(bim); gc(verbose = FALSE)

    for (lid in unique(tg$LOCUS_ID)) {
        lg <- tg[LOCUS_ID == lid]
        tn <- TRAIT_NS[[lg$trait[1]]]
        g  <- get_gwas(lid, CH)
        rds <- file.path(GDIR, paste0(lid, ".rds"))
        G <- if (file.exists(rds)) readRDS(rds) else NULL
        gwas_cs <- !is.null(G) && length(G$fit$sets$cs) > 0
        X <- NULL; RL <- list()
        if (gwas_cs) {
            U <- unique(unlist(lapply(Q, function(q)
                if (nrow(q)) q[gene %in% lg$gene & rsid %in% G$snps, rsid] else character())))
            if (length(U) >= MIN_SHARED) X <- read_dosage(CH, U, lid, alt_of)
        }

        for (gi in seq_len(nrow(lg))) {
            tr <- lg[gi]
            for (AX in names(AXES)) {
                q0 <- Q[[AX]]
                if (is.null(q0) || !nrow(q0)) next
                qg <- q0[gene == tr$gene][!duplicated(rsid)]
                if (!nrow(qg)) next
                base <- list(analysis = ANALYSIS, trait = tr$trait, LOCUS_ID = lid, chr = CH,
                             gene = tr$gene, symbol = tr$symbol, module_id = tr$module_id,
                             go_invisible = tr$go_invisible, tissue = REGION, region = REGION,
                             modality = AX, phenotype_id = qg$phenotype_id[1],
                             n_qtl_donors = N_QTL, qtl_min_p = min(qg$pval_nominal, na.rm = TRUE))

                if (nrow(g) && !is.null(tn)) {
                    a <- run_abf(qg, g, tn)
                    if (!is.null(a) && nrow(a))
                        abf_rows[[length(abf_rows) + 1]] <- cbind(as.data.table(base), a)
                }

                if (gwas_cs && !is.null(X) && is.null(RL[[AX]])) RL[[AX]] <- residual_ld(X, COV[[AX]])
                R <- RL[[AX]]
                outcome <- if (is.null(G)) {
                               if (lid %in% names(GWAS_REASON)) GWAS_REASON[[lid]] else "gwas_not_in_stage_a"
                           } else if (!gwas_cs) "gwas_no_credible_set"
                           else if (is.null(R)) "ld_unavailable" else NA_character_
                n_shared <- NA_integer_; n_cs_qtl <- NA_integer_
                if (is.na(outcome)) {
                    shared <- intersect(intersect(qg$rsid, G$snps), rownames(R))
                    n_shared <- length(shared)
                    if (n_shared < MIN_SHARED) {
                        outcome <- "too_few_shared_snps"
                    } else {
                        qp <- qg[match(shared, rsid)]
                        Rq <- R[shared, shared, drop = FALSE]
                        fq <- tryCatch(
                            susie_rss(z = qp$z, R = Rq, n = N_QTL, L = L_MAX,
                                      estimate_residual_variance = FALSE),
                            error = function(e) conditionMessage(e))
                        if (is.character(fq)) {
                            outcome <- "susie_rss_error"
                            susie_errors <- c(susie_errors, gsub("\\s+", " ", fq))
                        } else {
                            fq <- coloc::annotate_susie(fq, shared, Rq)
                            n_cs_qtl <- length(fq$sets$cs)
                            meta <- c(base, list(n_shared = n_shared, n_qtl_variants = nrow(qg),
                                                 n_cs_gwas = length(G$fit$sets$cs),
                                                 n_cs_qtl = n_cs_qtl, s_rss_gwas = G$kriging$s,
                                                 ld_source = "in_sample_residualized"))
                            r <- coloc_cell(G$fit, fq, meta)
                            outcome <- if (!n_cs_qtl) "no_qtl_credible_set"
                                       else if (is.null(r)) "coloc_susie_error" else "coloc_susie"
                            if (!is.null(r)) susie_rows[[length(susie_rows) + 1]] <- r
                        }
                        rm(fq, Rq)
                    }
                }
                fit_rows[[length(fit_rows) + 1]] <- data.table(
                    analysis = ANALYSIS, trait = tr$trait, LOCUS_ID = lid, gene = tr$gene,
                    tissue = REGION, modality = AX, outcome = outcome,
                    n_shared = n_shared, n_cs_qtl = n_cs_qtl)
            }
        }
        rm(X, RL, G); gc(verbose = FALSE)
        message(sprintf("  chr%-2d %-18s genes=%d  susie blocks=%d  abf blocks=%d",
                        CH, lid, nrow(lg), length(susie_rows), length(abf_rows)))
    }
    rm(Q, pv, alt_of); gc(verbose = FALSE)
}
unlink(TMPD, recursive = TRUE)

stem <- sprintf("%s__%s%s.parquet", ANALYSIS, REGION,
                if (is.na(ONE_CHR)) "" else sprintf(".chr%d", ONE_CHR))
write_kind <- function(rows, kind) {
    d <- file.path(OUT_BASE, kind); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(d, stem)
    if (length(rows)) {
        x <- rbindlist(rows, fill = TRUE)
        write_parquet(x, f)
        message(sprintf("  %-5s %8d rows -> %s", kind, nrow(x), f))
        return(x)
    }
    if (file.exists(f)) file.remove(f)    # never leave a stale shard from an earlier run
    message(sprintf("  %-5s no rows", kind))
    NULL
}
s <- write_kind(susie_rows, "susie")
a <- write_kind(abf_rows, "abf")
f <- write_kind(fit_rows, "fits")
if (!is.null(f)) print(f[, .N, by = .(modality, outcome)])
if (length(susie_errors)) {
    message("  susie_rss errors:")
    print(sort(table(substr(susie_errors, 1, 90)), decreasing = TRUE))
}
for (AX in names(AXES)) {
    ns <- if (is.null(s)) 0L else s[p12 == P12_PRIMARY & modality == AX & PP4 >= 0.8, uniqueN(gene)]
    na <- if (is.null(a)) 0L else a[p12 == P12_PRIMARY & modality == AX & nsnps >= 100 & PP4 >= 0.8, uniqueN(gene)]
    message(sprintf("  %s: genes with PP4 >= 0.8 at p12=%g -- coloc.susie %d, coloc.abf %d",
                    AX, P12_PRIMARY, ns, na))
}
cat("\nReproducibility\n"); print(Sys.time()); print(sessionInfo())
