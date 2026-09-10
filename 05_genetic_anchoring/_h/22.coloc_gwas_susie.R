#!/usr/bin/env Rscript
## Signal-level coloc, stage A -- fit and CACHE the per-locus GWAS SuSiE.
##
## `coloc.susie` needs a fine-mapped object on BOTH sides. The GWAS side depends only on
## (locus, trait), not on the QTL tissue or modality, so fitting it once here and caching
## it saves the 13 tissue x 2 modality tasks in stage B from re-deriving the identical
## fit 26 times per locus -- and, more importantly, guarantees every one of those tasks
## colocalizes against the *same* GWAS fit rather than 26 independent re-fits that could
## differ by convergence.
##
## The fit itself is the one `10.coloc_clpp.R` already performs and validates: identical
## LD matrices, identical allele re-signing to the panel REF, identical strand-ambiguity
## and long-range-LD exclusions, identical L and MAX_SNPS guards. This script is that
## code with the CLPP step removed and the fit object kept.
##
## An `estimate_s_rss` diagnostic is recorded per locus. It quantifies the
## consistency between the z-scores and the REFERENCE LD panel, which is the single
## largest threat to this whole layer: 1000G EUR is a proxy for both the GWAS and the
## GTEx sample, and a locus where that proxy is bad will fine-map to credible sets that
## are artefacts. Recording the diagnostic is what makes that checkable afterwards
## instead of assumed away.
##
## Usage:  Rscript 05_genetic_anchoring/_h/22.coloc_gwas_susie.R <analysis>
## Output: 05_genetic_anchoring/_m/coloc_signal_susie/gwas_susie/<analysis>/
##           <LOCUS_ID>.rds        -- annotated susie fit + snp set + n
##           _status.tsv           -- one row per locus: fitted / skipped and why
suppressPackageStartupMessages({
    library(data.table)
    library(susieR)
    library(coloc)
})

args     <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: 22.coloc_gwas_susie.R <analysis>")
ANALYSIS <- args[1]
ROOT     <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
CDIR     <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc", ANALYSIS)
SUSIE_D  <- file.path(CDIR, "susie")
OUTD     <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc_signal_susie",
                      "gwas_susie", ANALYSIS)
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
PANEL_DIR <- "/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink"

## Identical to 10.coloc_clpp.R, deliberately: the two estimators must not differ because
## one of them quietly fine-mapped a different SNP set.
L_MAX    <- 10L
MAX_SNPS <- 12000L
MIN_SNPS <- 20L

loci <- fread(file.path(SUSIE_D, "loci_testable.tsv"))
excl_f <- file.path(CDIR, "exclude_regions.tsv")
excl   <- if (file.exists(excl_f)) fread(excl_f) else
          data.table(chr = 6L, start = 25e6, stop = 34e6)

is_ambig <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

read_ld <- function(lid) {
    vf <- file.path(SUSIE_D, paste0(lid, ".unphased.vcor1.bin.vars"))
    bf <- file.path(SUSIE_D, paste0(lid, ".unphased.vcor1.bin"))
    if (!file.exists(vf) || !file.exists(bf)) return(NULL)
    vars <- readLines(vf); n <- length(vars)
    con <- file(bf, "rb"); m <- readBin(con, "numeric", n = n*n, size = 4); close(con)
    R <- matrix(m, n, n, byrow = TRUE); dimnames(R) <- list(vars, vars); R
}

status <- list()
for (i in seq_len(nrow(loci))) {
    L <- loci[i]; lid <- L$LOCUS_ID
    rec <- list(LOCUS_ID = lid, chr = L$chr, n_snp = NA_integer_,
                fitted = FALSE, n_cs = NA_integer_, reason = NA_character_,
                s_rss = NA_real_)

    R <- read_ld(lid)
    if (is.null(R)) {
        rec$reason <- "no LD matrix"; status[[length(status)+1]] <- rec; next
    }
    gwas <- fread(file.path(SUSIE_D, paste0(lid, ".gwas.tsv")))
    for (j in seq_len(nrow(excl)))
        if (L$chr == excl$chr[j])
            gwas <- gwas[!(pos >= excl$start[j] & pos <= excl$stop[j])]

    bim <- fread(file.path(PANEL_DIR, sprintf("1000G.EUR.QC.%d.bim", L$chr)),
                 header = FALSE, col.names = c("bchr","rsid","cm","bpos","A1","A2"))
    ## plink1 .bed convention: REF = A2 (col6), ALT = A1 (col5)
    bim <- bim[rsid %in% rownames(R), .(rsid, REF = A2, ALT = A1)]
    dt  <- merge(gwas, bim, by = "rsid")
    dt  <- dt[((a1==REF & a2==ALT) | (a1==ALT & a2==REF)) & !is_ambig(a1, a2)]
    dt[, z_ref := z * fifelse(a1==REF, 1, -1)]
    dt  <- dt[rsid %in% rownames(R)][!duplicated(rsid)]
    rec$n_snp <- nrow(dt)

    if (nrow(dt) < MIN_SNPS) {
        rec$reason <- sprintf("<%d usable SNPs", MIN_SNPS)
        status[[length(status)+1]] <- rec; rm(R, gwas, bim, dt); gc(verbose=FALSE); next
    }
    if (nrow(dt) > MAX_SNPS) {
        rec$reason <- sprintf("%d SNPs > MAX_SNPS (long-range LD)", nrow(dt))
        status[[length(status)+1]] <- rec; rm(R, gwas, bim, dt); gc(verbose=FALSE); next
    }

    Rk <- R[dt$rsid, dt$rsid, drop = FALSE]
    rm(R); gc(verbose = FALSE)
    nn <- max(1000, round(L$min_neff))

    fit <- tryCatch(
        susie_rss(z = dt$z_ref, R = Rk, n = nn, L = L_MAX,
                  estimate_residual_variance = FALSE),
        error = function(e) { message("  susie failed ", lid, ": ",
                                      conditionMessage(e)); NULL })
    if (is.null(fit)) {
        rec$reason <- "susie_rss error"
        status[[length(status)+1]] <- rec; rm(Rk, dt); gc(verbose=FALSE); next
    }
    if (!length(fit$sets$cs)) {
        ## No GWAS credible set means `coloc.susie` has nothing to colocalize against, so
        ## the locus is a fallback-to-abf cell downstream. Recorded and skipped WITHOUT
        ## paying for the diagnostic or the RDS -- this is common for the smaller GWAS
        ## (LBD n_eff ~6.6k), and eating an O(p^3) eigen on every such locus would
        ## dominate the run for a number nothing reads.
        rec$n_cs <- 0L; rec$reason <- "no GWAS credible set"
        status[[length(status)+1]] <- rec; rm(Rk, fit, dt); gc(verbose=FALSE); next
    }

    ## `estimate_s_rss` only, and only for loci that actually fine-mapped. It is one
    ## O(p^3) eigen decomposition yielding the standard z-vs-reference-LD consistency
    ## statistic. The full `kriging_rss` pass is deliberately NOT run: it repeats that
    ## cost to produce a per-SNP outlier table this layer does not act on.
    s_rss <- tryCatch(susieR::estimate_s_rss(z = dt$z_ref, R = Rk, n = nn),
                      error = function(e) NA_real_)
    kr <- list(s = s_rss)
    rec$s_rss <- s_rss

    ## coloc.bf_bf matches SNPs by the column names of lbf_variable, so the fit must be
    ## annotated before it is usable -- an unnamed fit silently colocalizes nothing.
    fit <- coloc::annotate_susie(fit, dt$rsid, Rk)
    rec$fitted <- TRUE
    rec$n_cs   <- length(fit$sets$cs)
    rec$reason <- "ok"

    ## Keep only what coloc.susie reads plus provenance. The LD matrix is never written:
    ## it is 100s of MB per locus and is reconstructible from the panel.
    saveRDS(list(fit = fit, snps = dt$rsid, n = nn, analysis = ANALYSIS,
                 LOCUS_ID = lid, chr = L$chr, kriging = kr),
            file.path(OUTD, paste0(lid, ".rds")), compress = "xz")
    message(sprintf("  %-18s snps=%5d  cs=%d  max_pip=%.3f  s_rss=%.3f",
                    lid, nrow(dt), rec$n_cs, max(fit$pip), kr$s))
    status[[length(status)+1]] <- rec
    rm(Rk, fit, dt); gc(verbose = FALSE)
}

st <- rbindlist(status, fill = TRUE)
fwrite(st, file.path(OUTD, "_status.tsv"), sep = "\t")
message(sprintf("\n[%s] %d/%d loci fitted; %d carry >=1 credible set",
                ANALYSIS, sum(st$fitted), nrow(st),
                sum(st$fitted & st$n_cs > 0, na.rm = TRUE)))
cat("\nReproducibility\n"); print(Sys.time()); print(sessionInfo())
