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
## The fit itself is the one `03b.coloc_clpp.R` already performs and validates: identical
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
## Usage:  Rscript 05_genetic_anchoring/_h/03c.coloc_gwas_susie.R <analysis>
## Output: 05_genetic_anchoring/_m/coloc_signal_susie/gwas_susie/<analysis>/
##         (or .../coloc_signal_susie/sensitivity/max_snps_<N>/gwas_susie/<analysis>/ when
##          COLOC_GWAS_MAX_SNPS is not the primary 12,000)
##           <LOCUS_ID>.rds        -- annotated susie fit + snp set + n
##           _status.tsv           -- one row per locus: fitted / skipped and why
suppressPackageStartupMessages({
    library(data.table)
    library(susieR)
    library(coloc)
})

args     <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: 03c.coloc_gwas_susie.R <analysis>")
ANALYSIS <- args[1]
ROOT     <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
CDIR     <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc", ANALYSIS)
SUSIE_D  <- file.path(CDIR, "susie")
PANEL_DIR <- "/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink"

## Identical to 03b.coloc_clpp.R, deliberately: the two estimators must not differ because
## one of them quietly fine-mapped a different SNP set.
L_MAX    <- 10L
MIN_SNPS <- 20L

## MAX_SNPS is a COMPUTE guard, not a statistical one: susie_rss is O(L*p^2) per
## iteration and the locus LD is held dense, so a big locus is expensive rather than
## invalid. The 12,000 default is inherited from 03b.coloc_clpp.R and must stay the
## default so the two estimators keep fine-mapping the same SNP sets.
##
## It is nonetheless load-bearing on the RESULT, not just the runtime: 61 of the 579
## loci exceed it (median 14,664 SNPs, max 28,677) and are dropped before any fit, so
## they are absent from the signal-level layer entirely -- silently, since a dropped
## locus looks exactly like a locus with no credible set unless _status.tsv is read.
## PICALM's AD locus (locus60_chr11, 15,713 SNPs) is one of them, which is why PICALM
## has no coloc.susie row in either sQTL arm.
##
## Raising it is therefore a real recovery run, and env-overridable so that run is
## explicit and reproducible instead of an edited constant:
##   COLOC_GWAS_MAX_SNPS=30000 sbatch ... 03c.coloc_gwas_susie.sh
## Two cautions when raising it. Memory goes as p^2 (30,000 SNPs -> 7.2 GB for the LD
## subset alone, before susie's own copies; the AD recovery peaked at 48.2 GB and needs
## --cpus-per-task=32).
##
## And the guard turns out not to be statistically neutral. MEASURED, aging__ad,
## 2026-09-10: loci exceeding the guard were EMPIRICALLY ENRICHED for greater
## GWAS-reference-LD inconsistency. All 11 recovered loci had s_rss >= 0.310 against a
## median of 0.255 over all 49 fitted, and the largest sat at the top of the whole
## distribution (locus88_chr17, 28,677 SNPs, s_rss 0.739; locus96_chr19, s_rss 0.699).
## That is a correlation between a compute threshold and statistical difficulty, not an
## explanation of it -- this run does not establish WHY the two travel together, and
## writing that they do because large loci sit in long-range LD would assert a mechanism
## the data does not support. Read s_rss per locus before believing any recovered
## credible set, and do not promote a high-s_rss recovered locus to a headline.
MAX_SNPS <- {
    v <- Sys.getenv("COLOC_GWAS_MAX_SNPS", "12000")
    n <- suppressWarnings(as.integer(v))
    if (is.na(n) || n < MIN_SNPS) stop("COLOC_GWAS_MAX_SNPS not a usable integer: ", v)
    n
}
message("MAX_SNPS = ", MAX_SNPS)

## A non-default guard is a scoped SENSITIVITY arm and gets its own root. Before this split
## existed, the 2026-09-10 aging__ad recovery at 30,000 overwrote the 12,000 primary cache
## (and, through stage B, the primary AD shards) in place. Stage B and
## `coloc_signal_susie --stage meta --max-snps` resolve the same root from the same value.
MAX_SNPS_PRIMARY <- 12000L
SIG_ROOT <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc_signal_susie")
if (MAX_SNPS != MAX_SNPS_PRIMARY)
    SIG_ROOT <- file.path(SIG_ROOT, "sensitivity", sprintf("max_snps_%d", MAX_SNPS))
OUTD <- file.path(SIG_ROOT, "gwas_susie", ANALYSIS)
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
## Every locus is refit on each run, so clear the previous run's fits first. Locus ids are
## positional: after a coloc_prep re-prep a locus this run skips (<MIN_SNPS, >MAX_SNPS,
## susie error) would otherwise keep an older fit for a DIFFERENT region, which stage B
## (06b) reads by id (2026-09-18: 50 of 89 aging__ad fits were left over).
unlink(list.files(OUTD, pattern = "\\.rds$", full.names = TRUE))
message("GWAS SuSiE cache -> ", OUTD)

loci <- fread(file.path(SUSIE_D, "loci_testable.tsv"))
excl_f <- file.path(CDIR, "exclude_regions.tsv")
excl   <- if (file.exists(excl_f)) fread(excl_f) else
          data.table(chr = 6L, start = 25e6, stop = 34e6)

is_ambig <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

## An LD matrix is opened lazily and only the block a fit uses is read, row by row: a
## whole-file readBin needs n*n in R's 32-bit integer range (<= 46,340 SNPs), which the
## recurrence-1 aging loci exceed (54k-82k), and loading them to fit <= MAX_SNPS is waste.
ld_open <- function(lid) {
    vf <- file.path(SUSIE_D, paste0(lid, ".unphased.vcor1.bin.vars"))
    bf <- file.path(SUSIE_D, paste0(lid, ".unphased.vcor1.bin"))
    if (!file.exists(vf) || !file.exists(bf)) return(NULL)
    list(bf = bf, vars = readLines(vf))
}
ld_block <- function(ld, keep) {
    n <- length(ld$vars); idx <- match(keep, ld$vars)
    stopifnot(!anyNA(idx))
    con <- file(ld$bf, "rb"); on.exit(close(con))
    if (n <= 46340L) {  # fits one readBin: the original whole-matrix read
        R <- matrix(readBin(con, "numeric", n = n * n, size = 4), n, n, byrow = TRUE)
        R <- R[idx, idx, drop = FALSE]; dimnames(R) <- list(keep, keep); return(R)
    }
    R <- matrix(NA_real_, length(idx), length(idx), dimnames = list(keep, keep))
    for (k in seq_along(idx)) {  # row-major float32; offsets as double past 2 GB
        seek(con, (as.numeric(idx[k]) - 1) * n * 4)
        R[k, ] <- readBin(con, "numeric", n = n, size = 4)[idx]
    }
    R
}

status <- list()
for (i in seq_len(nrow(loci))) {
    L <- loci[i]; lid <- L$LOCUS_ID
    rec <- list(LOCUS_ID = lid, chr = L$chr, n_snp = NA_integer_,
                fitted = FALSE, n_cs = NA_integer_, reason = NA_character_,
                s_rss = NA_real_)

    R <- ld_open(lid)
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
    bim <- bim[rsid %in% R$vars, .(rsid, REF = A2, ALT = A1)]
    dt  <- merge(gwas, bim, by = "rsid")
    dt  <- dt[((a1==REF & a2==ALT) | (a1==ALT & a2==REF)) & !is_ambig(a1, a2)]
    dt[, z_ref := z * fifelse(a1==REF, 1, -1)]
    dt  <- dt[rsid %in% R$vars][!duplicated(rsid)]
    rec$n_snp <- nrow(dt)

    if (nrow(dt) < MIN_SNPS) {
        rec$reason <- sprintf("<%d usable SNPs", MIN_SNPS)
        status[[length(status)+1]] <- rec; rm(R, gwas, bim, dt); gc(verbose=FALSE); next
    }
    if (nrow(dt) > MAX_SNPS) {
        rec$reason <- sprintf("%d SNPs > MAX_SNPS (long-range LD)", nrow(dt))
        status[[length(status)+1]] <- rec; rm(R, gwas, bim, dt); gc(verbose=FALSE); next
    }

    Rk <- ld_block(R, dt$rsid)
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
