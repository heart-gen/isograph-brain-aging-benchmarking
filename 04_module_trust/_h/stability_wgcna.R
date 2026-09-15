#!/usr/bin/env Rscript
# WGCNA within-cohort split-half stability (abundance-network reference ceiling).
#
# For one cohort+region: load the bundle-filtered gene expression (samples x genes),
# randomly split samples 50/50 over several seeds, run the SAME WGCNA module
# detection on each half (signed network, grey dropped — identical to the production
# baselines), and write each half's gene->module partition. The Python aggregator
# (isograph_benchmark.real_data.stability aggregate) then computes ARI/NMI uniformly
# for both methods.
#
# Usage:  Rscript stability_wgcna.R <cohort> <region> [seeds]
#   cohort in {brainseq, gtex};  region = bundle dir name;  seeds default 5.
suppressPackageStartupMessages({
    library(arrow)
    library(WGCNA)
    library(dplyr)
})

options(stringsAsFactors = FALSE)
WGCNA_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
if (is.na(WGCNA_THREADS) || WGCNA_THREADS < 2) WGCNA_THREADS <- 2L
enableWGCNAThreads(nThreads = WGCNA_THREADS)

project_root <- here::here()
if (!file.exists(file.path(project_root, ".here"))) {
    stop("Cannot locate project root (.here missing). Run from isograph-brain-aging-benchmarking/.")
}

MIN_MODULE_SIZE <- 30
R2_THRESHOLD <- 0.85
SEED_BASE <- 1000  # split seed = SEED_BASE + k, matching stability.py's SEED_BASE

# ── soft power for the signed topology blockwiseModules builds ───────────────────
select_soft_power <- function(datExpr) {
    sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed",
                             RsquaredCut = R2_THRESHOLD, verbose = 0)
    power <- sft$powerEstimate
    if (is.na(power)) power <- 12L
    as.integer(max(as.integer(power), 12L))
}

# ── one WGCNA partition (samples x genes -> gene/module df, grey excluded) ───────
detect_modules <- function(datExpr) {
    gs <- goodSamplesGenes(datExpr, verbose = 0)
    datExpr <- datExpr[gs$goodSamples, gs$goodGenes]
    power <- select_soft_power(datExpr)
    net <- blockwiseModules(
        datExpr, power = power, minModuleSize = MIN_MODULE_SIZE,
        TOMType = "signed", networkType = "signed", numericLabels = FALSE,
        mergeCutHeight = 0.25, verbose = 0
    )
    df <- data.frame(gene_id = colnames(datExpr), module_id = net$colors)
    df <- df[df$module_id != "grey", ]
    ms <- sort(table(df$module_id), decreasing = TRUE)
    lab <- setNames(sprintf("M%03d", seq_along(ms) - 1), names(ms))
    df$module_id <- lab[df$module_id]
    list(modules = df, power = power, n_modules = length(ms))
}

# ── cohort-specific expression loader -> datExpr (samples x genes) ───────────────
load_datExpr <- function(cohort, region) {
    if (cohort == "gtex") {
        proc_dir <- file.path(project_root, "inputs", "processed", "gtex_v11", region)
        bundle_dir <- file.path(project_root, "inputs", "bundles", "gtex_v11_brain", region)
        expr_tbl <- read_parquet(file.path(proc_dir, "gene_tpm.parquet"))
        id_col <- "Name"; to_log <- function(m) log2(m + 1)  # log2(TPM + 1)
    } else if (cohort == "brainseq") {
        proc_dir <- file.path(project_root, "inputs", "processed", "brainseq", region)
        bundle_dir <- file.path(project_root, "inputs", "bundles", "brainseq_v1", region)
        expr_tbl <- read_parquet(file.path(proc_dir, "gene_counts.parquet"))
        id_col <- "Geneid"
        to_log <- function(m) {                              # log2(CPM + 1)
            lib <- colSums(m)
            log2(sweep(m, 2, lib / 1e6, FUN = "/") + 1)
        }
    } else {
        stop(sprintf("unknown cohort '%s'", cohort))
    }
    bundle_samples <- read_parquet(file.path(bundle_dir, "samples.parquet"))
    bundle_genes <- read_parquet(file.path(bundle_dir, "genes.parquet"))
    gene_avail <- intersect(bundle_genes$gene_id, expr_tbl[[id_col]])
    expr_tbl <- expr_tbl[match(gene_avail, expr_tbl[[id_col]]), ]
    sample_cols <- intersect(bundle_samples$sample_id, names(expr_tbl))
    mat <- as.matrix(expr_tbl[, sample_cols])
    rownames(mat) <- gene_avail
    t(to_log(mat))  # samples x genes
}

# ── main ─────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: stability_wgcna.R <cohort> <region> [seeds]")
cohort <- args[1]; region <- args[2]
seeds <- if (length(args) >= 3) as.integer(args[3]) else 5L

out_dir <- file.path(project_root, "04_module_trust", "_m", "stability", "partitions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

datExpr <- load_datExpr(cohort, region)
n <- nrow(datExpr)
cat(sprintf("[%s/%s] %d samples x %d genes | %d split-half seeds\n",
            cohort, region, n, ncol(datExpr), seeds))

for (k in 0:(seeds - 1)) {
    set.seed(SEED_BASE + k)
    perm <- sample.int(n)
    h <- floor(n / 2)
    halves <- list(A = perm[1:h], B = perm[(h + 1):n])
    for (half in names(halves)) {
        idx <- halves[[half]]
        res <- tryCatch(detect_modules(datExpr[idx, , drop = FALSE]),
                        error = function(e) { message(sprintf("  seed%d %s FAILED: %s", k, half, conditionMessage(e))); NULL })
        if (is.null(res)) next
        df <- res$modules
        df$method <- "wgcna"; df$cohort <- cohort; df$region <- region
        df$seed <- k; df$half <- half
        fname <- sprintf("wgcna__%s__%s__seed%d__%s.parquet", cohort, region, k, half)
        write_parquet(df, file.path(out_dir, fname))
        cat(sprintf("  seed%d %s: %d modules, %d samples, power=%d\n",
                    k, half, res$n_modules, length(idx), res$power))
    }
}
cat("**** WGCNA split-half complete ****\n")
