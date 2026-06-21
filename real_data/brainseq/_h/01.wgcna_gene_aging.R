#!/usr/bin/env Rscript
# Gene-level WGCNA aging baseline on BrainSEQ control regions (caudate,
# hippocampus, dlpfc). Parallels real_data/gtex/_h/02.wgcna_gene.R and provides
# the BrainSEQ WGCNA arm needed for cross-cohort replication against GTEx.
# Outputs (per region): modules.parquet, age_linear.parquet, age_spline.parquet
# in real_data/brainseq/<region>/_m/wgcna_gene/
#
# Differences from the GTEx baseline: BrainSEQ provides gene *counts* (converted
# here to log2(CPM+1), the count-scale analogue of GTEx's log2(TPM+1)), and the
# covariates include categorical terms (Sex, MoD) plus SNP PCs that may be NA for
# a few samples, so the spline design is built on complete-covariate samples only.
suppressPackageStartupMessages({
    library(arrow)
    library(WGCNA)
    library(splines)
    library(dplyr)
})

options(stringsAsFactors = FALSE)
WGCNA_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
if (is.na(WGCNA_THREADS) || WGCNA_THREADS < 2) WGCNA_THREADS <- 2L  # enableWGCNAThreads needs >=2
enableWGCNAThreads(nThreads = WGCNA_THREADS)

# ── Paths ──────────────────────────────────────────────────────────────────────
project_root <- here::here()
if (!file.exists(file.path(project_root, ".here"))) {
    stop("Cannot locate project root (.here file missing). Run from isograph-brain-aging-benchmarking/.")
}

BRAINSEQ_REGIONS <- c("caudate", "hippocampus", "dlpfc")

COVARIATE_COLS <- c("Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
                    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5")
AGE_COL <- "Age"
MIN_MODULE_SIZE <- 30
R2_THRESHOLD <- 0.85
AGE_PROBS <- c(0.10, 0.25, 0.50, 0.75, 0.90)
AGE_LABELS <- c("p10", "p25", "p50", "p75", "p90")

# ── Helpers ────────────────────────────────────────────────────────────────────
standardize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

fdr_bh <- function(p) p.adjust(p, method = "BH")

linear_age_assoc <- function(eigengenes, sample_tbl, age_col) {
    merged <- merge(eigengenes, sample_tbl[, c("sample_id", age_col)], by = "sample_id")
    age <- as.numeric(merged[[age_col]])
    module_cols <- setdiff(names(merged), c("sample_id", age_col))
    rows <- lapply(module_cols, function(col) {
        eg <- as.numeric(merged[[col]])
        ok <- is.finite(eg) & is.finite(age)
        if (sum(ok) < 10) return(NULL)
        ct <- cor.test(age[ok], eg[ok], method = "pearson")
        data.frame(module_id = col, trait = "Age_linear",
                   effect = ct$estimate, pvalue = ct$p.value, n = sum(ok),
                   stringsAsFactors = FALSE)
    })
    result <- do.call(rbind, Filter(Negate(is.null), rows))
    if (!is.null(result) && nrow(result) > 0) {
        result$fdr <- fdr_bh(result$pvalue)
    }
    result
}

spline_age_assoc <- function(eigengenes, sample_tbl, covariate_cols, age_col) {
    available_covs <- intersect(covariate_cols, names(sample_tbl))
    available_covs <- available_covs[available_covs != age_col]
    keep <- c("sample_id", age_col, available_covs)
    merged <- merge(eigengenes, sample_tbl[, keep], by = "sample_id")
    # Restrict to samples with complete age + covariates so the spline design and
    # the covariate model matrix stay row-aligned (BrainSEQ SNP PCs can be NA).
    merged <- merged[complete.cases(merged[, c(age_col, available_covs)]), ]
    if (nrow(merged) < 10) return(NULL)

    age_z <- standardize(as.numeric(merged[[age_col]]))
    knots <- quantile(age_z, c(1/3, 2/3))

    cov_df <- merged[, available_covs, drop = FALSE]
    cov_df <- model.matrix(~ ., data = cov_df)[, -1, drop = FALSE]

    B_obs <- ns(age_z, knots = knots, Boundary.knots = range(age_z))
    age_eval_z <- qnorm(AGE_PROBS)
    age_eval_clipped <- pmax(pmin(age_eval_z, max(age_z)), min(age_z))
    B_proj <- predict(B_obs, newx = age_eval_clipped)

    intercept <- matrix(1, nrow = nrow(merged), ncol = 1)
    X <- cbind(intercept, B_obs, cov_df)
    n_spline <- ncol(B_obs)
    spline_idx <- 2:(1 + n_spline)

    module_cols <- setdiff(names(eigengenes), "sample_id")
    rows <- lapply(module_cols, function(col) {
        y <- as.numeric(merged[[col]])
        ok <- is.finite(y) & apply(is.finite(X), 1, all)
        if (sum(ok) < 10) return(NULL)
        Xf <- X[ok, , drop = FALSE]
        yf <- y[ok]
        fit <- lm.fit(Xf, yf)
        beta <- coef(fit)
        df_res <- max(fit$rank - 1, 1)
        sigma2 <- sum(fit$residuals^2) / df_res
        V <- sigma2 * MASS::ginv(t(Xf) %*% Xf)
        beta_spline <- beta[spline_idx]
        V_spline <- V[spline_idx, spline_idx, drop = FALSE]
        beta_proj <- as.vector(B_proj %*% beta_spline)
        V_proj <- B_proj %*% V_spline %*% t(B_proj)
        se_proj <- sqrt(pmax(diag(V_proj), 0))
        z_proj <- ifelse(se_proj > 0, beta_proj / se_proj, 0)
        p_proj <- 2 * pnorm(-abs(z_proj))
        data.frame(
            module_id = col, trait = "Age_spline",
            age_label = AGE_LABELS, age_prob = AGE_PROBS,
            effect = beta_proj, se = se_proj, z = z_proj, pvalue = p_proj,
            n = sum(ok), stringsAsFactors = FALSE
        )
    })
    result <- do.call(rbind, Filter(Negate(is.null), rows))
    if (!is.null(result) && nrow(result) > 0) {
        # BH within each age_label across modules (ave keeps row alignment).
        result$fdr <- ave(result$pvalue, result$age_label, FUN = fdr_bh)
    }
    result
}

select_soft_power <- function(datExpr) {
    powers <- 1:20
    # Select for the signed topology blockwiseModules builds; pickSoftThreshold's
    # default unsigned network underestimates the power (collapsing all modules).
    sft <- pickSoftThreshold(datExpr, powerVector = powers,
                             networkType = "signed",
                             RsquaredCut = R2_THRESHOLD, verbose = 0)
    power <- sft$powerEstimate
    if (is.na(power)) {
        warning("Scale-free fit below threshold; defaulting to power=12 (signed)")
        power <- 12L
    }
    power <- max(as.integer(power), 12L)
    as.integer(power)
}

# ── Per-region runner ──────────────────────────────────────────────────────────
run_brainseq_wgcna <- function(region) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), region))

    proc_dir <- file.path(project_root, "inputs", "processed", "brainseq", region)
    bundle_dir <- file.path(project_root, "inputs", "bundles", "brainseq_v1", region)
    out_dir <- file.path(project_root, "real_data", "brainseq", region, "_m", "wgcna_gene")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Sample + gene lists from the bundle so WGCNA matches IsoGraph exactly.
    bundle_samples <- read_parquet(file.path(bundle_dir, "samples.parquet"))
    bundle_genes <- read_parquet(file.path(bundle_dir, "genes.parquet"))
    keep_ids <- bundle_samples$sample_id
    keep_gene_ids <- bundle_genes$gene_id

    # Gene counts -> bundle-filtered genes -> log2(CPM + 1).
    gene_counts <- read_parquet(file.path(proc_dir, "gene_counts.parquet"))
    gene_ids_avail <- intersect(keep_gene_ids, gene_counts$Geneid)
    if (length(gene_ids_avail) < MIN_MODULE_SIZE) {
        warning(sprintf("  %s: only %d matching bundle-filtered genes; skipping.", region, length(gene_ids_avail)))
        return(invisible(NULL))
    }
    gene_counts <- gene_counts[match(gene_ids_avail, gene_counts$Geneid), ]
    gene_ids <- gene_counts$Geneid
    sample_cols <- intersect(keep_ids, names(gene_counts))
    if (length(sample_cols) < 30) {
        warning(sprintf("  %s: only %d matching samples; skipping.", region, length(sample_cols)))
        return(invisible(NULL))
    }

    counts <- as.matrix(gene_counts[, sample_cols])
    rownames(counts) <- gene_ids
    lib_sizes <- colSums(counts)
    cpm <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
    expr_mat <- log2(cpm + 1)
    cat(sprintf("  Bundle expression filter retained %d genes for WGCNA\n", nrow(expr_mat)))

    # Samples × genes for WGCNA
    datExpr <- t(expr_mat)
    goodSamples <- goodSamplesGenes(datExpr, verbose = 0)
    datExpr <- datExpr[goodSamples$goodSamples, goodSamples$goodGenes]

    power <- select_soft_power(datExpr)

    net <- blockwiseModules(
        datExpr,
        power = power,
        minModuleSize = MIN_MODULE_SIZE,
        TOMType = "signed",
        networkType = "signed",
        numericLabels = FALSE,
        mergeCutHeight = 0.25,
        verbose = 0
    )

    # Build modules table (exclude grey = unassigned)
    colors <- net$colors
    gene_names <- colnames(datExpr)
    modules_df <- data.frame(
        gene_id = gene_names,
        module_id = colors,
        stringsAsFactors = FALSE
    )
    modules_df <- modules_df[modules_df$module_id != "grey", ]

    # Size-order module labels to M000, M001, ...
    mod_sizes <- sort(table(modules_df$module_id), decreasing = TRUE)
    label_map <- setNames(
        sprintf("M%03d", seq_along(mod_sizes) - 1),
        names(mod_sizes)
    )
    modules_df$module_id <- label_map[modules_df$module_id]

    write_parquet(modules_df, file.path(out_dir, "modules.parquet"))

    # Eigengenes
    me <- moduleEigengenes(datExpr, colors = colors)$eigengenes
    me <- me[, names(me) != "MEgrey", drop = FALSE]
    colnames(me) <- label_map[sub("^ME", "", colnames(me))]
    me_df <- cbind(data.frame(sample_id = rownames(me), stringsAsFactors = FALSE), me)

    sample_tbl <- bundle_samples[bundle_samples$sample_id %in% rownames(me), ]

    linear <- linear_age_assoc(me_df, sample_tbl, AGE_COL)
    if (!is.null(linear)) write_parquet(linear, file.path(out_dir, "age_linear.parquet"))

    spline <- spline_age_assoc(me_df, sample_tbl, COVARIATE_COLS, AGE_COL)
    if (!is.null(spline)) write_parquet(spline, file.path(out_dir, "age_spline.parquet"))

    n_modules <- length(unique(modules_df$module_id))
    cat(sprintf("  -> %d modules | %d samples | power=%d\n", n_modules, nrow(datExpr), power))
}

# ── Main ───────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
regions_to_run <- if (length(args) > 0) args else BRAINSEQ_REGIONS
for (region in regions_to_run) {
    tryCatch(
        run_brainseq_wgcna(region),
        error = function(e) message(sprintf("ERROR in %s: %s", region, conditionMessage(e)))
    )
}
