#!/usr/bin/env Rscript
# Gene-level WGCNA on BrainSEQ caudate SCZD+Control bundle.
# Outputs: real_data/brainseq/caudate_sczd/_m/wgcna_gene/
#   modules.parquet, age_linear.parquet, age_spline.parquet
suppressPackageStartupMessages({
    library(arrow)
    library(WGCNA)
    library(splines)
    library(dplyr)
    library(MASS)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

# ── Paths ──────────────────────────────────────────────────────────────────────
script_dir <- dirname(normalizePath(if (interactive()) getwd() else commandArgs()[4], mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)

COVARIATE_COLS <- c("Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
                    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5")
AGE_COL <- "Age"
MIN_MODULE_SIZE <- 30
R2_THRESHOLD <- 0.85
AGE_PROBS <- c(0.10, 0.25, 0.50, 0.75, 0.90)
AGE_LABELS <- c("p10", "p25", "p50", "p75", "p90")

# ── Helpers (mirror run_models.py) ─────────────────────────────────────────────
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
    if (!is.null(result) && nrow(result) > 0) result$fdr <- fdr_bh(result$pvalue)
    result
}

spline_age_assoc <- function(eigengenes, sample_tbl, covariate_cols, age_col) {
    keep <- c("sample_id", age_col, intersect(covariate_cols, names(sample_tbl)))
    merged <- merge(eigengenes, sample_tbl[, keep], by = "sample_id")
    merged <- merged[!is.na(merged[[age_col]]), ]

    age_z <- standardize(as.numeric(merged[[age_col]]))
    knots <- quantile(age_z, c(1/3, 2/3))

    avail_covs <- intersect(covariate_cols, names(merged))
    avail_covs <- avail_covs[avail_covs != age_col]
    cov_mat <- model.matrix(~ ., data = merged[, avail_covs, drop = FALSE])[, -1, drop = FALSE]

    B_obs <- ns(age_z, knots = knots, Boundary.knots = range(age_z))
    age_eval_z <- qnorm(AGE_PROBS)
    age_eval_clipped <- pmax(pmin(age_eval_z, max(age_z)), min(age_z))
    B_proj <- predict(B_obs, newx = age_eval_clipped)

    X <- cbind(1, B_obs, cov_mat)
    n_spline <- ncol(B_obs)
    spline_idx <- 2:(1 + n_spline)

    module_cols <- setdiff(names(eigengenes), "sample_id")
    rows <- lapply(module_cols, function(col) {
        y <- as.numeric(merged[[col]])
        ok <- is.finite(y) & apply(is.finite(X), 1, all)
        if (sum(ok) < 10) return(NULL)
        Xf <- X[ok, , drop = FALSE]; yf <- y[ok]
        fit <- lm.fit(Xf, yf)
        df_res <- max(fit$rank - 1, 1)
        sigma2 <- sum(fit$residuals^2) / df_res
        V <- sigma2 * ginv(t(Xf) %*% Xf)
        beta_spline <- coef(fit)[spline_idx]
        V_spline <- V[spline_idx, spline_idx, drop = FALSE]
        beta_proj <- as.vector(B_proj %*% beta_spline)
        se_proj <- sqrt(pmax(diag(B_proj %*% V_spline %*% t(B_proj)), 0))
        z_proj <- ifelse(se_proj > 0, beta_proj / se_proj, 0)
        p_proj <- 2 * pnorm(-abs(z_proj))
        data.frame(module_id = col, trait = "Age_spline",
                   age_label = AGE_LABELS, age_prob = AGE_PROBS,
                   effect = beta_proj, se = se_proj, z = z_proj,
                   pvalue = p_proj, n = sum(ok), stringsAsFactors = FALSE)
    })
    result <- do.call(rbind, Filter(Negate(is.null), rows))
    if (!is.null(result) && nrow(result) > 0) {
        fdr_by_label <- tapply(result$pvalue, result$age_label, fdr_bh)
        result$fdr <- unlist(fdr_by_label[result$age_label])
    }
    result
}

# ── Main ───────────────────────────────────────────────────────────────────────
bundle_dir <- file.path(project_root, "inputs", "bundles", "brainseq_sczd", "caudate")
out_dir <- file.path(project_root, "real_data", "brainseq", "caudate_sczd", "_m", "wgcna_gene")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load gene counts from processed data (same source as bundle builder)
gene_counts_path <- file.path(project_root, "inputs", "processed", "brainseq", "caudate", "gene_counts.parquet")
gene_df <- read_parquet(gene_counts_path)

# Load sample metadata from bundle
sample_tbl <- read_parquet(file.path(bundle_dir, "samples.parquet"))
keep_ids <- sample_tbl$sample_id

# Subset to bundle samples (same filter as IsoGraph run)
id_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
id_cols_present <- intersect(id_cols, names(gene_df))
sample_cols_avail <- intersect(keep_ids, names(gene_df))
if (length(sample_cols_avail) < 30) stop("Too few samples found in gene counts matrix.")

gene_ids <- gene_df$Geneid
expr_mat <- as.matrix(gene_df[, sample_cols_avail])
rownames(expr_mat) <- gene_ids

# Library-size normalize and log2 transform
lib_sizes <- colSums(expr_mat)
expr_mat <- sweep(expr_mat, 2, lib_sizes / 1e6, "/")  # CPM
expr_mat <- log2(expr_mat + 1)

# Filter to top 50% genes by mean expression
gene_means <- rowMeans(expr_mat)
expr_mat <- expr_mat[gene_means >= median(gene_means), ]

datExpr <- t(expr_mat)
good <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[good$goodSamples, good$goodGenes]
cat(sprintf("Samples: %d | Genes: %d\n", nrow(datExpr), ncol(datExpr)))

# Soft-threshold selection
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = R2_THRESHOLD, verbose = 0)
power <- sft$powerEstimate
if (is.na(power)) { warning("R2 threshold not met; defaulting to 6"); power <- 6L }
cat(sprintf("Soft-threshold power: %d\n", power))

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

colors <- net$colors
modules_df <- data.frame(gene_id = colnames(datExpr), module_id = colors,
                         stringsAsFactors = FALSE)
modules_df <- modules_df[modules_df$module_id != "grey", ]

mod_sizes <- sort(table(modules_df$module_id), decreasing = TRUE)
label_map <- setNames(sprintf("M%03d", seq_along(mod_sizes) - 1), names(mod_sizes))
modules_df$module_id <- label_map[modules_df$module_id]
write_parquet(modules_df, file.path(out_dir, "modules.parquet"))

# Eigengenes
me <- moduleEigengenes(datExpr, colors = colors)$eigengenes
me <- me[, names(me) != "MEgrey", drop = FALSE]
colnames(me) <- label_map[sub("^ME", "", colnames(me))]
me_df <- cbind(data.frame(sample_id = rownames(me), stringsAsFactors = FALSE), me)

# Restrict sample table to samples that survived WGCNA QC
sample_tbl_use <- sample_tbl[sample_tbl$sample_id %in% rownames(me), ]

linear <- linear_age_assoc(me_df, sample_tbl_use, AGE_COL)
if (!is.null(linear)) write_parquet(linear, file.path(out_dir, "age_linear.parquet"))

spline <- spline_age_assoc(me_df, sample_tbl_use, COVARIATE_COLS, AGE_COL)
if (!is.null(spline)) write_parquet(spline, file.path(out_dir, "age_spline.parquet"))

cat(sprintf("Done: %d modules | power=%d → %s\n",
    length(unique(modules_df$module_id)), power, out_dir))
