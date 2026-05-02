#!/usr/bin/env Rscript
# Gene-level WGCNA on BrainSEQ caudate SCZD+Control bundle.
# Outputs: real_data/brainseq/caudate_sczd/_m/wgcna_gene/
#   modules.parquet, diagnosis_assoc.parquet
suppressPackageStartupMessages({
    library(data.table)
    library(WGCNA)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

# ── Paths ──────────────────────────────────────────────────────────────────────
script_dir <- dirname(normalizePath(if (interactive()) getwd() else commandArgs()[4], mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)

COVARIATE_COLS <- c("Sex", "MoD", "RIN", "mapping_rate", "mito_rate",
                    "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5")
DIAGNOSIS_COVARIATE_COLS <- c("Age", COVARIATE_COLS)
DX_COL <- "Dx"
CONTROL_LABEL <- "Control"
CASE_LABEL <- "SCZD"
MIN_MODULE_SIZE <- 30
R2_THRESHOLD <- 0.85
PYTHON_BIN <- path.expand(Sys.getenv(
    "ISOGRAPH_PYTHON",
    unset = file.path(Sys.getenv("HOME"), ".venvs", "isograph", "bin", "python")
))

# ── Helpers ───────────────────────────────────────────────────────────────────
fdr_bh <- function(p) p.adjust(p, method = "BH")

run_python <- function(args) {
    out <- system2(PYTHON_BIN, shQuote(args), stdout = TRUE, stderr = TRUE)
    status <- attr(out, "status")
    if (is.null(status)) status <- 0L
    if (!identical(status, 0L)) {
        stop(sprintf("Python command failed (%s):\n%s", PYTHON_BIN, paste(out, collapse = "\n")))
    }
    out
}

check_pyarrow <- function() {
    if (!file.exists(PYTHON_BIN) && !nzchar(Sys.which(PYTHON_BIN))) {
        stop(sprintf("Python executable not found: %s", PYTHON_BIN))
    }
    run_python(c("-c", "import pandas, pyarrow, pyarrow.parquet"))
    invisible(TRUE)
}

read_parquet_table <- function(path) {
    if (requireNamespace("arrow", quietly = TRUE)) {
        result <- tryCatch(
            arrow::read_parquet(path),
            error = function(e) {
                warning(sprintf("R arrow could not read %s (%s); trying Python pyarrow.", path, conditionMessage(e)))
                NULL
            }
        )
        if (!is.null(result)) return(as.data.frame(result, check.names = FALSE))
    }

    check_pyarrow()
    tmp <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    code <- paste(
        "import sys, pyarrow.parquet as pq",
        "pq.read_table(sys.argv[1]).to_pandas().to_csv(sys.argv[2], index=False)",
        sep = "; "
    )
    run_python(c("-c", code, path, tmp))
    fread(tmp, data.table = FALSE, check.names = FALSE)
}

write_parquet_table <- function(df, path) {
    if (requireNamespace("arrow", quietly = TRUE)) {
        ok <- tryCatch(
            {
                arrow::write_parquet(df, path, compression = "uncompressed")
                TRUE
            },
            error = function(e) {
                warning(sprintf("R arrow could not write %s (%s); trying Python pyarrow.", path, conditionMessage(e)))
                FALSE
            }
        )
        if (ok) return(invisible(path))
    }

    check_pyarrow()
    tmp <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    fwrite(df, tmp)
    code <- paste(
        "import sys, pandas as pd, pyarrow as pa, pyarrow.parquet as pq",
        "df = pd.read_csv(sys.argv[1])",
        "pq.write_table(pa.Table.from_pandas(df, preserve_index=False), sys.argv[2], compression=None)",
        sep = "; "
    )
    run_python(c("-c", code, tmp, path))
    invisible(path)
}

read_gene_counts <- function(parquet_path, raw_tsv_path) {
    if (requireNamespace("arrow", quietly = TRUE)) {
        result <- tryCatch(
            arrow::read_parquet(parquet_path),
            error = function(e) {
                warning(sprintf(
                    "R arrow could not read %s (%s); falling back to %s",
                    parquet_path, conditionMessage(e), raw_tsv_path
                ))
                NULL
            }
        )
        if (!is.null(result)) return(as.data.frame(result, check.names = FALSE))
    }
    fread(raw_tsv_path, data.table = FALSE, check.names = FALSE)
}

diagnosis_assoc <- function(eigengenes, sample_tbl, covariate_cols, dx_col,
                            control_label, case_label) {
    if (!dx_col %in% names(sample_tbl)) stop(sprintf("Missing diagnosis column: %s", dx_col))
    keep <- c("sample_id", dx_col, intersect(covariate_cols, names(sample_tbl)))
    merged <- merge(eigengenes, sample_tbl[, keep], by = "sample_id")
    merged <- merged[merged[[dx_col]] %in% c(control_label, case_label), ]
    merged[[dx_col]] <- factor(merged[[dx_col]], levels = c(control_label, case_label))

    if (length(unique(na.omit(merged[[dx_col]]))) < 2) {
        stop(sprintf("Diagnosis association requires both %s and %s samples.", control_label, case_label))
    }

    avail_covs <- intersect(covariate_cols, names(merged))
    module_cols <- setdiff(names(eigengenes), "sample_id")
    rows <- lapply(module_cols, function(col) {
        analysis <- data.frame(
            eigengene = as.numeric(merged[[col]]),
            Diagnosis = merged[[dx_col]],
            merged[, avail_covs, drop = FALSE],
            check.names = FALSE
        )
        analysis[] <- lapply(analysis, function(x) {
            if (is.numeric(x)) x[!is.finite(x)] <- NA_real_
            x
        })
        analysis <- na.omit(analysis)
        if (nrow(analysis) < 10 || length(unique(analysis$Diagnosis)) < 2) return(NULL)

        fit <- lm(reformulate(c("Diagnosis", avail_covs), response = "eigengene"),
                  data = analysis)
        coef_name <- sprintf("Diagnosis%s", case_label)
        coef_tbl <- coef(summary(fit))
        if (!coef_name %in% rownames(coef_tbl)) return(NULL)
        used <- model.frame(fit)
        data.frame(
            module_id = col,
            trait = sprintf("%s_%s_vs_%s", dx_col, case_label, control_label),
            effect = coef_tbl[coef_name, "Estimate"],
            se = coef_tbl[coef_name, "Std. Error"],
            t = coef_tbl[coef_name, "t value"],
            pvalue = coef_tbl[coef_name, "Pr(>|t|)"],
            n = nrow(used),
            n_control = sum(used$Diagnosis == control_label),
            n_case = sum(used$Diagnosis == case_label),
            stringsAsFactors = FALSE
        )
    })
    result <- do.call(rbind, Filter(Negate(is.null), rows))
    if (!is.null(result) && nrow(result) > 0) result$fdr <- fdr_bh(result$pvalue)
    result
}

# ── Main ───────────────────────────────────────────────────────────────────────
bundle_dir <- file.path(project_root, "inputs", "bundles", "brainseq_sczd", "caudate")
out_dir <- file.path(project_root, "real_data", "brainseq", "caudate_sczd", "_m", "wgcna_gene")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load gene counts from processed data, with raw TSV fallback for R-arrow builds
# that cannot read the processed parquet codec.
gene_counts_path <- file.path(project_root, "inputs", "processed", "brainseq", "caudate", "gene_counts.parquet")
gene_counts_raw_path <- file.path(project_root, "inputs", "raw", "brainseq", "counts", "caudate", "gene-counts.tsv")
gene_df <- read_gene_counts(gene_counts_path, gene_counts_raw_path)

# Load sample metadata from bundle
sample_tbl <- read_parquet_table(file.path(bundle_dir, "samples.parquet"))
gene_tbl <- read_parquet_table(file.path(bundle_dir, "genes.parquet"))
keep_ids <- sample_tbl$sample_id
keep_gene_ids <- gene_tbl$gene_id

# Subset to bundle samples (same filter as IsoGraph run)
id_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
id_cols_present <- intersect(id_cols, names(gene_df))
sample_cols_avail <- intersect(keep_ids, names(gene_df))
if (length(sample_cols_avail) < 30) stop("Too few samples found in gene counts matrix.")
gene_ids_avail <- intersect(keep_gene_ids, gene_df$Geneid)
if (length(gene_ids_avail) < MIN_MODULE_SIZE) stop("Too few bundle-filtered genes found in gene counts matrix.")
gene_df <- gene_df[match(gene_ids_avail, gene_df$Geneid), ]

gene_ids <- gene_df$Geneid
expr_mat <- as.matrix(gene_df[, sample_cols_avail])
rownames(expr_mat) <- gene_ids

# Library-size normalize and log2 transform. Low-expression genes were removed
# centrally while building the BrainSEQ SCZD bundle.
lib_sizes <- colSums(expr_mat)
cpm_mat <- sweep(expr_mat, 2, lib_sizes / 1e6, "/")
expr_mat <- log2(cpm_mat + 1)
cat(sprintf("Bundle expression filter retained %d genes for WGCNA\n", nrow(expr_mat)))

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
write_parquet_table(modules_df, file.path(out_dir, "modules.parquet"))

# Eigengenes
me <- moduleEigengenes(datExpr, colors = colors)$eigengenes
me <- me[, names(me) != "MEgrey", drop = FALSE]
colnames(me) <- label_map[sub("^ME", "", colnames(me))]
me_df <- cbind(data.frame(sample_id = rownames(me), stringsAsFactors = FALSE), me)

# Restrict sample table to samples that survived WGCNA QC
sample_tbl_use <- sample_tbl[sample_tbl$sample_id %in% rownames(me), ]

diagnosis <- diagnosis_assoc(me_df, sample_tbl_use, DIAGNOSIS_COVARIATE_COLS,
                             DX_COL, CONTROL_LABEL, CASE_LABEL)
if (!is.null(diagnosis)) {
    write_parquet_table(diagnosis, file.path(out_dir, "diagnosis_assoc.parquet"))
}

cat(sprintf("Done: %d modules | power=%d → %s\n",
    length(unique(modules_df$module_id)), power, out_dir))
