#!/usr/bin/env Rscript
# satuRn DTU switch test for the IsoGraph concordance CLI (isa_concordance.py).
# satuRn is the DTU test IsoformSwitchAnalyzeR v2 runs internally
# (isoformSwitchTestSatuRn); this is that test used exactly as published -- its
# empirical-null-recalibrated FDR (Gilis et al. 2022) is what makes the call honest.
#
# A single-coefficient contrast is tested via satuRn::testDTU (the only route that
# gets satuRn's empirical-null FDR; a multi-df spline omnibus would bypass it):
#   --exposure-kind numeric : formula ~ <col> + covs, test the <col> coefficient.
#                             Used for --trait age with a continuous z-scored age.
#   --exposure-kind factor  : formula ~ <col> + covs with <col> releveled to
#                             --ref-level, test the <test-coef> coefficient.
#                             Used for --trait dx (SCZD vs Control).
# Inputs are staged by the Python driver. NB satuRn's vcovUnsc is UNSCALED (true
# cov = dispersion * vcovUnsc, dispersion ~ O(10)); never Wald-test it directly.
suppressPackageStartupMessages({
  library(satuRn)
  library(SummarizedExperiment)
})

# ---- minimal arg parsing (avoid an optparse dependency) --------------------- #
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1L && i < length(args)) args[i + 1L] else default
}
counts_path   <- get_arg("--counts")
coldata_path  <- get_arg("--coldata")
txinfo_path   <- get_arg("--txinfo")
covariates    <- get_arg("--covariates", "")
exposure_col  <- get_arg("--exposure-col")
exposure_kind <- get_arg("--exposure-kind", "numeric")  # "numeric" | "factor"
ref_level     <- get_arg("--ref-level")                 # factor only
test_coef     <- get_arg("--test-coef")                 # design column to test
out_path      <- get_arg("--out")
seed          <- as.integer(get_arg("--seed", "13"))
cores         <- as.integer(get_arg("--cores", "1"))
set.seed(seed)

read_tsv <- function(p) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(p, sep = "\t", header = TRUE, data.table = FALSE))
  } else {
    read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# ---- load staged inputs ----------------------------------------------------- #
counts <- read_tsv(counts_path)
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

coldata <- read_tsv(coldata_path)
rownames(coldata) <- as.character(coldata$sample_id)

txinfo <- read_tsv(txinfo_path)
rownames(txinfo) <- txinfo$isoform_id

counts <- counts[rownames(txinfo), rownames(coldata), drop = FALSE]  # align

covs <- if (nchar(covariates) > 0) strsplit(covariates, ",")[[1]] else character(0)
for (cv in covs) {
  x <- coldata[[cv]]
  coldata[[cv]] <- if (is.character(x) || is.logical(x)) factor(x) else as.numeric(x)
}
if (identical(exposure_kind, "factor")) {
  coldata[[exposure_col]] <- relevel(factor(coldata[[exposure_col]]), ref = ref_level)
} else {
  coldata[[exposure_col]] <- as.numeric(coldata[[exposure_col]])
}

form <- as.formula(
  paste("~", exposure_col,
        if (length(covs)) paste("+", paste(covs, collapse = " + ")) else "")
)

# ---- build SummarizedExperiment + fit + test a single coefficient ----------- #
sumExp <- SummarizedExperiment(
  assays  = list(counts = counts),
  colData = S4Vectors::DataFrame(coldata),
  rowData = S4Vectors::DataFrame(txinfo)
)
metadata(sumExp)$formula <- form

bpparam <- if (cores > 1L) BiocParallel::MulticoreParam(cores) else BiocParallel::SerialParam()
sumExp <- satuRn::fitDTU(object = sumExp, formula = form,
                         parallel = (cores > 1L), BPPARAM = bpparam, verbose = TRUE)

design <- model.matrix(form, data = as.data.frame(colData(sumExp)))
stopifnot(test_coef %in% colnames(design))
L <- matrix(0, nrow = ncol(design), ncol = 1,
            dimnames = list(colnames(design), "Contrast"))
L[test_coef, 1] <- 1

sumExp <- tryCatch(
  satuRn::testDTU(object = sumExp, contrasts = L,
                  diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE),
  error = function(e) satuRn::testDTU(object = sumExp, contrasts = L, sort = FALSE)
)

res <- as.data.frame(rowData(sumExp)[["fitDTUResult_Contrast"]])
res$isoform_id <- rownames(res)
res$gene_id <- txinfo[res$isoform_id, "gene_id"]

pick <- function(df, name) if (name %in% names(df)) df[[name]] else NA_real_
out <- data.frame(
  isoform_id     = res$isoform_id,
  gene_id        = res$gene_id,
  estimate       = pick(res, "estimates"),
  pval           = if ("pval" %in% names(res)) res$pval else pick(res, "p.value"),
  regular_FDR    = pick(res, "regular_FDR"),
  empirical_pval = pick(res, "empirical_pval"),
  empirical_FDR  = pick(res, "empirical_FDR"),
  stringsAsFactors = FALSE
)

con <- gzfile(out_path, "wt")
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat(sprintf("[satuRn] tested coef '%s'; wrote %d transcript rows (%d genes) -> %s\n",
            test_coef, nrow(out), length(unique(out$gene_id)), out_path))
