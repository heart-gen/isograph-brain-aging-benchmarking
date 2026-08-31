#!/usr/bin/env Rscript
## Per-gene sQTL-vs-eQTL colocalization contrast, step 2: coloc.abf on GTEx v11
## cis all-pairs for ONE brain tissue.
##
## For every (analysis, locus, gene) target from coloc_modality_contrast.py --stage prep,
## this extracts the gene's full cis nominal statistics in this tissue for BOTH
## modalities, joins them to the per-locus GWAS summary statistics already built for the
## CLPP layer, and runs coloc::coloc.abf on each.
##
## The two arms are matched by construction: same donors (GTEx calls sQTL and eQTL in
## the same samples), same GWAS, same shared-variant set, same gene. So the module
## membership that selected the gene cancels, and the contrast is a within-gene one.
##
## sQTL phenotype choice is GTEx's own grouped-permutation representative intron
## (one row per gene in *.sGenes.txt.gz, chosen by group-corrected pval_beta). That
## choice never sees the GWAS, so it cannot be circular; using the best intron BY
## COLOCALIZATION would have handed the sQTL arm a free maximum over ~7 tests.
##
## coloc.abf is LD-free and its Bayes factors depend on z^2, so it needs no LD panel and
## is invariant to allele orientation -- the hg19/hg38 mismatch that constrains the CLPP
## path does not apply. Alleles are harmonized anyway, but only to report the DIRECTION
## of effect at the top shared SNP, which is not part of the primary test.
##
## Usage:  Rscript 05_genetic_anchoring/_h/20.coloc_modality_abf.R <Brain_Tissue>
## Output: 05_genetic_anchoring/_m/coloc_modality_contrast/abf/<tissue>.parquet
suppressPackageStartupMessages({
    library(data.table); library(arrow); library(coloc)
})

args   <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: 20.coloc_modality_abf.R <Brain_Tissue> [chr]")
TISSUE <- args[1]
## Optional single chromosome, for validating a change or re-running one chromosome
## without redoing the tissue. Writes to <tissue>.chr<N>.parquet so a scoped run can
## never be mistaken for, or overwrite, a whole-tissue result.
ONE_CHR <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_
ROOT   <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
GTEX   <- "/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11"
MDIR   <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc_modality_contrast")
BRIDGE <- file.path(ROOT, "inputs", "raw", "gtex_v11", "variant_bridge")
OUTD   <- file.path(MDIR, "abf"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

P12       <- c(1e-5, 5e-6, 1e-6)   # 1e-5 = coloc default = primary
MIN_SHARED <- 50L                  # emit the row; the meta stage applies the real floor
SDY       <- 1                     # GTEx phenotypes are inverse-normal transformed

## GTEx donors in this tissue, from the eQTL covariates header (columns - 1). BOTH
## modalities are called in exactly these samples -- that shared sample set is what
## makes the sQTL and eQTL arms a paired comparison rather than two separate studies.
COVF  <- file.path(GTEX, "GTEx_Analysis_v11_eQTL_covariates",
                   paste0(TISSUE, ".v11.covariates.txt"))
if (!file.exists(COVF)) stop("missing GTEx covariates file: ", COVF)
N_QTL <- length(strsplit(readLines(COVF, n = 1L), "\t")[[1]]) - 1L
message(sprintf("[%s] GTEx donors: %d", TISSUE, N_QTL))

targets <- as.data.table(read_parquet(file.path(MDIR, "targets.parquet")))
srep    <- as.data.table(read_parquet(file.path(MDIR, "sqtl_representative.parquet")))
srep    <- srep[tissue == TISSUE, .(gene, phenotype_id)]
setkey(srep, gene)
message(sprintf("[%s] %d targets, %d sQTL representative phenotypes",
                TISSUE, nrow(targets), nrow(srep)))

## ---- per-locus GWAS, cached once per analysis -------------------------------
## Same files and the same long-range-LD exclusions as 10.coloc_clpp.R, so the two
## estimators are computed on the same GWAS input.
gwas_cache <- new.env(parent = emptyenv())
get_gwas <- function(analysis, lid, chr) {
    key <- paste(analysis, lid)
    if (!is.null(gwas_cache[[key]])) return(gwas_cache[[key]])
    f <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc", analysis, "susie",
                   paste0(lid, ".gwas.tsv"))
    if (!file.exists(f)) { gwas_cache[[key]] <- data.table(); return(data.table()) }
    g <- fread(f)
    ef <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc", analysis,
                    "exclude_regions.tsv")
    if (file.exists(ef)) {
        ex <- fread(ef)
        for (j in seq_len(nrow(ex)))
            if (chr == ex$chr[j]) g <- g[!(pos >= ex$start[j] & pos <= ex$stop[j])]
    }
    g <- g[!is.na(beta) & !is.na(se) & se > 0]
    g <- g[!duplicated(rsid)]           # multi-allelic rsIDs are unusable as keys
    g <- g[, .(rsid, gwas_beta = beta, gwas_varbeta = se^2,
               gwas_a1 = toupper(a1), gwas_a2 = toupper(a2), gwas_p = p)]
    gwas_cache[[key]] <- g
    g
}

## ---- coloc for one (gene, locus, modality) ----------------------------------
## Returns one row per p12 prior. Allele orientation does not enter the posterior
## (coloc.abf uses z^2); it is used only for the reported direction.
run_one <- function(q, g, meta) {
    m <- merge(q, g, by = "rsid")
    m <- m[!is.na(qtl_beta) & !is.na(qtl_varbeta) & qtl_varbeta > 0]
    m <- m[!duplicated(rsid)]
    if (nrow(m) < MIN_SHARED) return(NULL)
    d_gwas <- list(beta = m$gwas_beta, varbeta = m$gwas_varbeta, snp = m$rsid,
                   type = "cc", s = meta$s, N = meta$n_gwas)
    d_qtl  <- list(beta = m$qtl_beta, varbeta = m$qtl_varbeta, snp = m$rsid,
                   type = "quant", sdY = SDY, N = meta$n_qtl)
    out <- vector("list", length(P12))
    for (i in seq_along(P12)) {
        ## coloc.abf reports each fit on stdout and warns on every dataset check.
        ## At ~45k fits that is multi-GB of log for no information, so both are
        ## captured; a genuine failure still surfaces as a NULL and a skipped row.
        r <- tryCatch({
                 invisible(capture.output(
                     z <- suppressWarnings(suppressMessages(
                         coloc.abf(d_gwas, d_qtl, p12 = P12[i])))))
                 z
             }, error = function(e) NULL)
        if (is.null(r)) next
        s <- as.list(r$summary)
        ## direction at the SNP carrying the most H4 posterior
        top <- as.data.table(r$results)[which.max(SNP.PP.H4)]
        tm  <- m[rsid == top$snp][1]
        flip <- NA_integer_
        if (nrow(tm)) {
            if (!is.na(tm$alt) && tm$gwas_a1 == tm$alt) flip <- 1L
            else if (!is.na(tm$ref) && tm$gwas_a1 == tm$ref) flip <- -1L
        }
        out[[i]] <- data.table(
            p12 = P12[i], nsnps = as.integer(s$nsnps),
            PP0 = s$PP.H0.abf, PP1 = s$PP.H1.abf, PP2 = s$PP.H2.abf,
            PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf,
            top_rsid = top$snp, top_pp4 = top$SNP.PP.H4,
            n_qtl_donors = meta$n_qtl,
            qtl_min_p = min(m$qtl_p, na.rm = TRUE),
            gwas_min_p = min(m$gwas_p, na.rm = TRUE),
            dir_concordant = if (is.na(flip)) NA
                             else sign(tm$gwas_beta) * sign(tm$qtl_beta) * flip > 0)
    }
    rbindlist(out)
}

## GWAS N / case fraction, derived from the sumstats by the prep stage rather than
## hardcoded here. Recorded for provenance only: with beta+varbeta supplied for a
## case-control dataset, coloc.abf's Bayes factors do not depend on N or s (asserted in
## tests/test_coloc_modality_contrast.py). They are passed so the dataset validates.
gm <- fread(file.path(MDIR, "gwas_meta.tsv"))
TRAIT_NS <- setNames(
    lapply(seq_len(nrow(gm)), function(i) list(s = gm$s[i], n_gwas = gm$n_gwas[i])),
    gm$trait)

results <- list()
chrs <- if (is.na(ONE_CHR)) sort(unique(targets$chr)) else ONE_CHR
for (CH in chrs) {
    tg <- targets[chr == CH]
    genes <- unique(tg$gene)
    ef <- file.path(GTEX, "GTEx_Analysis_v11_eQTL_all_associations",
                    sprintf("%s.v11.allpairs.chr%d.parquet", TISSUE, CH))
    sf <- file.path(GTEX, "GTEx_Analysis_v11_sQTL_all_associations",
                    sprintf("%s.v11.cis_sqtl.allpairs.chr%d.parquet", TISSUE, CH))
    if (!file.exists(ef) || !file.exists(sf)) { message("  chr", CH, ": missing all-pairs"); next }

    ## eQTL: resolve the versioned gene_id by bare match, then a single pushdown read
    ed <- open_dataset(ef)
    gid <- ed |> dplyr::distinct(gene_id) |> dplyr::collect() |> data.table::as.data.table()
    gid[, bare := sub("\\..*$", "", gene_id)]
    want_e <- gid[bare %in% genes, gene_id]
    eq <- if (length(want_e)) {
        as.data.table(ed |> dplyr::filter(gene_id %in% want_e) |>
            dplyr::select(gene_id, variant_id, pval_nominal, slope, slope_se) |>
            dplyr::collect())
    } else data.table()
    if (nrow(eq)) eq[, gene := sub("\\..*$", "", gene_id)][, gene_id := NULL]

    ## sQTL: GTEx's representative phenotype for each candidate gene
    want_s <- srep[gene %in% genes]
    sq <- if (nrow(want_s)) {
        as.data.table(open_dataset(sf) |>
            dplyr::filter(phenotype_id %in% want_s$phenotype_id) |>
            dplyr::select(phenotype_id, variant_id, pval_nominal, slope, slope_se) |>
            dplyr::collect())
    } else data.table()
    if (nrow(sq)) sq <- merge(sq, want_s, by = "phenotype_id")[, phenotype_id := NULL]

    ## variant_id (b38) -> rsID, and ref/alt for the direction report
    bf <- file.path(BRIDGE, sprintf("chr%d.parquet", CH))
    if (!file.exists(bf)) stop("missing variant bridge: ", bf,
                               " (run --stage bridge first)")
    br <- as.data.table(read_parquet(bf))
    setkey(br, variant_id)

    for (MOD in c("eQTL", "sQTL")) {
        q0 <- if (MOD == "eQTL") eq else sq
        if (!nrow(q0)) next
        q0 <- merge(q0, br, by = "variant_id")
        parts <- tstrsplit(q0$variant_id, "_", fixed = TRUE)
        q0[, `:=`(ref = parts[[3]], alt = parts[[4]],
                  qtl_beta = slope, qtl_varbeta = slope_se^2, qtl_p = pval_nominal)]
        q0 <- q0[, .(gene, rsid, ref, alt, qtl_beta, qtl_varbeta, qtl_p)]
        setkey(q0, gene)
        for (i in seq_len(nrow(tg))) {
            tr <- tg[i]
            meta <- TRAIT_NS[[tr$trait]]
            if (is.null(meta)) next
            meta$n_qtl <- N_QTL
            qg <- q0[.(tr$gene), nomatch = 0L]
            if (!nrow(qg)) next
            g <- get_gwas(tr$analysis, tr$LOCUS_ID, CH)
            if (!nrow(g)) next
            r <- run_one(qg, g, meta)
            if (is.null(r) || !nrow(r)) next
            results[[length(results) + 1]] <- cbind(
                data.table(analysis = tr$analysis, trait = tr$trait,
                           LOCUS_ID = tr$LOCUS_ID, chr = CH, gene = tr$gene,
                           symbol = tr$symbol, module_id = tr$module_id,
                           go_invisible = tr$go_invisible, tissue = TISSUE,
                           modality = MOD), r)
        }
    }
    message(sprintf("  chr%-2d  genes=%3d  eQTL_rows=%8d  sQTL_rows=%8d  cum_results=%d",
                    CH, length(genes), nrow(eq), nrow(sq), length(results)))
    rm(eq, sq, br); gc(verbose = FALSE)
}

if (!length(results)) { message("No coloc rows produced for ", TISSUE); quit(status = 0) }
res <- rbindlist(results, fill = TRUE)
stem <- if (is.na(ONE_CHR)) TISSUE else sprintf("%s.chr%d", TISSUE, ONE_CHR)
write_parquet(res, file.path(OUTD, paste0(stem, ".parquet")))
message(sprintf("\n[%s] %d rows -> %s", stem, nrow(res),
                file.path(OUTD, paste0(stem, ".parquet"))))
prim <- res[p12 == 1e-5 & nsnps >= 100]
message(sprintf("  primary-prior cells: %d | sQTL PP4>=0.8: %d | eQTL PP4>=0.8: %d",
                nrow(prim), nrow(prim[modality == "sQTL" & PP4 >= 0.8]),
                nrow(prim[modality == "eQTL" & PP4 >= 0.8])))
