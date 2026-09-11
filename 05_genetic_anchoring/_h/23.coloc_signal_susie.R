#!/usr/bin/env Rscript
## Signal-level coloc, stage B -- QTL SuSiE + `coloc.susie`, one (analysis, tissue) task.
##
## For every (locus, gene) target this re-fits the gene's cis QTL with `susie_rss` on the
## GTEx v11 all-pairs nominal statistics, then colocalizes SIGNAL AGAINST SIGNAL against
## the cached GWAS SuSiE from stage A. Both modalities (sQTL, eQTL) are run on the same
## cell, so the modality contrast stays a within-gene comparison.
##
## WHY THE QTL SIDE IS RE-FIT
##   GTEx v11 ships `*.SuSiE_summary.parquet` -- a PIP per credible-set variant. It does
##   not ship `lbf_variable`, which is what `coloc.susie` colocalizes on. So the fit has
##   to be redone here from the nominal betas.
##
## THE LD PROXY, AND THE THREE THINGS DONE ABOUT IT
##   The only LD available for GTEx is the 1000G EUR panel already built per locus for
##   the CLPP layer, and GTEx brain donors are not a pure EUR sample. Under LD mismatch
##   `susie_rss` will happily return credible sets that are artefacts. So:
##     1. GWAS and QTL are fine-mapped on the SAME panel and the same SNP universe, so
##        the mismatch is common to the two traits rather than differential.
##     2. Every re-fit QTL credible set is checked against GTEx's OWN shipped credible
##        set for that phenotype (`cs_matches_gtex`). GTEx fit theirs with in-sample
##        genotypes; agreement is the external check. The Python meta stage makes the
##        agreeing set primary and reports the rest as a sensitivity arm.
##     3. `estimate_s_rss` is recorded per LOCUS by stage A and rides on every row as
##        `s_rss_gwas`, so a locus whose z-scores are inconsistent with the panel is
##        identifiable rather than silently trusted.
##
## ALLELE ORIENTATION MATTERS HERE, UNLIKE IN THE abf LAYER
##   `coloc.abf` depends on z^2 and is invariant to allele orientation, which is why
##   `20.coloc_modality_abf.R` harmonizes only to report a direction. `susie_rss` is NOT:
##   fine-mapping reads the z-vector jointly against the LD matrix, so the RELATIVE signs
##   across SNPs must match the panel's allele coding or the credible sets are wrong. The
##   QTL z is therefore re-signed to the panel REF allele, exactly as stage A does for
##   the GWAS z, and variants that cannot be harmonized are dropped from both arms alike.
##
## PRIOR SENSITIVITY IS FREE HERE
##   The SuSiE fits are the expensive part; `coloc.bf_bf` on the resulting log Bayes
##   factors is not. So the p12 sweep is run for EVERY cell rather than only for the
##   strong ones, and a hit can be reported as "colocalized for p12 in [a, b]" instead of
##   as a bare posterior at one arbitrary prior.
##
## Usage:  Rscript 05_genetic_anchoring/_h/23.coloc_signal_susie.R <analysis> <tissue> [chr]
## Env:    COLOC_SIGNAL_ARM        gene pool (default "switch"), mirrors the abf arms
##         COLOC_SIGNAL_SQTL       "representative" (default) | "all"
##         COLOC_GWAS_MAX_SNPS     GWAS SNP guard of the stage-A cache to use (default
##                                 12000 = primary; any other value is a sensitivity arm)
## Output: <signal dir>[/sensitivity/max_snps_<N>][/arms/<arm>][/all_introns]/susie/
##         <analysis>__<tissue>[.chr<N>].parquet
suppressPackageStartupMessages({
    library(data.table); library(arrow); library(susieR); library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: 23.coloc_signal_susie.R <analysis> <tissue> [chr]")
ANALYSIS <- args[1]
TISSUE   <- args[2]
ONE_CHR  <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
ROOT   <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
GTEX   <- "/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11"
PANEL_DIR <- "/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink"
ARM    <- Sys.getenv("COLOC_SIGNAL_ARM", unset = "switch")
SQTL_MODE <- Sys.getenv("COLOC_SIGNAL_SQTL", unset = "representative")
stopifnot(SQTL_MODE %in% c("representative", "all"))
## Same env var, same rule as stage A: a non-default GWAS SNP guard is a SENSITIVITY arm,
## read from and written to its own root, so it cannot overwrite the primary shards. Prep
## artefacts (targets, representative map, GTEx credible sets) do not depend on the guard
## and are always read from the primary arm dir.
MAX_SNPS_PRIMARY <- 12000L
MAX_SNPS <- suppressWarnings(as.integer(Sys.getenv("COLOC_GWAS_MAX_SNPS", "12000")))
if (is.na(MAX_SNPS)) stop("COLOC_GWAS_MAX_SNPS not an integer")

BASE <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc_signal_susie")
MDIR <- if (ARM == "switch") BASE else file.path(BASE, "arms", ARM)
if (!dir.exists(MDIR)) stop("missing arm dir ", MDIR, " (run --stage prep --arm ", ARM, ")")
SIG  <- if (MAX_SNPS == MAX_SNPS_PRIMARY) BASE else
        file.path(BASE, "sensitivity", sprintf("max_snps_%d", MAX_SNPS))
SDIR <- if (ARM == "switch") SIG else file.path(SIG, "arms", ARM)
GDIR   <- file.path(SIG, "gwas_susie", ANALYSIS)
if (!dir.exists(GDIR)) stop("missing GWAS SuSiE cache ", GDIR, " (run 22.coloc_gwas_susie.sh)")
CDIR   <- file.path(ROOT, "05_genetic_anchoring", "_m", "coloc", ANALYSIS, "susie")
BRIDGE <- file.path(ROOT, "inputs", "raw", "gtex_v11", "variant_bridge")
## The sqtl phenotype choice gets its own directory, exactly as the gene-pool arm does.
## Both modes write <analysis>__<tissue>.parquet, so without this split an all-introns
## run would silently overwrite the representative-intron shards it is meant to be
## compared against.
RDIR   <- if (SQTL_MODE == "representative") SDIR else file.path(SDIR, "all_introns")
OUTD   <- file.path(RDIR, "susie"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

L_MAX      <- 10L
MIN_SHARED <- 50L        # emit the row; the meta stage applies the real (100) floor
P12_SWEEP  <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
P12_PRIMARY <- 1e-5

## GTEx donors in this tissue -- both modalities are called in exactly these samples.
COVF <- file.path(GTEX, "GTEx_Analysis_v11_eQTL_covariates",
                  paste0(TISSUE, ".v11.covariates.txt"))
if (!file.exists(COVF)) stop("missing GTEx covariates file: ", COVF)
N_QTL <- length(strsplit(readLines(COVF, n = 1L), "\t")[[1]]) - 1L
message(sprintf("[%s / %s] arm=%s sqtl=%s GTEx donors=%d",
                ANALYSIS, TISSUE, ARM, SQTL_MODE, N_QTL))

targets <- as.data.table(read_parquet(file.path(MDIR, "targets.parquet")))
targets <- targets[analysis == ANALYSIS]
if (!nrow(targets)) { message("no targets for this analysis; nothing to do"); quit(status = 0) }

srep <- as.data.table(read_parquet(file.path(MDIR, "sqtl_representative.parquet")))
srep <- srep[tissue == TISSUE, .(gene, phenotype_id)]

## GTEx's own credible sets, for the agreement filter. Kept as variant_id; bridged to
## rsID per chromosome alongside the QTL statistics.
gcs <- as.data.table(read_parquet(file.path(MDIR, "gtex_credible_sets.parquet")))
gcs <- gcs[tissue == TISSUE]

is_ambig <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

read_ld <- function(lid) {
    vf <- file.path(CDIR, paste0(lid, ".unphased.vcor1.bin.vars"))
    bf <- file.path(CDIR, paste0(lid, ".unphased.vcor1.bin"))
    if (!file.exists(vf) || !file.exists(bf)) return(NULL)
    vars <- readLines(vf); n <- length(vars)
    con <- file(bf, "rb"); m <- readBin(con, "numeric", n = n*n, size = 4); close(con)
    R <- matrix(m, n, n, byrow = TRUE); dimnames(R) <- list(vars, vars); R
}

## One (gwas signal x qtl signal) posterior table for one cell, swept over p12.
## Returns NULL when either side has no credible set -- those cells fall back to
## coloc.abf in the Python meta stage and are flagged there, never dropped.
coloc_cell <- function(fit_g, fit_q, meta) {
    if (!length(fit_g$sets$cs) || !length(fit_q$sets$cs)) return(NULL)
    out <- vector("list", length(P12_SWEEP))
    for (k in seq_along(P12_SWEEP)) {
        r <- tryCatch(
            suppressWarnings(coloc::coloc.susie(fit_g, fit_q, p12 = P12_SWEEP[k])),
            error = function(e) NULL)
        if (is.null(r) || is.null(r$summary) || !nrow(r$summary)) next
        s <- as.data.table(r$summary)
        if (!"PP.H4.abf" %in% names(s)) next
        out[[k]] <- data.table(
            p12 = P12_SWEEP[k],
            idx1 = s$idx1, idx2 = s$idx2,
            nsnps = s$nsnps, hit1 = s$hit1, hit2 = s$hit2,
            PP0 = s$PP.H0.abf, PP1 = s$PP.H1.abf, PP2 = s$PP.H2.abf,
            PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf)
    }
    out <- rbindlist(out, fill = TRUE)
    if (!nrow(out)) return(NULL)
    cbind(as.data.table(meta), out)
}

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

    ## eQTL: resolve the versioned gene_id by bare match, then one pushdown read.
    ed  <- open_dataset(ef)
    gid <- as.data.table(ed |> dplyr::distinct(gene_id) |> dplyr::collect())
    gid[, bare := sub("\\..*$", "", gene_id)]
    want_e <- gid[bare %in% genes, gene_id]
    eq <- if (length(want_e)) as.data.table(
            ed |> dplyr::filter(gene_id %in% want_e) |>
            dplyr::select(gene_id, variant_id, pval_nominal, slope, slope_se) |>
            dplyr::collect()) else data.table()
    if (nrow(eq)) {
        eq[, gene := sub("\\..*$", "", gene_id)]
        eq[, phenotype_id := gene_id][, gene_id := NULL]
    }

    ## sQTL: GTEx's own representative intron per gene, or every intron of the gene.
    sd_ <- open_dataset(sf)
    if (SQTL_MODE == "representative") {
        want_s <- srep[gene %in% genes]
        sq <- if (nrow(want_s)) as.data.table(
                sd_ |> dplyr::filter(phenotype_id %in% want_s$phenotype_id) |>
                dplyr::select(phenotype_id, variant_id, pval_nominal, slope, slope_se) |>
                dplyr::collect()) else data.table()
        if (nrow(sq)) sq <- merge(sq, want_s, by = "phenotype_id")
    } else {
        ## phenotype_id is `chr:start:end:clu_N_strand:ENSG...`, so the gene is its last
        ## colon field. Every intron of the gene is tested; the meta stage must then
        ## account for the maximum over introns, which the representative arm avoids.
        pid <- as.data.table(sd_ |> dplyr::distinct(phenotype_id) |> dplyr::collect())
        pid[, gene := sub("\\..*$", "", sub("^.*:", "", phenotype_id))]
        want_s <- pid[gene %in% genes]
        sq <- if (nrow(want_s)) as.data.table(
                sd_ |> dplyr::filter(phenotype_id %in% want_s$phenotype_id) |>
                dplyr::select(phenotype_id, variant_id, pval_nominal, slope, slope_se) |>
                dplyr::collect()) else data.table()
        if (nrow(sq)) sq <- merge(sq, want_s, by = "phenotype_id")
    }

    bf <- file.path(BRIDGE, sprintf("chr%d.parquet", CH))
    if (!file.exists(bf)) stop("missing variant bridge: ", bf)
    br <- as.data.table(read_parquet(bf))

    ## panel REF/ALT for this chromosome: the allele coding the LD matrix is built on.
    bim <- fread(file.path(PANEL_DIR, sprintf("1000G.EUR.QC.%d.bim", CH)),
                 header = FALSE, col.names = c("bchr","rsid","cm","bpos","A1","A2"))
    bim <- bim[, .(rsid, REF = A2, ALT = A1)]

    prep_qtl <- function(q0) {
        if (!nrow(q0)) return(q0)
        q0 <- merge(q0, br, by = "variant_id")
        parts <- tstrsplit(q0$variant_id, "_", fixed = TRUE)
        q0[, `:=`(g_ref = parts[[3]], g_alt = parts[[4]])]
        q0 <- merge(q0, bim, by = "rsid")
        ## Harmonize the GTEx b38 alleles onto the hg19 panel coding, dropping
        ## strand-ambiguous and non-matching variants exactly as the GWAS side does.
        q0 <- q0[((g_ref == REF & g_alt == ALT) | (g_ref == ALT & g_alt == REF))
                 & !is_ambig(g_ref, g_alt)]
        ## slope is per the GTEx ALT allele; z_ref is per the panel REF allele, which is
        ## the orientation stage A put the GWAS z in.
        q0[, z_ref := (slope / slope_se) * fifelse(g_alt == REF, 1, -1)]
        q0[is.finite(z_ref)]
    }
    eq <- prep_qtl(eq); sq <- prep_qtl(sq)

    ## GTEx's own credible-set variants for this chromosome, as rsIDs.
    gcs_ch <- merge(gcs, br, by = "variant_id")
    gtex_cs_rsid <- split(gcs_ch$rsid, paste(gcs_ch$phenotype_id, gcs_ch$cs_id, sep = "|"))
    gtex_by_pheno <- split(gcs_ch$rsid, gcs_ch$phenotype_id)

    for (lid in unique(tg$LOCUS_ID)) {
        rds <- file.path(GDIR, paste0(lid, ".rds"))
        if (!file.exists(rds)) next
        G <- readRDS(rds)
        if (!length(G$fit$sets$cs)) { rm(G); next }
        R <- read_ld(lid)
        if (is.null(R)) { rm(G); next }
        Rg <- R[G$snps, G$snps, drop = FALSE]
        rm(R); gc(verbose = FALSE)

        lg <- tg[LOCUS_ID == lid]
        for (gi in seq_len(nrow(lg))) {
            tr <- lg[gi]
            for (MOD in c("eQTL", "sQTL")) {
                q0 <- if (MOD == "eQTL") eq else sq
                if (!nrow(q0)) next
                qg <- q0[gene == tr$gene]
                if (!nrow(qg)) next
                for (PH in unique(qg$phenotype_id)) {
                    qp <- qg[phenotype_id == PH][!duplicated(rsid)]
                    shared <- intersect(qp$rsid, G$snps)
                    if (length(shared) < MIN_SHARED) next
                    qp <- qp[match(shared, rsid)]
                    Rq <- Rg[shared, shared, drop = FALSE]

                    ## NB no per-cell `estimate_s_rss`: it is an O(p^3) eigen
                    ## decomposition and there are ~46k cells, so it would dominate the
                    ## entire run. The QTL-side LD check is `cs_matches_gtex` below --
                    ## agreement with GTEx's own in-sample fine-mapping, which is a
                    ## stronger check than a reference-LD consistency statistic anyway.
                    ## The locus-level `s_rss_gwas` from stage A rides on every row.
                    fq <- tryCatch(
                        susie_rss(z = qp$z_ref, R = Rq, n = N_QTL, L = L_MAX,
                                  estimate_residual_variance = FALSE),
                        error = function(e) NULL)
                    if (is.null(fq)) next
                    fq <- coloc::annotate_susie(fq, shared, Rq)

                    ## Does each re-fit credible set agree with GTEx's own?
                    ref_rsid <- gtex_by_pheno[[PH]]
                    cs_match <- vapply(fq$sets$cs, function(ix)
                        length(intersect(names(ix), ref_rsid)) > 0,
                        logical(1))
                    match_by_idx <- setNames(cs_match, as.character(fq$sets$cs_index))

                    meta <- list(
                        analysis = ANALYSIS, trait = tr$trait, LOCUS_ID = lid,
                        chr = CH, gene = tr$gene,
                        symbol = if ("symbol" %in% names(tr)) tr$symbol else NA_character_,
                        module_id = if ("module_id" %in% names(tr)) tr$module_id else NA_character_,
                        go_invisible = if ("go_invisible" %in% names(tr)) tr$go_invisible else NA,
                        tissue = TISSUE, modality = MOD, phenotype_id = PH,
                        n_shared = length(shared), n_qtl_variants = nrow(qp),
                        n_cs_gwas = length(G$fit$sets$cs), n_cs_qtl = length(fq$sets$cs),
                        n_qtl_donors = N_QTL,
                        s_rss_gwas = G$kriging$s,
                        qtl_min_p = min(qp$pval_nominal, na.rm = TRUE))
                    r <- coloc_cell(G$fit, fq, meta)
                    if (is.null(r)) next
                    r[, cs_matches_gtex := unname(match_by_idx[as.character(idx2)])]
                    r[is.na(cs_matches_gtex), cs_matches_gtex := FALSE]
                    results[[length(results) + 1]] <- r
                    rm(fq, Rq); gc(verbose = FALSE)
                }
            }
        }
        rm(Rg, G); gc(verbose = FALSE)
        message(sprintf("  %-18s done (%d result blocks so far)", lid, length(results)))
    }
    rm(eq, sq, br, bim); gc(verbose = FALSE)
}

if (!length(results)) {
    message("no coloc.susie rows produced for ", ANALYSIS, " / ", TISSUE)
    quit(status = 0)
}
res <- rbindlist(results, fill = TRUE)
suffix <- if (is.na(ONE_CHR)) "" else sprintf(".chr%d", ONE_CHR)
outf <- file.path(OUTD, sprintf("%s__%s%s.parquet", ANALYSIS, TISSUE, suffix))
write_parquet(res, outf)
message(sprintf("\n[%s / %s] %d signal-pair rows over %d genes -> %s",
                ANALYSIS, TISSUE, nrow(res), uniqueN(res$gene), outf))
prim <- res[p12 == P12_PRIMARY & cs_matches_gtex == TRUE]
if (nrow(prim))
    message(sprintf("  primary (p12=%g, GTEx-matched): %d rows, PP4>=0.8 in %d",
                    P12_PRIMARY, nrow(prim), sum(prim$PP4 >= 0.8)))
cat("\nReproducibility\n"); print(Sys.time()); print(sessionInfo())
