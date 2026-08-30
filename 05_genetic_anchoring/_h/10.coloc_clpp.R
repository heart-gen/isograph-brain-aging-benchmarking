#!/usr/bin/env Rscript
## Coloc capstone, step 3 — GWAS fine-mapping (SuSiE) + eCAVIAR CLPP of IsoGraph
## switch genes vs GTEx brain sQTL AND eQTL. Adapts the organoid pipeline
## (colocalization _h/05.susie_clpp.R), extended to the paired sQTL/eQTL contrast
## that is the point of the IsoGraph capstone.
##
## Method (per locus from coloc_prep.py):
##   1. Read the plink2 REF-based LD matrix + per-locus PGC3 z, re-sign each z to
##      the panel REF allele (drop allele-mismatch / strand-ambiguous SNPs), and
##      run susie_rss for a GWAS causal PIP per SNP.
##   2. For each switch gene in the locus, each QTL kind (sQTL, eQTL), tissue and
##      credible set, map QTL variants (b38) -> rsID and compute
##        CLPP = sum_snp PIP_gwas * PIP_qtl   over the shared SNPs.
##      CLPP >= 0.01 = colocalized (eCAVIAR convention); also flag >= 0.05.
## coloc.susie is not run: GTEx v11 ships only credible-set summaries (PIP per
## variant), not full SuSiE objects, so eCAVIAR CLPP is the appropriate estimator.
##
## Why paired: sQTL CLPP > (or comparable to) eQTL CLPP for the same switch genes,
## concentrated in GO-invisible modules, is splicing-specific genetic anchoring of
## the switch direction. GO-invisible tags come from coloc_prep candidate_genes.tsv.
##
## Usage:  Rscript 05_genetic_anchoring/_h/10.coloc_clpp.R <analysis>   (default brainseq-sczd)
## Output (05_genetic_anchoring/_m/coloc/<analysis>/coloc/):
##   clpp_results.tsv       — gene x kind x tissue x QTL-CS, CLPP + top shared SNP
##   coloc_gene_summary.tsv — per gene x kind: best tissue/CS by CLPP
##   gwas_susie_cs.tsv      — GWAS credible sets per locus
suppressPackageStartupMessages({
    library(data.table)
    library(susieR)
})

args     <- commandArgs(trailingOnly = TRUE)
ANALYSIS <- if (length(args) >= 1) args[1] else "brainseq-sczd"
ROOT     <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
MDIR     <- file.path(ROOT, "real_data", "coloc", "_m", ANALYSIS)
SUSIE_D  <- file.path(MDIR, "susie")
OUTD     <- file.path(MDIR, "coloc")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
PANEL_DIR <- "/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink"
CLPP_MIN <- 0.01
L_MAX    <- 10
## Skip pathologically large loci: an N x N float64 LD matrix at N ~ 18k is ~2.6 GB and
## susie_rss holds several copies (OOM risk); such wide merged loci sit in long-range LD
## where fine-mapping is unreliable anyway. Reported as skipped.
MAX_SNPS <- 12000

loci  <- fread(file.path(SUSIE_D, "loci_testable.tsv"))
cs    <- fread(file.path(MDIR, "qtl_credible_sets.tsv"))
vmap  <- fread(file.path(MDIR, "variant_rsid_map.tsv"),
               header = FALSE, col.names = c("variant_id", "rsid"))
cs    <- merge(cs, vmap, by = "variant_id", all.x = TRUE)
cand  <- fread(file.path(MDIR, "candidate_genes.tsv"))
gotag <- unique(cand[, .(gene, go_invisible)])

## hg19 regions to drop from fine-mapping (long-range LD): MHC for every trait, plus
## APOE for AD/LBD. Written per-trait by coloc_prep.py; falls back to MHC if absent.
excl_f <- file.path(MDIR, "exclude_regions.tsv")
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

all_clpp <- list(); all_cs <- list()
for (i in seq_len(nrow(loci))) {
    L <- loci[i]; lid <- L$LOCUS_ID
    R <- read_ld(lid)
    if (is.null(R)) { message("  no LD for ", lid, "; skipping"); next }
    gwas <- fread(file.path(SUSIE_D, paste0(lid, ".gwas.tsv")))
    ## Drop SNPs in the trait's long-range-LD exclusion regions (MHC always; APOE for
    ## AD/LBD) at the SNP level, matching the MAGMA pipeline: long-range LD makes
    ## fine-mapping / coloc unreliable. Loci fully inside a region collapse below the
    ## SNP floor and are skipped; loci that only clip an edge keep their outside SNPs.
    for (j in seq_len(nrow(excl)))
        if (L$chr == excl$chr[j])
            gwas <- gwas[!(pos >= excl$start[j] & pos <= excl$stop[j])]
    bim  <- fread(file.path(PANEL_DIR, sprintf("1000G.EUR.QC.%d.bim", L$chr)),
                  header = FALSE, col.names = c("bchr","rsid","cm","bpos","A1","A2"))
    ## plink1 .bed convention: REF = A2 (col6), ALT = A1 (col5)
    bim <- bim[rsid %in% rownames(R), .(rsid, REF = A2, ALT = A1)]
    dt  <- merge(gwas, bim, by = "rsid")
    dt  <- dt[((a1==REF & a2==ALT) | (a1==ALT & a2==REF)) & !is_ambig(a1, a2)]
    dt[, z_ref := z * fifelse(a1==REF, 1, -1)]
    dt  <- dt[rsid %in% rownames(R)][!duplicated(rsid)]
    if (nrow(dt) < 20) { message("  ", lid, ": <20 usable SNPs; skipping"); next }
    if (nrow(dt) > MAX_SNPS) {
        message("  ", lid, ": ", nrow(dt), " SNPs > MAX_SNPS (", MAX_SNPS,
                "); long-range-LD locus, skipping")
        rm(R, gwas, bim, dt); gc(verbose = FALSE); next
    }
    Rk  <- R[dt$rsid, dt$rsid, drop = FALSE]
    nn  <- max(1000, round(L$min_neff))
    fit <- tryCatch(susie_rss(z = dt$z_ref, R = Rk, n = nn, L = L_MAX,
                              estimate_residual_variance = FALSE),
                    error = function(e) { message("  susie failed ", lid, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    dt[, pip_gwas := fit$pip]
    gp <- dt[, .(rsid, pos, pip_gwas, beta, a1_gwas = a1, a2_gwas = a2, p)]
    if (length(fit$sets$cs)) for (csn in names(fit$sets$cs)) {
        idx <- fit$sets$cs[[csn]]
        all_cs[[length(all_cs)+1]] <- data.table(
            LOCUS_ID = lid, gwas_cs = csn, n = length(idx),
            lead_rsid = dt$rsid[idx][which.max(dt$pip_gwas[idx])],
            lead_pip = max(dt$pip_gwas[idx]))
    }
    locus_genes <- unlist(strsplit(L$genes, ","))
    for (g in intersect(unique(cs$gene), locus_genes)) {
        eg <- cs[gene == g & !is.na(rsid)]
        if (!nrow(eg)) next
        for (kd in unique(eg$kind)) for (rg in unique(eg[kind==kd, tissue]))
            for (cid in unique(eg[kind==kd & tissue==rg, cs_id])) {
                e <- eg[kind==kd & tissue==rg & cs_id==cid, .(rsid, pip_qtl = pip)]
                m <- merge(e, gp, by = "rsid")
                if (!nrow(m)) {
                    all_clpp[[length(all_clpp)+1]] <- data.table(
                        gene=g, kind=kd, tissue=rg, qtl_cs=cid, LOCUS_ID=lid,
                        n_qtl=nrow(e), n_shared=0L, clpp=0, best_rsid=NA_character_,
                        best_pos=NA_integer_, pip_gwas=NA_real_, pip_qtl=NA_real_, gwas_p=NA_real_)
                    next
                }
                m[, prod := pip_gwas * pip_qtl]
                b <- m[which.max(prod)]
                all_clpp[[length(all_clpp)+1]] <- data.table(
                    gene=g, kind=kd, tissue=rg, qtl_cs=cid, LOCUS_ID=lid,
                    n_qtl=nrow(e), n_shared=nrow(m), clpp=sum(m$prod),
                    best_rsid=b$rsid, best_pos=b$pos,
                    pip_gwas=b$pip_gwas, pip_qtl=b$pip_qtl, gwas_p=b$p)
            }
    }
    message(sprintf("  %-16s susie_snps=%5d  max_pip=%.3f  gwas_cs=%d",
                    lid, nrow(dt), max(dt$pip_gwas), length(fit$sets$cs)))
    rm(R, Rk, fit, dt, bim, gwas); gc(verbose = FALSE)
}

clpp <- rbindlist(all_clpp)
if (!nrow(clpp)) { message("No CLPP rows produced."); quit(status = 0) }
clpp <- merge(clpp, gotag, by = "gene", all.x = TRUE)
clpp[, colocalized := clpp >= CLPP_MIN]
setorder(clpp, -clpp)
fwrite(clpp, file.path(OUTD, "clpp_results.tsv"), sep = "\t")

## per gene x kind: best tissue/CS by CLPP
summ <- clpp[order(-clpp), .SD[1], by = .(gene, kind)]
summ <- merge(summ, loci[, .(LOCUS_ID, true_lead_snp, true_lead_p)], by = "LOCUS_ID", all.x = TRUE)
setorder(summ, -clpp)
fwrite(summ, file.path(OUTD, "coloc_gene_summary.tsv"), sep = "\t")
if (length(all_cs)) fwrite(rbindlist(all_cs), file.path(OUTD, "gwas_susie_cs.tsv"), sep = "\t")

## quick console contrast: colocalized (CLPP>=0.01) gene counts by kind x GO tag
message("\n================ COLOCALIZED (CLPP >= ", CLPP_MIN, ") ================")
tab <- clpp[colocalized == TRUE, .(genes = uniqueN(gene)), by = .(kind, go_invisible)]
print(tab[order(kind, go_invisible)])
message("\nGenes tested: ", uniqueN(clpp$gene),
        " | colocalized (sQTL): ", uniqueN(clpp[kind=="sQTL" & colocalized==TRUE, gene]),
        " | colocalized (eQTL): ", uniqueN(clpp[kind=="eQTL" & colocalized==TRUE, gene]),
        " | tests: ", nrow(clpp))

cat("\nReproducibility\n"); Sys.time(); print(sessionInfo())
