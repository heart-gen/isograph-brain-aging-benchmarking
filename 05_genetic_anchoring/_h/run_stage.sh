#!/usr/bin/env bash
## Stage 05 run order, as a SLURM dependency graph. The tier is the leading number of each
## wrapper: steps in one tier run in parallel once the steps they wait on have finished. A few
## same-tier steps wait on a sibling because they write a shared file (noted inline).
##
##   bash 05_genetic_anchoring/_h/run_stage.sh --dry-run
##   bash 05_genetic_anchoring/_h/run_stage.sh --after <stage 04 job ids>
## Options: scripts/slurm_dag.sh. Inputs: stage-02 fits (incl. isograph_vae_res5 and the matched
## WGCNA baselines), stage-03 module_enrichment and structure_switch_pairs.
##
## Re-running over an existing tree: several steps SKIP outputs that already exist or glob-collect
## whatever is on disk (09 per-locus LD matrices, the 03c GWAS SuSiE cache, 02f LD-score annotations,
## the SMR task tree). Clear those first -- see "Re-running" in 05_genetic_anchoring/README.md.
##
## Not submitted here: 08d.locus_ld_robustness (per locus, arguments chosen from the nominations),
## and 10a/11a, which 09a submits once prep has sized the SMR arrays.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/scripts/slurm_dag.sh"
dag_init 05_genetic_anchoring "$@"
H=05_genetic_anchoring/_h
COLOC=(aging__ad aging__als aging__lbd aging__pd aging__scz brainseq-sczd__scz)
ARMS=(switch background wgcna_switch wgcna_multiplex)
RES5=--export=ALL,MAGMA_ISOGRAPH_BACKEND=isograph_vae_res5

## 01 -- everything that reads only stages 02-04
step 01a "" $H/01a.qtl_anchoring.sh
step 01b "" $H/01b.qtl_anchoring_matched.sh
step 01b.constraint "" $H/01b.qtl_anchoring_matched.sh --outcome binary --covariate-set constraint
step 01b.continuous "" $H/01b.qtl_anchoring_matched.sh --outcome continuous
step 01c "" $H/01c.qtl_anchoring_sensitivity.sh
step 01d "" $H/01d.sqtl_concordance.sh
step 01e "" $H/01e.module_anchoring.sh
step 01f "" $H/01f.prep_module_gene_sets.sh
step 01f.res5 "" "${RES5}" $H/01f.prep_module_gene_sets.sh
for t in ad als lbd pd scz; do
    step "01g.aging__${t}" "" --time=03:00:00 $H/01g.coloc_prep.sh --gene-source aging --trait ${t} --min-recurrence 1
done
# shares coloc/_tmp/_scz_rsids.txt with aging__scz
step 01g.brainseq-sczd__scz "01g.aging__scz" --time=03:00:00 $H/01g.coloc_prep.sh --gene-source brainseq-sczd --trait scz
step 01h.brainseq-sczd "" $H/01h.ldsc_annot_prep.sh --analysis brainseq-sczd
# shares ldsc/_tmp/_want_variants.txt with the brainseq-sczd prep
step 01h.aging "01h.brainseq-sczd" $H/01h.ldsc_annot_prep.sh --bundle aging --min-recurrence 1
step 01i.all_samples "" $H/01i.brainseq_switch_qtl.sh
step 01i.ea_only "" --export=ALL,SWQTL_ARM=ea_only $H/01i.brainseq_switch_qtl.sh

## 02
step 02a "01a 01b" $H/02a.qtl_anchoring_meta.sh
step 02a.constraint "01c 01b.constraint" $H/02a.qtl_anchoring_meta.sh --outcome binary --covariate-set constraint
step 02a.continuous "01c 01b.continuous" $H/02a.qtl_anchoring_meta.sh --outcome continuous
step 02a.dose "01c" $H/02a.qtl_anchoring_meta.sh --outcome dose
step 02b "01d" $H/02b.sqtl_concordance_meta.sh
step 02c "01e" $H/02c.module_anchoring_meta.sh
step 02d "01f" $H/02d.run_magma.sh
# waits on the canonical run so the shared SNP p-value and gene-analysis caches exist
step 02d.res5 "01f.res5 02d" "${RES5}" $H/02d.run_magma.sh
# memory tracks the largest locus: a square float32 LD matrix is 4n^2 bytes, and the aging
# recurrence-1 gene set reaches 54k SNPs (AD, 11.7 GB) and 82k (SCZ, 27 GB)
for a in "${COLOC[@]}"; do
    cpus=16; [[ "${a}" == aging__scz ]] && cpus=32
    step "02e.${a}" "01g.${a}" --time=06:00:00 --cpus-per-task=${cpus} $H/02e.locus_ld.sh "${a}"
done
step 02f.brainseq-sczd "01h.brainseq-sczd" $H/02f.ldsc_make_annot_ldscores.sh brainseq-sczd
step 02f.aging "01h.aging" $H/02f.ldsc_make_annot_ldscores.sh aging
step 02g.all_samples "01i.all_samples" $H/02g.brainseq_qtl_checks.sh
step 02g.ea_only "01i.ea_only" --export=ALL,SWQTL_ARM=ea_only $H/02g.brainseq_qtl_checks.sh
step 02h "01i" $H/02h.brainseq_switch_qtl_meta.sh

## 03
step 03a "02d" $H/03a.plot_magma.sh
step 03a.res5 "02d.res5 02d" "${RES5}" $H/03a.plot_magma.sh
for a in "${COLOC[@]}"; do
    step "03b.${a}" "02e.${a}" --time=08:00:00 $H/03b.coloc_clpp.sh "${a}"
done
step 03c.primary "02e" $H/03c.coloc_gwas_susie.sh
# scoped SNP-guard sensitivity arm (2026-09-10): aging__ad only, at 30,000 SNPs
step 03c.max_snps_30000 "02e.aging__ad" --array=0 --cpus-per-task=32 \
    --export=ALL,COLOC_GWAS_MAX_SNPS=30000 $H/03c.coloc_gwas_susie.sh
step 03d.brainseq-sczd__scz "02f.brainseq-sczd" $H/03d.ldsc_munge_h2.sh scz brainseq-sczd
for t in ad als lbd pd scz; do
    step "03d.aging__${t}" "02f.aging" $H/03d.ldsc_munge_h2.sh ${t} aging
done
step 03e "02h" $H/03e.brainseq_effect_size.sh

## 04
step 04a "03b" $H/04a.coloc_meta.sh
step 04b "03b" $H/04b.coloc_direction.sh
step 04c "03b 03a" $H/04c.module_coloc_convergence.sh
step 04d.switch "03b" --export=ALL,COLOC_MODALITY_ARM=switch $H/04d.coloc_modality_prep.sh
for arm in background wgcna_switch wgcna_multiplex; do
    # after the switch arm, which builds the shared variant bridge
    step "04d.${arm}" "04d.switch" --export=ALL,COLOC_MODALITY_ARM=${arm} $H/04d.coloc_modality_prep.sh
done
step 04e "03d" $H/04e.ldsc_summary.sh

## 05
for arm in "${ARMS[@]}"; do
    step "05a.${arm}" "04d.${arm}" --export=ALL,COLOC_MODALITY_ARM=${arm} $H/05a.coloc_modality_abf.sh
done
step 05b "04d.switch" $H/05b.coloc_signal_susie_prep.sh
step 05c "04b" $H/05c.deep_dive_events.sh

## 06
for arm in "${ARMS[@]}"; do
    step "06a.${arm}" "05a.${arm}" $H/06a.coloc_modality_meta.sh --arm ${arm}
done
step 06b.representative "03c.primary 05b" $H/06b.coloc_signal_susie.sh
step 06b.all_introns "03c.primary 05b" --time=24:00:00 --cpus-per-task=24 \
    --export=ALL,COLOC_SIGNAL_SQTL=all $H/06b.coloc_signal_susie.sh
step 06b.max_snps_30000_representative "03c.max_snps_30000 05b" --array=0-12 --cpus-per-task=32 \
    --export=ALL,COLOC_GWAS_MAX_SNPS=30000 $H/06b.coloc_signal_susie.sh
step 06b.max_snps_30000_all_introns "03c.max_snps_30000 05b" --array=0-12 --cpus-per-task=32 --time=24:00:00 \
    --export=ALL,COLOC_GWAS_MAX_SNPS=30000,COLOC_SIGNAL_SQTL=all $H/06b.coloc_signal_susie.sh
step 06c "03c.primary 05b 02g.ea_only" $H/06c.coloc_brainseq_prep.sh

## 07
step 07a "06a" $H/07a.coloc_modality_compare.sh
step 07b.representative "06b.representative" $H/07b.coloc_signal_susie_meta.sh
step 07b.all_introns "06b.all_introns" $H/07b.coloc_signal_susie_meta.sh --sqtl all
step 07b.max_snps_30000_representative "06b.max_snps_30000_representative" \
    $H/07b.coloc_signal_susie_meta.sh --max-snps 30000
step 07b.max_snps_30000_all_introns "06b.max_snps_30000_all_introns" \
    $H/07b.coloc_signal_susie_meta.sh --sqtl all --max-snps 30000
step 07c "06c" $H/07c.coloc_brainseq_susie.sh

## 08 -- the locus audits run one after another: they share the curated-event testability cache
step 08a "07b.all_introns 04b" $H/08a.coloc_isoform_events_signal.sh
step 08b.abf "06a.switch 05b 04b" $H/08b.locus_event_audit.sh abf
step 08b.susie "08b.abf 07b.representative" $H/08b.locus_event_audit.sh susie
step 08b.susie_all_introns "08b.susie 07b.all_introns" $H/08b.locus_event_audit.sh susie --sqtl-arm all
step 08b.susie_all_introns_max_snps_30000 "08b.susie_all_introns 07b.max_snps_30000_all_introns" \
    $H/08b.locus_event_audit.sh susie --sqtl-arm all --max-snps 30000
step 08c "07c 07b.all_introns" $H/08c.coloc_brainseq_meta.sh
manual_step 08d "08d.locus_ld_robustness.sh <analysis> <LOCUS_ID> <gene> <tissue[,...]>, for each locus that will carry a biological claim"

## 09 -- SMR + HEIDI: each submitter runs prep, then submits its own 10a arrays and 11a metas
step 09a.gtex "07b.all_introns 03b" --export=ALL,SMR_SUBMIT_LEDGER=${DAG_LEDGER} $H/09a.smr_heidi_submit.sh
step 09a.brainseq "08c 07b.all_introns 02g.ea_only" \
    --export=ALL,SMR_QTL_SOURCE=brainseq,SMR_ARM=ea_only,SMR_SUBMIT_LEDGER=${DAG_LEDGER} $H/09a.smr_heidi_submit.sh
manual_step 10a "10a.smr_heidi.sh + 11a.smr_heidi_meta.sh are submitted by 09a (array size = work_list.tsv)"

dag_finish
