#!/usr/bin/env bash
## Stage 04 run order, as a SLURM dependency graph. The tier is the leading number of each
## wrapper: steps in one tier run in parallel once the steps they wait on have finished.
##
##   bash 04_module_trust/_h/run_stage.sh --dry-run
##   bash 04_module_trust/_h/run_stage.sh --after <stage 03 job ids>
## Options: 00_scripts/slurm_dag.sh. Inputs: stage-02 fits; stage-03 interpretation, enrichment and
## composition-unique genes (read by 03c-03h).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/00_scripts/slurm_dag.sh"
dag_init 04_module_trust "$@"
H=04_module_trust/_h
TRUST=(brainseq:caudate brainseq:hippocampus brainseq:dlpfc
       gtex:caudate_basal_ganglia gtex:hippocampus gtex:frontal_cortex_ba9)

## 01 -- split-half fits (graphs saved for the resolution sweep), cross-cohort matching, curvature
step 01a "" --export=ALL,STABILITY_SAVE_EDGES=1 $H/01a.stability_isograph.sh
step 01b "" $H/01b.stability_wgcna.sh
step 01c "" $H/01c.replication.sh
step 01d "" $H/01d.age_model_curvature.sh

## 02
step 02a "01a" $H/02a.stability_resolution_sweep.sh
for cr in "${TRUST[@]}"; do
    c=${cr%%:*}; r=${cr##*:}
    for m in isograph wgcna; do
        step "02b.${c}_${r}_${m}" "01a 01b" --export=ALL,COHORT=${c},REGION=${r},METHOD=${m} $H/02b.module_meta.sh
    done
done
step 02c "01a 01b" $H/02c.trust_gate.sh
step 02d "01c"     $H/02d.replication_go.sh
step 02e "01a"     $H/02e.dtu_added_value_saturn.sh

## 03
step 03a "01a 01b 02a" $H/03a.stability_aggregate.sh
step 03b.isograph "02b" --export=ALL,METHOD=isograph $H/03b.within_cohort.sh
step 03b.wgcna    "02b" --export=ALL,METHOD=wgcna    $H/03b.within_cohort.sh
step 03c.linear "02c" $H/03c.module_trust_replication.sh
step 03c.spline "02c" $H/03c.module_trust_replication.sh --model spline
for cov in full complement none; do
    step "03d.${cov}" "02c" $H/03d.replication_permutation.sh --covariates ${cov}
done
step 03e "02c" $H/03e.replication_pooled.sh
step 03f "02c" $H/03f.eigengene_projection.sh
step 03g "02c" $H/03g.complementarity.sh
step 03h "02d" $H/03h.baseline_comparison.sh
step 03i "01a 02e" $H/03i.dtu_added_value_analyze.sh

## 04
step 04a "03c" $H/04a.replication_model_contrast.sh
step 04b "03d" $H/04b.replication_permutation_report.sh
step 04c "03f" $H/04c.eigengene_projection_aggregate.sh
step 04d "03i" $H/04d.dtu_added_value_summarize.sh

dag_finish
