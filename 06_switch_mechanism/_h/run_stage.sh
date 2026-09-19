#!/usr/bin/env bash
## Stage 06 run order, as a SLURM dependency graph. The tier is the leading number of each
## wrapper: steps in one tier run in parallel once the steps they wait on have finished.
##
##   bash 06_switch_mechanism/_h/run_stage.sh --dry-run
##   bash 06_switch_mechanism/_h/run_stage.sh --after <stage 05 job ids>
## Options: scripts/slurm_dag.sh. Inputs: stage-03 structure_switch_pairs, stage-05 coloc events
## (04b), deep_dive_events (05c) and signal-layer events (08a); GTEx junction usage from
## inputs/_h/build_gtex_junction_usage.sh.
##
## 01f downloads ClinVar/gnomAD. It runs as an ordinary batch step -- compute nodes do reach
## the network -- and skips files already present, so re-running the stage is cheap.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/scripts/slurm_dag.sh"
dag_init 06_switch_mechanism "$@"
H=06_switch_mechanism/_h

## 01
step 01a "" $H/01a.switch_consequence.sh
step 01b "" $H/01b.validate_switch_splicing_brainseq.sh
step 01c "" $H/01c.validate_switch_splicing_gtex.sh
step 01d "" $H/01d.isa_concordance.sh
step 01e "01a" $H/01e.longread_switch_confirm.sh   # reads 01a's pair_consequence
step 01f "" $H/01f.download_clinical.sh
step 01g "" $H/01g.scz_confound_sensitivity.sh
step 01h "" $H/01h.switch_feature_sensitivity.sh --cohort brainseq --region caudate
step 01i "" $H/01i.switch_feature_sensitivity_gtex.sh
step 01j "" $H/01j.switch_feature_refit.sh
step 01k "" $H/01k.junction_coloc_confirm.sh

## 02 -- the three orthogonal-confirmation modes share an output directory, so they run in turn
step 02a "01a" $H/02a.switch_consequence_meta.sh
step 02b.anchored "01e" $H/02b.switch_orthogonal_confirm.sh --mode anchored
step 02b.global_null "02b.anchored" $H/02b.switch_orthogonal_confirm.sh --mode global-null
step 02b.signal "02b.global_null" $H/02b.switch_orthogonal_confirm.sh --events signal
step 02c "01a 01f" $H/02c.clinical_consequence.sh
step 02d "01h 01i" $H/02d.switch_feature_sensitivity_aggregate.sh
step 02e "01j" $H/02e.switch_feature_refit_aggregate.sh

## 03
step 03a "02c" $H/03a.clinical_consequence_meta.sh

## Quest only (PI item 10a): the WASP BAMs exist only on Quest, so these four are never
## submitted from Bridges-2. On Quest: bash 06_switch_mechanism/_h/ase_junction_quest.sh
manual_step 01l "Quest only: 01l.ase_junction_targets.sh (see the stage README)"
manual_step 02f "Quest only: 02f.ase_junction_count.sh, one task per WASP BAM, after 01l"
manual_step 03b "Quest only: 03b.ase_junction_screen.sh, after 02f"
manual_step 04a "Quest only: 04a.ase_junction_allelic.sh, after 03b (only if its gate passed)"

dag_finish
