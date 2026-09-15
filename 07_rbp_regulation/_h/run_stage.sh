#!/usr/bin/env bash
## Stage 07 run order, as a SLURM dependency graph. The tier is the leading number of each
## wrapper: steps in one tier run in parallel once the steps they wait on have finished.
##
##   bash 07_rbp_regulation/_h/run_stage.sh --dry-run
##   bash 07_rbp_regulation/_h/run_stage.sh --after <stage 06 job ids>
## Options: scripts/slurm_dag.sh. Inputs: stage-03 structure_switch_pairs, GENCODE v47, ATtRACT,
## and the neuronal CLIP downloads from inputs/_h/download_neuronal_clip.sh.
##
## 03b fetches ENCODE eCLIP peaks for the RBPs 02a nominates and needs a login node. The first
## submission therefore holds 04a; once 02a has finished, run 03b on a login node and resubmit
## the held step:
##   bash 07_rbp_regulation/_h/03b.rbp_binding_fetch.sh
##   bash 07_rbp_regulation/_h/run_stage.sh --only 04a --login-done 03b
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/scripts/slurm_dag.sh"
dag_init 07_rbp_regulation "$@"
H=07_rbp_regulation/_h

step 01a "" $H/01a.rbp_motif_families.sh
step 02a "01a" $H/02a.rbp_regulon.sh
step 03a "02a" $H/03a.rbp_regulon_intronic.sh
login_step 03b "02a" $H/03b.rbp_binding_fetch.sh
step 04a "01a 02a 03b" $H/04a.rbp_binding.sh
step 04b "03a" $H/04b.neuronal_clip_freeze.sh
step 05a "04b" $H/05a.neuronal_clip_windows.sh
step 06a "05a" $H/06a.neuronal_clip_motif_qc.sh
step 07a "06a" $H/07a.neuronal_clip_overlap.sh
step 07b "06a" $H/07b.nova_family_renomination.sh
step 08a "07b" $H/08a.nova2_ctag_clip.sh
step 09a "08a" $H/09a.nova2_perturbation.sh

dag_finish
