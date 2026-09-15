#!/usr/bin/env bash
## Stage 08 run order, as a SLURM dependency graph. Every step here reads more than one of stages
## 04-07, which is why these analyses sit after all of them.
##
##   bash 08_integration/_h/run_stage.sh --dry-run
##   bash 08_integration/_h/run_stage.sh --after <stage 07 job ids>
## Options: scripts/slurm_dag.sh.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/scripts/slurm_dag.sh"
dag_init 08_integration "$@"
H=08_integration/_h

## 01
step 01a "" $H/01a.build_deep_dive.sh
step 01b "" $H/01b.scz_age_projection.sh
step 01c.isograph_linear "" --array=0 $H/01c.replication_functional.sh
step 01c.wgcna_linear    "" --array=1 $H/01c.replication_functional.sh
step 01c.isograph_spline "" --array=2 $H/01c.replication_functional.sh
step 01c.wgcna_spline    "" --array=3 $H/01c.replication_functional.sh

## 02 -- one wet-lab panel per RBP arm of WETLAB_PERTURBATION_DESIGN.md
for rbp in KHDRBS1 NONO ELAVL1; do
    step "02a.${rbp}" "01a" $H/02a.build_rbp_target_panel.sh --rbp ${rbp}
done

## 03
step 03a "02a" $H/03a.rbp_pair_assayability.sh

dag_finish
