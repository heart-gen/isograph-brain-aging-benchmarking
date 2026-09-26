#!/usr/bin/env bash
## Stage 03 run order, as a SLURM dependency graph. The tier is the leading number of each
## wrapper: steps in one tier run in parallel once the steps they wait on have finished.
##
##   bash 03_module_characterization/_h/run_stage.sh --dry-run
##   bash 03_module_characterization/_h/run_stage.sh --after <stage 02 job ids>
## Options: scripts/slurm_dag.sh. Inputs: every stage-02 fit.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/scripts/slurm_dag.sh"
dag_init 03_module_characterization "$@"
H=03_module_characterization/_h

## 01 -- per-fit characterization. --force: interpret_modules keeps existing outputs otherwise.
step 01a "" $H/01a.interpret_modules_brainseq.sh --force
step 01b "" $H/01b.interpret_modules_gtex.sh --force
step 01c "" $H/01c.module_enrichment_brainseq.sh
step 01d "" $H/01d.module_enrichment_gtex.sh
step 01e "" $H/01e.incremental_association_brainseq.sh
step 01f "" $H/01f.incremental_association_gtex.sh
step 01g "" $H/01g.celltype_composition_brainseq.sh
step 01h "" $H/01h.celltype_composition_gtex.sh

## 02 -- joins over tier 01
step 02a "01a 01c" $H/02a.go_invisible_gate.sh
step 02b "01c 01e" $H/02b.characterize_composition_unique.sh
step 02c "01g 01h" $H/02c.composition_meta.sh

## 03
step 03a "01e 01f 02b" $H/03a.abundance_structure.sh
step 03b "03a" $H/03b.abundance_structure_all.sh   # orthogonality for all 17 stores; pool with _h/03c (local)

dag_finish
