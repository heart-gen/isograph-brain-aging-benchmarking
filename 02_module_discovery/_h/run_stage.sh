#!/usr/bin/env bash
## Stage 02 run order, as a SLURM dependency graph. The tier is the leading number of each
## wrapper: steps in one tier run in parallel once the steps they wait on have finished.
##
##   bash 02_module_discovery/_h/run_stage.sh --dry-run     # print the plan, submit nothing
##   bash 02_module_discovery/_h/run_stage.sh               # submit
##   bash 02_module_discovery/_h/run_stage.sh --only 02     # e.g. only the tier-02 steps
## Options (--after, --from/--to, --only/--skip, --login-done): scripts/slurm_dag.sh.
## Inputs: the IsoGraph bundles from inputs/_h/build_data_pipeline.sh.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/scripts/slurm_dag.sh"
dag_init 02_module_discovery "$@"
H=02_module_discovery/_h

## 00 -- the shared gene universe (production transcript filter applied to each bundle).
## The gene-level WGCNA baselines read it so they are fit on the genes IsoGraph models.
step 00a "" $H/00a.production_gene_universe.sh

## 01 -- fits and the three WGCNA baselines, straight from the bundles
step 01a ""  $H/01a.run_isograph_brainseq_aging.sh
step 01b ""  $H/01b.run_isograph_brainseq_sczd.sh
step 01c ""  $H/01c.run_isograph_gtex.sh
# resolution-2.0 refits (isograph_vae_res2): the disclosed resolution comparison behind S-real-2
step 01a.res2 "" $H/01a.run_isograph_brainseq_aging.sh --leiden-resolution 2.0
step 01b.res2 "" $H/01b.run_isograph_brainseq_sczd.sh --leiden-resolution 2.0
step 01c.res2 "" $H/01c.run_isograph_gtex.sh --leiden-resolution 2.0
step 01d "00a"  $H/01d.wgcna_gene_brainseq_aging.sh
step 01e "00a"  $H/01e.wgcna_gene_brainseq_sczd.sh
step 01f "00a"  $H/01f.wgcna_gene_gtex.sh
step 01g ""  $H/01g.wgcna_matched_features_brainseq.sh
step 01h ""  $H/01h.wgcna_matched_features_sczd.sh
step 01i ""  $H/01i.wgcna_matched_features_gtex.sh

## 02 -- re-clustering sweeps on the saved graphs, and the module-size tables over every fit
step 02a "01a 01b" $H/02a.sweep_leiden_brainseq.sh
step 02b "01c"     $H/02b.sweep_leiden_gtex.sh
step 02c "01a 01b 01c 01d 01e 01f 01g 01h 01i" $H/02c.module_sizes.sh

dag_finish
