#!/usr/bin/env bash
## QUEST ONLY (PI item 10a): submit 01l -> 02f (array) -> 03b -> 04a as a dependency chain for one
## region. The Bridges-2 DAG (run_stage.sh) lists these as manual steps and never submits
## them. The region's BAMs must be staged and listed under quest.ase_junction.regions in
## configs/data_sources.yaml.
##
##   bash 06_switch_mechanism/_h/ase_junction_quest.sh --region caudate
##   bash 06_switch_mechanism/_h/ase_junction_quest.sh --region caudate --dry-run
set -euo pipefail

cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
H=06_switch_mechanism/_h
mkdir -p 06_switch_mechanism/_m/logs

REGION=dlpfc
DRY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --region) REGION=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

## The array is sized from the manifest, an upper bound on the samples 01l keeps (it drops
## samples without a BAM, index or phASER VCF, or failing quickcheck); tasks past the end
## of samples.tsv exit cleanly.
N=$(python -c "
import pandas as pd
from isograph_benchmark.real_data.ase_junction_switch import quest_config
print(len(pd.read_csv(quest_config('${REGION}')['manifest'], sep='\t', usecols=['sample_id'])))")
ARRAY="1-${N}%60"

if (( DRY )); then
    echo "sbatch ${H}/01l.ase_junction_targets.sh --region ${REGION}"
    echo "sbatch --array=${ARRAY} --dependency=afterok:<01l> ${H}/02f.ase_junction_count.sh --region ${REGION}"
    echo "sbatch --dependency=afterok:<02f> ${H}/03b.ase_junction_screen.sh --region ${REGION}"
    echo "sbatch --dependency=afterok:<03b> ${H}/04a.ase_junction_allelic.sh --region ${REGION}"
    exit 0
fi

T=$(sbatch --parsable --job-name="ase-junc-targets-${REGION}" \
    "${H}/01l.ase_junction_targets.sh" --region "${REGION}")
C=$(sbatch --parsable --job-name="ase-junc-count-${REGION}" --array="${ARRAY}" \
    --dependency="afterok:${T}" "${H}/02f.ase_junction_count.sh" --region "${REGION}")
S=$(sbatch --parsable --job-name="ase-junc-screen-${REGION}" \
    --dependency="afterok:${C}" "${H}/03b.ase_junction_screen.sh" --region "${REGION}")
A=$(sbatch --parsable --job-name="ase-junc-allelic-${REGION}" \
    --dependency="afterok:${S}" "${H}/04a.ase_junction_allelic.sh" --region "${REGION}")
echo -e "step\tjob_id\tregion\n01l\t${T}\t${REGION}\n02f\t${C}\t${REGION}\n03b\t${S}\t${REGION}\n04a\t${A}\t${REGION}"
