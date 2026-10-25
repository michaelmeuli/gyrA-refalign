#!/usr/bin/env bash
#SBATCH --job-name=gyrA-refalign
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/gyrA_refalign/logs/%x-%j.log

set -euo pipefail

K="/shares/sander.imm.uzh/MM/kansasii"
REPO_ROOT="$K/repos/gyrA-refalign"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

cd "$REPO_ROOT"
python 01_make_gyrA_references.py
python 02_call_project.py
python 03_call_public.py

for d in project public; do
  bash "$K/repos/kansasii-docs/bin/write_manifest.sh" "$K/output/gyrA_refalign/$d" "$REPO_ROOT" "submit_gyrA.sbatch" \
    "$K/data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii" "$K/output/screening_map_results.csv" \
    "$K/runs/mkan329/assembly/results" "$K/runs/public_PRJNA1048499/assembly/results"
done
