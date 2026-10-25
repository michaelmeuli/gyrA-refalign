#!/usr/bin/env bash
# gyrA of the 72 public genomes (01), of the project isolates (02) and of the public samples (03).
# Output: output/gyrA_refalign/{project,public}/ (kept apart on purpose)
#SBATCH --job-name=gyrA-refalign
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/gyrA_refalign/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/gyrA-refalign
OUTPUT=$OUTPUT_ROOT/gyrA_refalign
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

cd "$REPO"
python scripts/01_make_gyrA_references.py
python scripts/02_call_project.py
python scripts/03_call_public.py

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/project" "$REPO" submit_gyrA.sbatch \
  "$GTDB_COMPLEX_GENOMES" "$GYRA_REFERENCES" "$SCREENING_MAP_RESULTS" "$ASSEMBLIES" "$OUTPUT_ROOT/hsp65_refalign/isolates.tsv"
bash "$MANIFEST_SCRIPT" "$OUTPUT/public" "$REPO" submit_gyrA.sbatch \
  "$GTDB_COMPLEX_GENOMES" "$GYRA_REFERENCES" "$RUNS/public_*/assembly/results"
