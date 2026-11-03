#!/usr/bin/env bash
# Read-level hsp65 amplicon (441 bp) of the public samples whose assembly gyrA disagrees with GTDB (05): reads -> consensus -> 72 amplicons.
# Needs bwa (snippy env); run after submit_gyrA_reads.sbatch (04), whose read_check.tsv lists the samples.
# Output: output/gyrA-refalign/results/public_reads/
#SBATCH --job-name=hsp65-reads
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/gyrA-refalign/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/gyrA-refalign
OUTPUT=$OUT_GYRA/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snippy

cd "$REPO"
python scripts/05_hsp65_read_check.py

# provenance; MANIFEST_NAME keeps it apart from submit_gyrA_reads's MANIFEST.txt in the same directory
MANIFEST_NAME=MANIFEST_hsp65_reads.txt bash "$MANIFEST_SCRIPT" "$OUTPUT/public_reads" "$REPO" submit_hsp65_reads.sbatch \
  "$OUTPUT/public_reads/read_check.tsv" "$OUTPUT/public/samples.tsv" "$OUTPUT/public/reference_genomes.tsv" \
  "$HSP65_AMPLICONS" "$PUBLIC_READS/PRJNA1078897/reads"
