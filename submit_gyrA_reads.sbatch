#!/usr/bin/env bash
# Read-level check of the public samples whose assembly gyrA disagrees with GTDB (04): reads mapped to gyrA, hsp65, rpoB, recA, tuf.
# Needs bwa, samtools and mafft (snippy env); the script itself uses the standard library only.
# Output: output/gyrA-refalign/results/public_reads/
#SBATCH --job-name=gyrA-reads
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
export MAFFT_BINARIES="$CONDA_PREFIX/libexec/mafft"  # the conda mafft wrapper otherwise looks in /usr/local and fails with "v0.000 != v7.520"

cd "$REPO"
python scripts/04_gyra_read_check.py

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT/public_reads" "$REPO" submit_gyrA_reads.sbatch \
  "$OUTPUT/public/samples.tsv" "$OUTPUT/public/reference_genomes.tsv" "$GYRA_REFERENCES" \
  "$OUT_MLSA/results/locus_discovery/alignments" "$PUBLIC_READS/PRJNA1078897/reads"
