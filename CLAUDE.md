# CLAUDE.md

See README.md. Easy to get wrong:

- **Project and public stay separate:** `02` writes only `output/gyrA-refalign/results/project/`, `01` and `03` only `public/`, `04` (read-level check of the public gyrA disagreements) only `public_reads/`. Do not merge them into one table.
- **References are re-extracted** by laying the persicum gyrA CDS over each genome (`gyrA_refalign/gyra.py`); do not reuse `locus_discovery/alignments/gyrA.raw.fasta` (inconsistent start calls).
- **Paths:** `gyrA_refalign/kansasii_paths.py` is a vendored copy of `kansasii-docs/templates/kansasii_paths.py` (never edit it here; `bash repos/check_sync.sh --fix`). Repo-specific constants (`PROBE_ACCESSION`, `OUTPUT_*`) are in `gyrA_refalign/__init__.py`.
- **Independent of other repos:** `refalign.py` is vendored from hsp65-refalign; do not import `mlsa` or `hsp65_refalign`.
- **Environment:** `env_mlsa`; run via `sbatch submit_gyrA.sbatch`, not the login node. Type check: `~/data/conda/envs/env_typecheck/bin/mypy` (strict).
- **Public samples with CheckM contamination above `PUBLIC_MAX_CONTAMINATION` (5 %) are left out of `03`, `04` and `05`** (`results/public/excluded_contaminated.tsv`); samples without a CheckM result are kept. Project isolates are not filtered (flagged isolates stay in, as in hsp65-refalign).
