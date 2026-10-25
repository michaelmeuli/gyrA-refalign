# gyrA-refalign

Is any *M. kansasii* isolate persicum-like (or another subtype) at **gyrA**, the one locus that separates all seven
complex species? Same idea as `hsp65-refalign`, but the sequence comes from the WGS assembly, not a Sanger read.
Project isolates and public data are analysed and written separately.

| Step | Script | Output |
|---|---|---|
| 1 | `scripts/01_make_gyrA_references.py`: gyrA of the 72 public GTDB genomes (probe = *M. persicum* GCF_002086675.1 gyrA CDS laid over each genome), leave-one-out check | `output/gyrA_refalign/public/reference_genomes.tsv`, `REFERENCES_SUMMARY.txt`; references in `data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gyrA/` |
| 2 | `scripts/02_call_project.py`: Mkan329-* assemblies (reference strains 126-133 excluded) | `output/gyrA_refalign/project/isolates.tsv`, `SUMMARY.txt` |
| 3 | `scripts/03_call_public.py`: `runs/public_PRJNA*` samples that GTDB puts in the complex | `output/gyrA_refalign/public/samples.tsv`, `SAMPLES_SUMMARY.txt` |

Run: `sbatch submit_gyrA.sbatch` (~30 min, 8 CPUs, env_mlsa). Status per sample: `ok`, `ambiguous` (second species within 2
differences), `divergent` (< 98% identity to every reference), `low_coverage` (< 95% of the probe placed), `no_hit`.

## Results (2026-10-25)
- References: 72/72 genomes give a gyrA; leave-one-out 72/72 correct; the 3 Korean kansasii genomes (persicum-like hsp65) are kansasii at gyrA.
- Project: gyrA species = WGS species in 171/171 callable uncontaminated isolates; no kansasii isolate is persicum-like at gyrA (28 persicum-like are all WGS persicum). Not callable: -010, -021, -055, -113 (low coverage; -010/-021/-113 contaminated), -184 (WGS *M. avium*, divergent).
- Public: 1026/1033 callable samples agree with GTDB. 7 disagree, all PRJNA1078897: 5 GTDB-persicum samples have a kansasii gyrA (identical to a kansasii reference, ANI to persicum 99.5%), 2 GTDB-kansasii samples (SRR28038314, SRR28038483) have a persicum gyrA although their hsp65 is kansasii. Possible mosaic/recombinant or mixed cultures; not investigated.
- 6 public samples are too fragmented for gyrA (low_coverage).

## Notes
- The mlsa-kansasii gyrA extraction is not used: it returns a 3774 bp gene model for kansasii/ostraviense/gastri/innocens (1260 bp longer than the 2514 bp of the others).
- `gyrA_refalign/refalign.py` is a vendored copy of hsp65-refalign's; keep matching logic in sync by hand.
- Reads `data/imm/screening_map_results.csv` (`SCREENING_MAP_RESULTS`) like every repo; the copy in `output/` is for download only.
- Layout: `gyrA_refalign/` (package, vendored `kansasii_paths.py`), `scripts/01-03`, `submit_gyrA.sbatch`. Paths are flags: `--complex-genomes`, `--gyra-references`, `--assemblies`, `--screening-map-results`, `--hsp65-isolates`, `--out`.
- The samples that differ from their species label at gyrA (or at hsp65) are collected in `output/atypical_samples/` and explained in `kansasii-docs/docs/ATYPICAL.md`.
