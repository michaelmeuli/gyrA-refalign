# gyrA-refalign

Is any *M. kansasii* isolate persicum-like (or another subtype) at **gyrA**, the one locus that separates all seven
complex species? Same idea as `hsp65-refalign`, but the sequence comes from the WGS assembly, not a Sanger read.
Project isolates and public data are analysed and written separately.

| Step | Script | Output |
|---|---|---|
| 1 | `scripts/01_make_gyrA_references.py`: gyrA of the 72 public GTDB genomes (probe = *M. persicum* GCF_002086675.1 gyrA CDS laid over each genome), leave-one-out check | `output/gyrA-refalign/results/public/reference_genomes.tsv`, `REFERENCES_SUMMARY.txt`; references in `data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gyrA/` |
| 2 | `scripts/02_call_project.py`: Mkan329-* assemblies (reference strains 126-133 excluded) | `output/gyrA-refalign/results/project/isolates.tsv`, `SUMMARY.txt` |
| 3 | `scripts/03_call_public.py`: `runs/public_PRJNA*` samples that GTDB puts in the complex | `output/gyrA-refalign/results/public/samples.tsv`, `SAMPLES_SUMMARY.txt` |
| 4 | `scripts/04_gyra_read_check.py` (`sbatch submit_gyrA_reads.sbatch`, snippy env): raw reads of the public samples whose gyrA disagrees with GTDB (+ 2 concordant controls) mapped to gyrA, rpoB, recA, tuf; allele counts at sites fixed between kansasii and persicum | `output/gyrA-refalign/results/public_reads/read_check.tsv`, `diagnostic_sites.tsv`, `read_haplotypes.tsv`, `READ_CHECK_SUMMARY.txt` |
| 5 | `scripts/05_hsp65_read_check.py` (`sbatch submit_hsp65_reads.sbatch`): read consensus of the 441 bp hsp65 amplicon compared with the 72 amplicons | `public_reads/hsp65_amplicon_check.tsv`, `hsp65_consensus.fasta`, `HSP65_AMPLICON_SUMMARY.txt` |

Run: `sbatch submit_gyrA.sbatch` (~30 min, 8 CPUs, env_mlsa). Status per sample: `ok`, `ambiguous` (second species within 2
differences), `divergent` (< 98% identity to every reference), `low_coverage` (< 95% of the probe placed), `no_hit`.

## Results (2026-10-25)
- References: 72/72 genomes give a gyrA; leave-one-out 72/72 correct; the 3 Korean kansasii genomes (persicum-like hsp65) are kansasii at gyrA.
- Project: gyrA species = WGS species in 171/171 callable uncontaminated isolates; no kansasii isolate is persicum-like at gyrA (28 persicum-like are all WGS persicum). Not callable: -010, -021, -055, -113 (low coverage; -010/-021/-113 contaminated), -184 (WGS *M. avium*, divergent).
- Public (2026-11-03): public samples with CheckM contamination above 5 % (27 of 1039) are left out (`results/public/excluded_contaminated.tsv`); 29 samples without a CheckM result are kept. Of the 1012 left, 1006 are callable and 1000 agree with GTDB. 6 disagree, all PRJNA1078897: 4 GTDB-persicum samples have a kansasii gyrA (identical to a kansasii reference, ANI to persicum 99.5%), 2 GTDB-kansasii samples (SRR28038314, SRR28038483) have a persicum gyrA although their hsp65 is kansasii.
- 6 public samples are too fragmented for gyrA (low_coverage).
- Public samples are filtered by `PUBLIC_MAX_CONTAMINATION` (5 %, CheckM `3_quality/checkM/checkm_output.tsv`) in 03, 04 and 05; project isolates are not (see `hsp65-refalign/CLAUDE.md`).

### hsp65 amplicon from the reads (05)
Differences (of 441 bp) between the read consensus and the closest reference amplicon of each group; contaminated samples are left out.

| Sample(s) | GTDB | Typical kansasii | Korean-type | Persicum | hsp65 call |
|---|---|---|---|---|---|
| SRR28038242, -316, -317, -318 | persicum | 10 | 0 | 1 | persicum with the 7 reference genomes; ambiguous (margin 1) with the 72 |
| SRR28038314, -483 | kansasii | 0 | 11 | 12 | kansasii, margin 8 |
| SRR28038241 (control) | kansasii | 0 | 11 | 12 | kansasii, margin 8 |
| SRR28038325 (control) | persicum | 12 | 2 | 0 | persicum, margin 2 |
| SRR28038296 (known Korean-type) | kansasii | 10 | 0 | 1 | ambiguous (Korean-type, as expected) |

## Notes
- hsp65 is not in 04 (only one site is fixed between kansasii and persicum once the Korean-type genomes are left out); 05 compares the read consensus over the whole amplicon instead.
- The mlsa-kansasii gyrA extraction is not used: it returns a 3774 bp gene model for kansasii/ostraviense/gastri/innocens (1260 bp longer than the 2514 bp of the others).
- `gyrA_refalign/refalign.py` is a vendored copy of hsp65-refalign's; keep matching logic in sync by hand.
- Reads `data/imm/screening_map_results.csv` (`SCREENING_MAP_RESULTS`) like every repo; the copy in `output/` is for download only.
- Layout: `gyrA_refalign/` (package, vendored `kansasii_paths.py`), `scripts/01-05`, `submit_gyrA.sbatch`, `submit_gyrA_reads.sbatch`, `submit_hsp65_reads.sbatch`. Paths are flags: `--complex-genomes`, `--gyra-references`, `--assemblies`, `--screening-map-results`, `--hsp65-isolates`, `--out`.
- The samples that differ from their species label at gyrA (or at hsp65) are collected in `output/kansasii-docs/results/atypical_samples/` and explained in `kansasii-docs/docs/ATYPICAL.md`.
