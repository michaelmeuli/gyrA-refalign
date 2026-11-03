# gyrA-refalign

Is any *M. kansasii* isolate persicum-like (or another subtype) at **gyrA**, the one locus that separates all seven
complex species? Same idea as `hsp65-refalign`, but the sequence comes from the WGS assembly, not a Sanger read.
Project isolates and public data are analysed and written separately.

| Step | Script | Output |
|---|---|---|
| 1 | `scripts/01_make_gyrA_references.py`: gyrA of the 72 public GTDB genomes (probe = *M. persicum* GCF_002086675.1 gyrA CDS laid over each genome), leave-one-out check | `output/gyrA-refalign/results/public/reference_genomes.tsv`, `REFERENCES_SUMMARY.txt`; references in `data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_gyrA/` |
| 2 | `scripts/02_call_project.py`: Mkan329-* assemblies (reference strains 126-133 excluded) | `output/gyrA-refalign/results/project/isolates.tsv`, `SUMMARY.txt` |
| 3 | `scripts/03_call_public.py`: `runs/public_PRJNA*` samples that GTDB puts in the complex | `output/gyrA-refalign/results/public/samples.tsv`, `SAMPLES_SUMMARY.txt` |
| 4 | `scripts/04_gyra_read_check.py` (`sbatch submit_gyrA_reads.sbatch`, snippy env): raw reads of the public samples whose gyrA disagrees with GTDB (+ 2 concordant controls) mapped to gyrA, rpoB, recA, tuf (and hsp65, see note); allele counts at sites fixed between kansasii and persicum | `output/gyrA-refalign/results/public_reads/read_check.tsv`, `diagnostic_sites.tsv`, `read_haplotypes.tsv`, `READ_CHECK_SUMMARY.txt` |
| 5 | `scripts/05_hsp65_read_check.py` (`sbatch submit_hsp65_reads.sbatch`): read consensus of the 441 bp hsp65 amplicon compared with the 72 amplicons | `public_reads/hsp65_amplicon_check.tsv`, `hsp65_consensus.fasta`, `HSP65_AMPLICON_SUMMARY.txt` |

Run: `sbatch submit_gyrA.sbatch` (~30 min, 8 CPUs, env_mlsa). Status per sample: `ok`, `ambiguous` (second species within 2
differences), `divergent` (< 98% identity to every reference), `low_coverage` (< 95% of the probe placed), `no_hit`.

## Results (2026-10-25)
- References: 72/72 genomes give a gyrA; leave-one-out 72/72 correct; the 3 Korean kansasii genomes (persicum-like hsp65) are kansasii at gyrA.
- Project: gyrA species = WGS species in 171/171 callable uncontaminated isolates; no kansasii isolate is persicum-like at gyrA (28 persicum-like are all WGS persicum). Not callable: -010, -021, -055, -113 (low coverage; -010/-021/-113 contaminated), -184 (WGS *M. avium*, divergent).
- Public: 1026/1033 callable samples agree with GTDB. 7 disagree, all PRJNA1078897: 5 GTDB-persicum samples have a kansasii gyrA (identical to a kansasii reference, ANI to persicum 99.5%), 2 GTDB-kansasii samples (SRR28038314, SRR28038483) have a persicum gyrA although their hsp65 is kansasii. Read-level check (04, 05; 2026-11-03): in all 7 the reads carry the other species' gyrA (99.6-99.98 % of bases at 182 sites), with no persicum-kansasii mixture and no breakpoint inside the gene, while rpoB, recA and tuf carry the GTDB species, so these are gyrA mosaics. The 5 persicum samples share two gyrA haplotypes (SRR28038252/-316/-318 and -242/-317), i.e. probably two events. hsp65 (441 bp read consensus): the 5 persicum samples are identical to the Korean-type amplicon (1 diff to persicum, 10 to typical kansasii; persicum with the 7 representatives, `ambiguous` with the 72), the 2 kansasii samples are typical kansasii (margin 8). SRR28038252's assembly is contaminated (CheckM 100 %, 10.5 Mb) but its reads show no second complex species at gyrA (minor alleles 1-3 % at rpoB/recA/tuf). Short reads only; no long-read confirmation.
- 6 public samples are too fragmented for gyrA (low_coverage).

## Notes
- In 04's `read_check.tsv` the `hsp65` rows rest on a single fixed kansasii/persicum site (the Korean-type genomes are left out) and read "kansasii" for the Korean-type persicum samples; use 05 for hsp65.
- The mlsa-kansasii gyrA extraction is not used: it returns a 3774 bp gene model for kansasii/ostraviense/gastri/innocens (1260 bp longer than the 2514 bp of the others).
- `gyrA_refalign/refalign.py` is a vendored copy of hsp65-refalign's; keep matching logic in sync by hand.
- Reads `data/imm/screening_map_results.csv` (`SCREENING_MAP_RESULTS`) like every repo; the copy in `output/` is for download only.
- Layout: `gyrA_refalign/` (package, vendored `kansasii_paths.py`), `scripts/01-05`, `submit_gyrA.sbatch`, `submit_gyrA_reads.sbatch`, `submit_hsp65_reads.sbatch`. Paths are flags: `--complex-genomes`, `--gyra-references`, `--assemblies`, `--screening-map-results`, `--hsp65-isolates`, `--out`.
- The samples that differ from their species label at gyrA (or at hsp65) are collected in `output/kansasii-docs/results/atypical_samples/` and explained in `kansasii-docs/docs/ATYPICAL.md`.
