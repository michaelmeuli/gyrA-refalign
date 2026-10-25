"""Repo-specific constants; project-wide paths, species and colours are in `kansasii_paths` (shared, vendored)."""
from .kansasii_paths import GYRA_REFERENCES, OUTPUT_ROOT

PROBE_ACCESSION = "GCF_002086675.1"  # M. persicum AFPC-000227; its annotated gyrA is the standard 2514 bp CDS
PROBE_FASTA = GYRA_REFERENCES / f"probe_{PROBE_ACCESSION}_gyrA.fasta"
HSP65_ISOLATES = OUTPUT_ROOT / "hsp65_refalign" / "isolates.tsv"  # output of hsp65-refalign 01 (hsp65 call per isolate)
OUTPUT = OUTPUT_ROOT / "gyrA_refalign"
OUTPUT_PROJECT = OUTPUT / "project"  # Mkan329-* isolates only
OUTPUT_PUBLIC = OUTPUT / "public"    # public reference genomes and public WGS runs only
