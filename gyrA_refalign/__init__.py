import os
from pathlib import Path

SPECIES = ["kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"]

# Mkan329-126 to -133 are reference strains, not patient isolates; they are excluded from every analysis
REFERENCE_STRAINS = frozenset(f"Mkan329-{i:03d}" for i in range(126, 134))

# fixed colour per species (same as mlsa-kansasii) so figures stay comparable
SPECIES_COLORS = {"kansasii": "#1b9e77", "persicum": "#d95f02", "pseudokansasii": "#7570b3", "innocens": "#e7298a",
                  "attenuatum": "#66a61e", "ostraviense": "#e6ab02", "gastri": "#a6761d"}

# Everything lives under the shared kansasii tree; repos/<this repo> sits two levels below it.
ROOT = Path(os.environ.get("KANSASII_ROOT", Path(__file__).resolve().parent.parent.parent.parent))
DATA = ROOT / "data"
GTDB = DATA / "gtdb_genomes" / "Mycobacteriaceae"
# the 72 public GTDB complex genomes (+ .gff), one directory per species
GTDB_GENOMES = GTDB / "mlsa-kansasii"
# gyrA of those genomes (<accession>_gyrA.fasta, 2514-2517 bp CDS) and the probe they were cut out with; made by 01
GYRA_REFERENCES = GTDB / "kansasii_complex_gyrA"
PROBE_ACCESSION = "GCF_002086675.1"  # M. persicum AFPC-000227; its annotated gyrA is the standard 2514 bp CDS
PROBE_FASTA = GYRA_REFERENCES / f"probe_{PROBE_ACCESSION}_gyrA.fasta"
SCREENING_RESULTS = ROOT / "output" / "screening_map_results.csv"
HSP65_ISOLATES = ROOT / "output" / "hsp65_refalign" / "isolates.tsv"
RUNS = ROOT / "runs"
ASSEMBLIES = RUNS / "mkan329" / "assembly" / "results"
PUBLIC_RUNS = ["public_PRJNA1048499", "public_PRJNA1078897", "public_PRJNA323639", "public_PRJNA845115"]
OUTPUT = ROOT / "output" / "gyrA_refalign"
OUTPUT_PROJECT = OUTPUT / "project"  # Mkan329-* isolates only
OUTPUT_PUBLIC = OUTPUT / "public"    # public reference genomes and public WGS runs only
