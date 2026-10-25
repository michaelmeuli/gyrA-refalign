"""Shared constants of the kansasii project: paths, species, colours, sample identifiers.

CANONICAL COPY: kansasii-docs/templates/kansasii_paths.py. Every repo vendors this file unchanged
(next to its code, so it imports as `kansasii_paths`); `repos/check_sync.sh` flags drift and `--fix` restores it.
Edit the template, never a vendored copy. The shell twin with the same names is templates/kansasii_paths.sh.

Rule: a path that is the same in two repos has the same constant name in all of them. A script that needs one
of these paths takes it as a command-line flag whose default is the constant (e.g. `--gtdb-representatives`),
so nothing below is hard-wired into the logic. Scripts read from DATA and RUNS; OUTPUT_ROOT receives results.
Override the project root with the KANSASII_ROOT environment variable.
"""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------------------------------------- root
_DEFAULT_ROOT = Path("/shares/sander.imm.uzh/MM/kansasii")


def _find_root() -> Path:
    """KANSASII_ROOT, else the nearest parent of this file that holds `repos/` and `data/`, else the shares path."""
    env = os.environ.get("KANSASII_ROOT")
    if env:
        return Path(env)
    for parent in Path(__file__).resolve().parents:
        if (parent / "repos").is_dir() and (parent / "data").is_dir():
            return parent
    return _DEFAULT_ROOT


ROOT = _find_root()
REPOS = ROOT / "repos"
DATA = ROOT / "data"
RUNS = ROOT / "runs"
OUTPUT_ROOT = ROOT / "output"  # one directory per repo or analysis below it; a repo's own is `OUTPUT`
MANIFEST_SCRIPT = REPOS / "kansasii-docs" / "bin" / "write_manifest.sh"

# ---------------------------------------------------------------------------------------------- reference genomes
GTDB_DIR = DATA / "gtdb_genomes" / "Mycobacteriaceae"
NCBI_DATASET = GTDB_DIR / "ncbi_dataset" / "data"  # <accession>/ downloaded by immensekansasii/scripts/01
# one GTDB representative genome per complex species, <accession>_query.fasta (references of plain read alignment)
GTDB_REPRESENTATIVES = GTDB_DIR / "kansasii_complex_gtdb_representatives"
# the 72 public GTDB genomes (+ .gff) of the complex, one directory per species (mlsa-kansasii/scripts/02)
GTDB_COMPLEX_GENOMES = GTDB_DIR / "mlsa-kansasii"
# in-silico hsp65 amplicon (441 bp) of each of the 72 genomes, <accession>_hsp65.fasta (mlsa-kansasii/scripts/06)
HSP65_AMPLICONS = GTDB_DIR / "kansasii_complex_hsp65_amplicons"
# gyrA CDS of each of the 72 genomes, <accession>_gyrA.fasta (gyrA-refalign/scripts/01)
GYRA_REFERENCES = GTDB_DIR / "kansasii_complex_gyrA"
# GTDB release tables (downloaded by immensekansasii/scripts/01 and mlsa-kansasii/scripts/01) and the accession lists made from them
GTDB_RELEASE = "232"
GTDB_METADATA_DIR = DATA / "immensekansasii" / "gtdb" / f"gtdb{GTDB_RELEASE}"
GTDB_ACCESSIONS_DIR = OUTPUT_ROOT / "immensekansasii" / "gtdb" / f"gtdb{GTDB_RELEASE}"
# M. kansasii ATCC 12478, the snippy / cgMLST / marker-screen reference
ATCC12478_ACCESSION = "GCF_000157895.3"
ATCC12478_DIR = NCBI_DATASET / ATCC12478_ACCESSION
ATCC12478_FASTA = ATCC12478_DIR / f"{ATCC12478_ACCESSION}_ASM15789v2_genomic.fna"
ATCC12478_GFF = ATCC12478_DIR / "genomic.gff"

# -------------------------------------------------------------------------------------------- lab tables and reads
IMM_DIR = DATA / "imm"
MIC_CSV = IMM_DIR / "mic.csv"
SCREENING_MAP_LINK = IMM_DIR / "screening_map_link.csv"  # immensekansasii/scripts/02
SCREENING_MAP_RESULTS = IMM_DIR / "screening_map_results.csv"  # immensekansasii/scripts/06
SCREENING_MAP_AST_JC = IMM_DIR / "screening_map_ast_jc.csv"
SANGER_LINK_DIR = DATA / "sanger" / "seq_kansasii"  # .ab1 symlinks (mlsa-kansasii/scripts/03)
ILLUMINA_DIR = DATA / "illumina"
PUBLIC_READS = ILLUMINA_DIR / "public"  # <BioProject>/metadata.tsv, selected.tsv and reads (immensekansasii/scripts/09-11)

# ----------------------------------------------------------------------------------------- pipeline (runs/) results
ASSEMBLIES = RUNS / "mkan329" / "assembly" / "results"  # <PROBENNUMMER>/1_unicycler/<PROBENNUMMER>.fasta
PUBLIC_PROJECTS = ["PRJNA1078897", "PRJNA1048499", "PRJNA845115", "PRJNA323639"]


def public_assemblies(project: str) -> Path:
    """Pipeline results of a public BioProject: <sample>/1_unicycler/<sample>.fasta, 3_quality/GTDB/..."""
    return RUNS / f"public_{project}" / "assembly" / "results"


# ------------------------------------------------------------------------------------ protocols and external software
HSP65_PROTOCOL = DATA / "protocols" / "32.43505-4-Mykobakterien NTM 65kDa.docx"  # TB-11 / TB-12(w) primers
# not part of the project: pipeline databases and the container cache of the IMMENSE pipeline
IMMENSE_DEPENDENCIES = Path("/shares/sander.imm.uzh/software/pipelines/IMMense/IMMense_dependencies")
IMMENSE_CONTAINERS = IMMENSE_DEPENDENCIES / "containers"
IMMENSE_DB_ROOT = IMMENSE_DEPENDENCIES / "databases" / "kansasii_complex"  # snippy databases per run

# ------------------------------------------------------------------------------------------- species and samples
# order of the species in matrices and legends; the 7th, gastri, is public references only (no project isolate)
SPECIES = ["kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"]
# species that count as the complex in cohort analyses (MST, cgMLST, markers)
COMPLEX = ["kansasii", "persicum", "pseudokansasii", "attenuatum", "innocens", "ostraviense"]
# fixed colour per species, so figures of all repos stay comparable
SPECIES_COLORS = {"kansasii": "#1b9e77", "persicum": "#d95f02", "pseudokansasii": "#7570b3", "innocens": "#e7298a",
                  "attenuatum": "#66a61e", "ostraviense": "#e6ab02", "gastri": "#a6761d"}
# Mkan329-126 to -133 are control / reference strains, not patient isolates
REFERENCE_STRAINS = frozenset(f"Mkan329-{i:03d}" for i in range(126, 134))
# every column of screening_map_link.csv that holds a TNR of the row's isolate (TNR is the primary one)
TNR_COLUMNS = ("TNR", "TNR_JC", "TNR_MM", "TNR3", "TNR4", "TNR5", "TNR6")
# screening_map_results.csv adds TNR_MLSA, a TNR found only through its Sanger read (immensekansasii/scripts/06)
TNR_COLUMNS_RESULTS = TNR_COLUMNS + ("TNR_MLSA",)
# the three Korean M. kansasii GTDB genomes with a persicum-like hsp65 (docs/ATYPICAL.md)
ATYPICAL_HSP65_GENOMES = ("GCF_002705785.1", "GCF_002705825.1", "GCF_002705865.1")
