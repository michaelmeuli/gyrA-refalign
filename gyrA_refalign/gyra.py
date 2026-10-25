"""gyrA of an assembly by laying a fixed gyrA CDS (the probe) over it, and the species call on that sequence.

The probe is a 2514 bp gyrA CDS. Laying it over a genome yields that genome's gyrA in probe coordinates, so the
same call works for GTDB genomes, project assemblies and public assemblies without any annotation. (The existing
mlsa-kansasii gyrA extraction cannot be used as references: for kansasii, ostraviense, gastri and innocens it
returns a 3774/3777 bp gene model, 1260 bp longer than in the other species, because of a different start call.)
"""
from __future__ import annotations

import hashlib
import multiprocessing
from pathlib import Path

from Bio import SeqIO

from gyrA_refalign import PROBE_FASTA
from gyrA_refalign.kansasii_paths import GYRA_REFERENCES
from gyrA_refalign.refalign import Reference, align_to_reference, identify_read, load_references

MIN_COVERAGE = 0.95  # fraction of probe columns that must be ACGT in the assembly for a call


def read_probe(probe_fasta: Path = PROBE_FASTA) -> str:
    return str(next(SeqIO.parse(probe_fasta, "fasta")).seq).upper()


def extract_gyra(fasta: Path, probe: str) -> tuple[str, float] | None:
    """(gyrA in probe orientation, fraction of probe columns covered by ACGT), None if the probe does not place."""
    hit = align_to_reference(Reference("own", "own", "", fasta), probe)
    if hit is None:
        return None
    seq = hit.ref_cols.replace("-", "")
    covered = sum(c in "ACGT" for c in hit.ref_cols) / len(probe)
    return seq, covered


def load_gyra_references(gyra_references: Path = GYRA_REFERENCES) -> list[Reference]:
    refs = load_references(gyra_references, None)
    return [r for r in refs if not r.path.name.startswith("probe_")]


def _hap(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()[:8]


def call_gyra(name: str, fasta: Path, probe: str, refs: list[Reference]) -> dict[str, object]:
    """Species call of one assembly's gyrA; one row for the output tables."""
    row: dict[str, object] = {"gyrA_bp": 0, "gyrA_coverage": 0.0, "status": "no_hit", "closest_species": "NA",
                              "closest_accession": "", "atypical": False, "diffs": None, "identity": None,
                              "second_species": "", "margin": None, "persicum_like": False, "haplotype": "",
                              "exact_match": ""}
    ext = extract_gyra(fasta, probe)
    if ext is None:
        return row
    seq, cov = ext
    row["gyrA_bp"], row["gyrA_coverage"] = len(seq), round(cov, 4)
    if cov < MIN_COVERAGE:
        row["status"] = "low_coverage"
        return row
    res = identify_read(refs, name, "gyrA", seq)
    b = res.best
    if b is None:
        return row
    row["status"] = res.status
    row["closest_species"] = res.closest_species
    row["closest_accession"] = res.hits[b].ref.accession
    row["atypical"] = res.atypical
    row["diffs"] = res.diffs[b]
    row["identity"] = round(res.identity(b), 4)
    r = res.runner_up
    if r is not None:
        row["second_species"] = res.hits[r].ref.species
        row["margin"] = res.diffs[r] - res.diffs[b]
    row["persicum_like"] = res.closest_species == "persicum"
    row["haplotype"] = _hap(seq)
    exact = [h.ref.accession for h, d in zip(res.hits, res.diffs) if d == 0]
    row["exact_match"] = ";".join(exact[:3]) + (f";+{len(exact) - 3}" if len(exact) > 3 else "")
    return row


_probe = ""
_refs: list[Reference] = []


def _init(gyra_references: Path) -> None:
    global _probe, _refs
    _probe, _refs = read_probe(gyra_references / PROBE_FASTA.name), load_gyra_references(gyra_references)


def _work(item: tuple[str, Path]) -> tuple[str, dict[str, object]]:
    return item[0], call_gyra(item[0], item[1], _probe, _refs)


def call_many(items: list[tuple[str, Path]], ncpu: int,
              gyra_references: Path = GYRA_REFERENCES) -> dict[str, dict[str, object]]:
    """call_gyra for (name, assembly) pairs on ncpu processes."""
    with multiprocessing.Pool(ncpu, initializer=_init, initargs=(gyra_references,)) as pool:
        return dict(pool.imap_unordered(_work, items, chunksize=4))
