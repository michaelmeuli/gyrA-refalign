"""Small helpers shared by the read-level check scripts (04, 05): FASTA/TSV reading and external tools."""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    """Run a tool and keep its stderr in the error, so a failed bwa or mafft call says why."""
    res = subprocess.run(cmd, text=True, capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd[:3])} ... exited {res.returncode}: {res.stderr.strip()[:500]}")
    return res


def map_reads(ref: Path, r1: Path, r2: Path, threads: int = 4) -> str:
    """SAM text of the reads that map to `ref` (bwa mem | samtools view -F 0x904 -q 20).

    The filter runs while bwa streams: a whole-genome read set is mostly unmapped against a single locus and
    would not fit in memory as SAM text."""
    bwa = subprocess.Popen(["bwa", "mem", "-t", str(threads), str(ref), str(r1), str(r2)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    view = subprocess.run(["samtools", "view", "-h", "-F", "0x904", "-q", "20", "-"], stdin=bwa.stdout, capture_output=True, text=True)
    if bwa.stdout:
        bwa.stdout.close()
    if bwa.wait() != 0 or view.returncode != 0:
        raise RuntimeError(f"bwa mem | samtools view failed (bwa {bwa.returncode}, samtools {view.returncode}): {view.stderr.strip()[:300]}")
    return view.stdout
