"""Species identification of a Sanger read by plain alignment to one reference
genome per species (no MLSA, no locus extraction from the references).

Each read is seeded against every reference genome with exact 15-mers (both
strands), the best-supported genomic window is aligned end to end with
Biopython's PairwiseAligner, and the reference that differs at the fewest
positions of the read is the closest species.

Vendored from hsp65-refalign (hsp65_refalign/refalign.py, itself vendored from mlsa-kansasii) so that this repo
does not depend on either. Matching logic is unchanged; here the "read" is a gyrA sequence taken from an assembly.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import numpy.typing as npt
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

K = 15
# GTDB kansasii genomes with a persicum-like hsp65 (11 SNPs in the 441 bp amplicon) that are
# M. kansasii by ANI and at every other locus (LIT.md). hsp65 alone would call them persicum, so a read
# matching them is reported as kansasii with an "atypical hsp65" note, not silently as persicum.
# The reference strain of each species (= the 7 genomes of kansasii_complex_gtdb_representatives): ATCC 12478,
# AFPC-000227, MK142, MK13, MK41, 241/15, DSM 43505. Figures show these, whatever the larger reference set is.
REPRESENTATIVES = frozenset({"GCF_000157895.3", "GCF_002086675.1", "GCF_900566075.1", "GCF_900566055.1",
                             "GCF_900566085.1", "GCF_002705925.1", "GCF_002102175.1"})
ATYPICAL_HSP65 = frozenset({"GCF_002705785.1", "GCF_002705825.1", "GCF_002705865.1"})
PAD = 40  # bases of reference added on each side of the seeded window
MIN_SEEDS = 5
SPECIES_ORDER = ["kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"]
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _LUT[_b + 32] = _i
_SENTINEL = np.uint32(0xFFFFFFFF)


def _kmer_codes(seq: str | bytes) -> npt.NDArray[np.uint32]:
    """2-bit code of every K-mer (uint32; sentinel where it holds a non-ACGT base)."""
    raw = seq.encode() if isinstance(seq, str) else seq
    b = _LUT[np.frombuffer(raw, dtype=np.uint8)]
    n = len(b) - K + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint32)
    bad = np.cumsum(np.concatenate([[0], b == 255]))
    invalid = (bad[K:] - bad[:-K]) > 0
    code = np.zeros(n, dtype=np.uint32)
    for j in range(K):
        code = (code << np.uint32(2)) | (b[j:n + j] & np.uint8(3)).astype(np.uint32)
    code[invalid] = _SENTINEL
    return code


_COMP = str.maketrans("ACGTN-", "TGCAN-")


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


_NAME_RE = re.compile(
    r"^\S+ Mycobacterium(?:\s+(?P<species>[a-z]+))(?:\s+(?:strain|subsp\.\s+\S+)?\s*(?P<strain>.+?))?"
    r"(?:\s+(?:isolate\b|plasmid\b|chromosome\b|NODE|contig|Scaffold|scaffold).*|,.*)?$")


@dataclass
class Reference:
    accession: str
    species: str
    strain: str
    path: Path
    atypical: bool = False
    contigs: dict[str, str] = field(default_factory=dict, repr=False)
    _ids: list[str] = field(default_factory=list, repr=False)
    # offset of each contig in _codes
    _starts: npt.NDArray[np.intp] = field(default_factory=lambda: np.empty(0, dtype=np.intp), repr=False)
    # all contigs' k-mer codes, sentinel-separated
    _codes: npt.NDArray[np.uint32] = field(default_factory=lambda: np.empty(0, dtype=np.uint32), repr=False)

    @property
    def label(self) -> str:
        return f"Mycobacterium {self.species} {self.strain}{' (atypical hsp65)' if self.atypical else ''}".strip()

    def load(self) -> None:
        if not self.contigs:
            self.contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(self.path, "fasta")}
            gap = np.full(K, _SENTINEL, dtype=np.uint32)
            parts: list[npt.NDArray[np.uint32]] = []
            starts: list[int] = []
            pos = 0
            self._ids = []
            for cid, seq in self.contigs.items():
                c = _kmer_codes(seq)
                parts += [c, gap]
                self._ids.append(cid)
                starts.append(pos)
                pos += len(c) + K
            self._starts = np.array(starts, dtype=np.intp)
            self._codes = np.concatenate(parts)


def load_references(ref_dir: Path, species_order: list[str] | None = None) -> list[Reference]:
    """One Reference per *.fasta/*.fna/*.fa in ref_dir, named from the first header."""
    refs = []
    for path in sorted(p for p in ref_dir.iterdir() if p.suffix in {".fasta", ".fna", ".fa"}):
        with open(path, encoding="utf-8") as fh:
            header = fh.readline().strip().lstrip(">")
        accession = re.match(r"(GC[AF]_\d+\.\d)", path.name)
        m = _NAME_RE.match(header)
        acc = accession.group(1) if accession else path.stem
        refs.append(Reference(
            accession=acc,
            species=m.group("species") if m else "unknown",
            strain=(m.group("strain") or "").strip() if m else header[:30],
            path=path, atypical=acc in ATYPICAL_HSP65))
    if species_order:
        refs.sort(key=lambda r: (species_order.index(r.species) if r.species in species_order else 99, r.species))
    return refs


_aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                           open_gap_score=-5, extend_gap_score=-2)
_aligner.end_deletion_score = 0  # the reference window may overhang the read freely


@dataclass
class RefHit:
    """The read laid over one reference: ref_cols[i] is the reference base at read position i
    ('-' = deletion in the reference or not covered, 'N' = ambiguous)."""
    ref: Reference
    ref_cols: str
    contig: str
    ref_start: int
    strand: str
    seeds: int


def _seed_clusters(ref: Reference, read: str) -> list[tuple[int, str, int, int]]:
    """(seed count, contig, diagonal, strand +1/-1) of the best-supported diagonals, best first."""
    out = []
    for strand, q in ((1, read), (-1, revcomp(read))):
        qc = _kmer_codes(q)
        valid = np.flatnonzero(qc != _SENTINEL)
        if not len(valid):
            continue
        uniq, idx, counts = np.unique(qc[valid], return_index=True, return_counts=True)
        keep = counts == 1
        qcodes, qoff = uniq[keep], valid[idx[keep]]
        gc = ref._codes
        pos = np.searchsorted(qcodes, gc)
        hit = np.flatnonzero(qcodes[np.minimum(pos, len(qcodes) - 1)] == gc)
        if len(hit) < MIN_SEEDS:
            continue
        ci = np.searchsorted(ref._starts, hit, side="right") - 1
        rel = hit - ref._starts[ci]  # position within the contig
        offs = qoff[pos[hit]]
        key = ci.astype(np.int64) * (1 << 32) + (rel - offs + (1 << 31))  # (contig, diagonal), sortable
        order = np.argsort(key)
        key = key[order]
        cuts = np.flatnonzero(np.diff(key) > 20) + 1
        for grp in np.split(key, cuts):
            if len(grp) >= MIN_SEEDS:
                med = int(np.median(grp))
                out.append((len(grp), ref._ids[med >> 32], (med & 0xFFFFFFFF) - (1 << 31), strand))
    return sorted(out, reverse=True)[:3]


def align_to_reference(ref: Reference, read: str) -> RefHit | None:
    """Best end-to-end placement of the read on the reference, None if no seeds."""
    ref.load()
    best = None
    for seeds, cid, diag, strand in _seed_clusters(ref, read):
        contig = ref.contigs[cid]
        lo, hi = max(0, diag - PAD), min(len(contig), diag + len(read) + PAD)
        q = read if strand == 1 else revcomp(read)
        aln = _aligner.align(contig[lo:hi], q)[0]
        if best is None or aln.score > best[0]:
            best = (aln.score, aln, cid, lo, strand, seeds)
    if best is None:
        return None
    _, aln, cid, lo, strand, seeds = best
    cols = ["-"] * len(read)  # in orientation of q
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(q1 - q0):
            cols[q0 + k] = ref.contigs[cid][lo + t0 + k]
    ref_cols = "".join(cols)
    ref_cols = "".join(c if c in "ACGT-" else "N" for c in ref_cols)
    if strand == -1:  # report in the orientation of the read as given
        ref_cols = ref_cols[::-1].translate(_COMP)
    return RefHit(ref, ref_cols, cid, lo, "+" if strand == 1 else "-", seeds)


def _diff(a: str, b: str) -> tuple[int, int]:
    """(differences, compared columns), over columns where both have an ACGT base."""
    m = d = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            m += 1
            d += x != y
    return d, m


@dataclass
class ReadResult:
    name: str
    locus: str | None
    read: str
    hits: list[RefHit]
    diffs: list[int]     # per hit: differences vs the read
    compared: list[int]  # per hit: columns compared
    min_identity: float = 0.98
    min_margin: int = 2

    @property
    def order(self) -> list[int]:
        return sorted(range(len(self.hits)), key=lambda i: (self.diffs[i] / max(self.compared[i], 1), -self.compared[i]))

    @property
    def best(self) -> int | None:
        return self.order[0] if self.hits else None

    def identity(self, i: int) -> float:
        return 1 - self.diffs[i] / self.compared[i] if self.compared[i] else 0.0

    @property
    def status(self) -> str:
        """ok: >= min_identity and clearly closest; ambiguous: >= min_identity but second
        reference within min_margin differences; divergent: closest is < min_identity."""
        if not self.hits:
            return "no_hit"
        o = self.order
        if self.identity(o[0]) < self.min_identity:
            return "divergent"
        r = self.runner_up
        if r is not None and self.diffs[r] - self.diffs[o[0]] < self.min_margin:
            return "ambiguous"
        return "ok"

    @property
    def closest_species(self) -> str:
        best = self.best
        return self.hits[best].ref.species if best is not None else "NA"

    @property
    def runner_up(self) -> int | None:
        """Best hit of another species than the closest (several genomes per species are fine)."""
        o = self.order
        return next((i for i in o[1:] if self.hits[i].ref.species != self.hits[o[0]].ref.species), None)

    @property
    def atypical(self) -> bool:
        """The closest reference is a kansasii genome with a persicum-like hsp65."""
        b = self.best
        return b is not None and self.hits[b].ref.atypical


def identify_read(refs: list[Reference], name: str, locus: str | None, seq: str,
                  min_identity: float = 0.98, min_margin: int = 2) -> ReadResult:
    hits = [h for h in (align_to_reference(r, seq) for r in refs) if h is not None]
    dc = [_diff(seq, h.ref_cols) for h in hits]
    return ReadResult(name, locus, seq, hits, [d for d, _ in dc], [c for _, c in dc], min_identity, min_margin)


