#!/usr/bin/env python
"""Read-level hsp65 amplicon of the public samples whose gyrA call disagrees with GTDB (05).

04 found the gyrA alleles of these samples to be real at read level, but only one hsp65 site separates kansasii
from persicum once the Korean-type genomes are left out, so 04 cannot say what hsp65 does. Here the raw reads are
mapped to a kansasii hsp65 amplicon (441 bp, TB-11/TB-12w), a consensus is called at every position and compared
with the 72 per-genome amplicons of the routine (HSP65_AMPLICONS). Per sample: differences to the best typical
kansasii, Korean-type, persicum and other-species amplicon, and the species call with its margin.

Standard library only (bwa and samtools from the snippy env). Writes OUTPUT_PUBLIC_READS/{hsp65_amplicon_check.tsv,
hsp65_consensus.fasta, HSP65_AMPLICON_SUMMARY.txt}. Public only.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the gyrA_refalign package

from gyrA_refalign import OUTPUT_PUBLIC, OUTPUT_PUBLIC_READS  # noqa: E402
from gyrA_refalign.kansasii_paths import ATYPICAL_HSP65_GENOMES, HSP65_AMPLICONS, PUBLIC_READS  # noqa: E402
from gyrA_refalign.tools import map_reads, read_fasta, read_tsv, run  # noqa: E402

MAP_REF = "GCF_000157895.3"  # M. kansasii ATCC 12478, the kansasii reference strain
MIN_BASEQ = 20
MIN_MAPQ = 20
MIN_DEPTH = 10  # reads at a position; fewer is not called
CALL_FRAC = 0.90  # major base share for a called position; below it the position is ambiguous (possible mixture)
PERSICUM_REP = "GCF_002086675.1"  # M. persicum AFPC-000227


def pileup(sam: str, length: int) -> list[Counter[str]]:
    """Base counts per reference position from SAM text (primary alignments, MAPQ and base quality filtered)."""
    cols: list[Counter[str]] = [Counter() for _ in range(length)]
    for line in sam.splitlines():
        if not line or line[0] == "@":
            continue
        f = line.split("\t")
        flag, rpos, mapq, cigar, seq, qual = int(f[1]), int(f[3]) - 1, int(f[4]), f[5], f[9], f[10]
        if flag & 0x904 or mapq < MIN_MAPQ or cigar == "*":
            continue
        q, r = 0, rpos
        for n_s, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            n = int(n_s)
            if op in "M=X":
                for j in range(n):
                    if 0 <= r + j < length and ord(qual[q + j]) - 33 >= MIN_BASEQ:
                        cols[r + j][seq[q + j]] += 1
                q += n
                r += n
            elif op in "IS":
                q += n
            elif op in "DN":
                r += n
    return cols


def consensus(cols: list[Counter[str]]) -> tuple[str, list[int], list[int]]:
    """Consensus string ('N' where not called), positions with too little depth, ambiguous positions (1-based)."""
    seq: list[str] = []
    low: list[int] = []
    amb: list[int] = []
    for i, c in enumerate(cols):
        tot = sum(c.values())
        if tot < MIN_DEPTH:
            seq.append("N")
            low.append(i + 1)
            continue
        base, n = c.most_common(1)[0]
        if n / tot >= CALL_FRAC:
            seq.append(base)
        else:
            seq.append("N")
            amb.append(i + 1)
    return "".join(seq), low, amb


def diffs(cons: str, ref: str) -> tuple[int, int]:
    """Differences and compared columns, over the positions where the consensus is called."""
    cmp_cols = [(a, b) for a, b in zip(cons, ref) if a != "N"]
    return sum(a != b for a, b in cmp_cols), len(cmp_cols)


def best(d: dict[str, tuple[int, int]], accs: list[str]) -> tuple[str, str]:
    """('diffs/columns', accession) of the closest of the given references; '' if there is none."""
    cand = [(d[a][0], a) for a in accs if a in d]
    if not cand:
        return "", ""
    n, a = min(cand)
    return f"{n}/{d[a][1]}", a


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--read-check", type=Path, default=OUTPUT_PUBLIC_READS / "read_check.tsv", help="samples of 04")
    ap.add_argument("--samples", type=Path, default=OUTPUT_PUBLIC / "samples.tsv")
    ap.add_argument("--references", type=Path, default=OUTPUT_PUBLIC / "reference_genomes.tsv")
    ap.add_argument("--extra", nargs="*", default=["SRR28038296"], help="further public samples (default: the Korean-type kansasii SRR28038296)")
    ap.add_argument("--out", type=Path, default=OUTPUT_PUBLIC_READS)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    meta = {r["sample"]: r for r in read_tsv(args.samples)}
    roles: dict[str, str] = {}
    for r in read_tsv(args.read_check):
        roles.setdefault(r["sample"], r["role"])
    for s in args.extra:
        roles.setdefault(s, "extra")
    species = {r["accession"]: r["gtdb_species"] for r in read_tsv(args.references)}
    amplicons = {p.name.removesuffix("_hsp65.fasta"): next(iter(read_fasta(p).values())) for p in sorted(HSP65_AMPLICONS.glob("*_hsp65.fasta"))}
    ref_len = len(amplicons[MAP_REF])
    refs = {a: s for a, s in amplicons.items() if len(s) == ref_len and a in species}
    print(f"{len(refs)} reference amplicons of {ref_len} bp ({len(amplicons) - len(refs)} left out: other length)", flush=True)
    korean = [a for a in ATYPICAL_HSP65_GENOMES if a in refs]
    typical_kans = [a for a, sp in species.items() if sp == "kansasii" and a in refs and a not in korean]
    pers = [a for a, sp in species.items() if sp == "persicum" and a in refs]

    rows: list[dict[str, str]] = []
    cons_out: list[str] = []
    with tempfile.TemporaryDirectory(prefix="hsp65_reads_", dir=args.out) as td:
        ref_fa = Path(td) / "ref.fasta"
        ref_fa.write_text(f">{MAP_REF}\n{amplicons[MAP_REF]}\n")
        run(["bwa", "index", str(ref_fa)])
        for s, role in roles.items():
            proj = meta[s]["bioproject"]
            rd = PUBLIC_READS / proj / "reads"
            sam = map_reads(ref_fa, rd / f"{s}_R1.fastq.gz", rd / f"{s}_R2.fastq.gz")
            cols = pileup(sam, ref_len)
            cons, low, amb = consensus(cols)
            d = {a: diffs(cons, seq) for a, seq in refs.items()}
            per_sp: dict[str, tuple[int, str]] = {}
            for a, (n, _) in d.items():
                sp = species[a]
                if sp not in per_sp or n < per_sp[sp][0]:
                    per_sp[sp] = (n, a)
            ranked = sorted((n, sp, a) for sp, (n, a) in per_sp.items())
            closest_n, closest_sp, closest_acc = ranked[0]
            second = next((x for x in ranked[1:]), None)
            bk, bk_a = best(d, typical_kans)
            bko, bko_a = best(d, korean)
            bp, bp_a = best(d, pers)
            rows.append({
                "sample": s, "role": role, "gtdb_species": meta[s]["gtdb_species"], "assembly_gyrA_call": meta[s]["closest_species"],
                "called_positions": str(ref_len - len(low) - len(amb)), "low_depth_positions": str(len(low)),
                "ambiguous_positions": ",".join(map(str, amb)), "median_depth": str(sorted(sum(c.values()) for c in cols)[ref_len // 2]),
                "diffs_typical_kansasii": bk, "typical_kansasii_acc": bk_a, "diffs_korean_type": bko, "korean_acc": bko_a,
                "diffs_persicum": bp, "persicum_acc": bp_a,
                "diffs_kansasii_ref": f"{d[MAP_REF][0]}/{d[MAP_REF][1]}", "diffs_persicum_ref": f"{d[PERSICUM_REP][0]}/{d[PERSICUM_REP][1]}",
                "closest_species_72": closest_sp, "closest_acc_72": closest_acc, "closest_diffs_72": str(closest_n),
                "second_species_72": second[1] if second else "", "margin_72": str(second[0] - closest_n) if second else "",
            })
            cons_out.append(f">{s} {role} hsp65 consensus (N = not called)\n{cons}\n")
            print(f"{s} ({role}) GTDB {meta[s]['gtdb_species']}: closest {closest_sp} {closest_acc} {closest_n} diffs, margin "
                  f"{rows[-1]['margin_72']}; typical kansasii {bk}, Korean {bko}, persicum {bp}; ambiguous {len(amb)}, low depth {len(low)}", flush=True)

    with (args.out / "hsp65_amplicon_check.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (args.out / "hsp65_consensus.fasta").write_text("".join(cons_out), encoding="utf-8")
    lines = [f"hsp65 amplicon ({ref_len} bp) from reads, consensus: base quality >= {MIN_BASEQ}, MAPQ >= {MIN_MAPQ}, depth >= {MIN_DEPTH}, "
             f"called if major base >= {CALL_FRAC}; {len(refs)} reference amplicons; Korean-type = {', '.join(korean)}",
             "diffs = differences/compared columns to the closest reference of that group",
             "sample\trole\tGTDB\ttypical kansasii\tKorean-type\tpersicum\tclosest species (72), margin"]
    for r in rows:
        lines.append(f"{r['sample']}\t{r['role']}\t{r['gtdb_species']}\t{r['diffs_typical_kansasii']}\t{r['diffs_korean_type']}\t"
                     f"{r['diffs_persicum']}\t{r['closest_species_72']} {r['closest_acc_72']}, margin {r['margin_72']}")
    (args.out / "HSP65_AMPLICON_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
