#!/usr/bin/env python
"""Read-level check of the public samples whose gyrA call disagrees with the GTDB species (04).

For each disagreeing sample (and two concordant controls) the raw reads are mapped to a kansasii reference of
each locus (gyrA, hsp65, rpoB, recA, tuf) and the base of every read is read out at the sites that are fixed
between kansasii and persicum in the public reference genomes. This separates three explanations:
  mosaic        reads carry the other species' allele at gyrA only, the other loci agree with GTDB
  mixed         both species' alleles are present in the reads, at several loci
  mislabelled   gyrA and the other loci all carry the other species (GTDB label and reads disagree)

Standard library only (run in the snippy env, which provides bwa, samtools and mafft).
Reads: PUBLIC_READS/<project>/reads/<sample>_R{1,2}.fastq.gz. Writes OUTPUT_PUBLIC_READS/{read_check,diagnostic_sites,
read_haplotypes}.tsv and READ_CHECK_SUMMARY.txt. Public only: the project isolates are not touched.
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

from gyrA_refalign import OUTPUT_PUBLIC, OUTPUT_PUBLIC_READS, PUBLIC_MAX_CONTAMINATION  # noqa: E402
from gyrA_refalign.kansasii_paths import GYRA_REFERENCES, OUT_MLSA, PUBLIC_READS, public_assemblies  # noqa: E402
from gyrA_refalign.tools import checkm_contamination, read_fasta, read_tsv, run  # noqa: E402

LOCI = ["gyrA", "rpoB", "recA", "tuf"]  # hsp65 is done in 05 on the whole amplicon; one fixed site is too few
LOCUS_ALIGNMENTS = OUT_MLSA / "results" / "locus_discovery" / "alignments"  # <locus>.raw.fasta, headers species__accession
MIN_BASEQ = 20
MIN_MAPQ = 20
MIN_DEPTH = 10  # reads at a site; fewer is not called
PURE = 0.95  # share of informative bases for a species to call the locus pure
MIXED = 0.10  # minor share above this at a locus counts as a second species present


def locus_sequences(locus: str, atypical: set[str], gyra_refs: Path) -> dict[str, tuple[str, str]]:
    """accession -> (species, sequence) for kansasii and persicum genomes; Korean-type genomes left out.

    gyrA comes from this repo's own re-extracted references; the other loci from the mlsa locus discovery alignments
    (read as data, not imported)."""
    out: dict[str, tuple[str, str]] = {}
    if locus == "gyrA":
        species = {r["accession"]: r["gtdb_species"] for r in read_tsv(OUTPUT_PUBLIC / "reference_genomes.tsv")}
        for acc, sp in species.items():
            f = gyra_refs / f"{acc}_gyrA.fasta"
            if sp in ("kansasii", "persicum") and acc not in atypical and f.exists():
                out[acc] = (sp, next(iter(read_fasta(f).values())))
    else:
        for name, seq in read_fasta(LOCUS_ALIGNMENTS / f"{locus}.raw.fasta").items():
            sp, acc = name.split("__", 1)
            if sp in ("kansasii", "persicum") and acc not in atypical:
                out[acc] = (sp, seq.replace("-", ""))
    return out


def diagnostic_sites(seqs: dict[str, tuple[str, str]], ref_acc: str, tmp: Path) -> list[tuple[int, str, str]]:
    """Reference positions (0-based on the reference) where all kansasii and all persicum sequences are fixed and differ."""
    fa = tmp / "msa_in.fasta"
    fa.write_text("".join(f">{a}\n{s}\n" for a, (_, s) in seqs.items()))
    aln = run(["mafft", "--quiet", "--auto", str(fa)]).stdout
    msa: dict[str, list[str]] = {}
    cur = ""
    for line in aln.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            msa[cur] = []
        else:
            msa[cur].append(line.strip().upper())
    cols = {a: "".join(v) for a, v in msa.items()}
    ref_row = cols[ref_acc]
    sites: list[tuple[int, str, str]] = []
    pos = -1
    for i, rb in enumerate(ref_row):
        if rb != "-":
            pos += 1
        k = {cols[a][i] for a, (sp, _) in seqs.items() if sp == "kansasii"}
        p = {cols[a][i] for a, (sp, _) in seqs.items() if sp == "persicum"}
        if len(k) == 1 and len(p) == 1 and k != p and "-" not in k | p and rb != "-" and "N" not in k | p:
            sites.append((pos, next(iter(k)), next(iter(p))))
    return sites


def read_bases(sam: str, sites: dict[int, tuple[str, str]]) -> tuple[Counter[tuple[int, str]], Counter[str]]:
    """Per-site base counts and per-read species strings from SAM text (primary, MAPQ and base quality filtered)."""
    site_counts: Counter[tuple[int, str]] = Counter()
    read_hap: Counter[str] = Counter()
    order = sorted(sites)
    for line in sam.splitlines():
        if not line or line[0] == "@":
            continue
        f = line.split("\t")
        flag, rpos, mapq, cigar, seq, qual = int(f[1]), int(f[3]) - 1, int(f[4]), f[5], f[9], f[10]
        if flag & 0x904 or mapq < MIN_MAPQ or cigar == "*":
            continue
        base_at: dict[int, str] = {}
        q = 0
        r = rpos
        for n_s, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            n = int(n_s)
            if op in "M=X":
                for j in range(n):
                    if r + j in sites and ord(qual[q + j]) - 33 >= MIN_BASEQ:
                        base_at[r + j] = seq[q + j]
                q += n
                r += n
            elif op in "IS":
                q += n
            elif op in "DN":
                r += n
        letters = ""
        for p in order:
            if p in base_at:
                b = base_at[p]
                site_counts[(p, b)] += 1
                letters += "K" if b == sites[p][0] else "P" if b == sites[p][1] else "x"
        if letters:
            read_hap["".join(sorted(set(letters)))] += 1
    return site_counts, read_hap


def check_sample(sample: str, project: str, locus: str, seqs: dict[str, tuple[str, str]], tmp: Path
                 ) -> tuple[list[tuple[int, str, str]], Counter[tuple[int, str]], Counter[str]]:
    ref_acc = next(a for a, (sp, _) in seqs.items() if sp == "kansasii")
    ref_fa = tmp / f"{locus}.ref.fasta"
    ref_fa.write_text(f">{ref_acc}\n{seqs[ref_acc][1]}\n")
    if not (tmp / f"{locus}.ref.fasta.bwt").exists():
        run(["bwa", "index", str(ref_fa)])
    sites = diagnostic_sites(seqs, ref_acc, tmp)
    rd = PUBLIC_READS / project / "reads"
    sam = run(["bwa", "mem", "-t", "4", str(ref_fa), str(rd / f"{sample}_R1.fastq.gz"), str(rd / f"{sample}_R2.fastq.gz")]).stdout
    counts, haps = read_bases(sam, {p: (k, pe) for p, k, pe in sites})
    return sites, counts, haps


def classify(frac_k: float | None, frac_p: float | None) -> str:
    if frac_k is None or frac_p is None:
        return "no_data"
    if frac_k >= PURE:
        return "kansasii"
    if frac_p >= PURE:
        return "persicum"
    if min(frac_k, frac_p) >= MIXED:
        return "mixed"
    return "mostly_kansasii" if frac_k > frac_p else "mostly_persicum"


def verdict(gtdb: str, gyra: str, loci: dict[str, str]) -> str:
    other = [c for lc, c in loci.items() if lc != "gyrA"]
    other_gtdb = sum(c == gtdb for c in other)
    n_mixed = sum(c == "mixed" for c in loci.values())
    if "no_data" in loci.values() and gyra == "no_data":
        return "low_depth"
    if n_mixed >= 2:
        return "mixed"
    if gyra not in (gtdb, "no_data") and gyra not in ("mixed",) and other_gtdb >= 3:
        return "mosaic"
    if gyra not in (gtdb, "mixed", "no_data") and all(c in (gyra, "no_data") for c in other):
        return "mislabelled"
    if gyra == gtdb and other_gtdb == len([c for c in other if c != "no_data"]):
        return "concordant"
    return "unclear"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", type=Path, default=OUTPUT_PUBLIC / "samples.tsv", help="gyrA call table from 03")
    ap.add_argument("--references", type=Path, default=OUTPUT_PUBLIC / "reference_genomes.tsv")
    ap.add_argument("--gyra-references", type=Path, default=GYRA_REFERENCES)
    ap.add_argument("--controls", type=int, default=1, help="concordant controls per species (default 1)")
    ap.add_argument("--max-contamination", type=float, default=PUBLIC_MAX_CONTAMINATION, help="CheckM contamination (%%) above which a sample is excluded")
    ap.add_argument("--out", type=Path, default=OUTPUT_PUBLIC_READS)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows = read_tsv(args.samples)
    refs = read_tsv(args.references)
    atypical = {r["accession"] for r in refs if r.get("atypical") == "True" or r["accession"] in {r2["closest_accession"] for r2 in rows if r2.get("atypical") == "True"}}
    excluded: list[str] = []

    def clean(r: dict[str, str], record: bool) -> bool:
        """Contamination filter; only the disagreeing samples are recorded as excluded (control candidates are just skipped)."""
        c = checkm_contamination(public_assemblies(r["bioproject"]), r["sample"])
        if c is None or c > args.max_contamination:
            if record:
                excluded.append(f"{r['sample']} (CheckM contamination {'n/a' if c is None else f'{c:.2f}'} %)")
            return False
        return True

    targets = [r for r in rows if r["status"] == "ok" and r["gyrA_matches_gtdb"] == "False" and clean(r, True)]
    # controls: concordant, exact reference haplotype, same project as the targets, one per species
    projects = {r["bioproject"] for r in targets}
    controls: list[dict[str, str]] = []
    for sp in ("kansasii", "persicum"):
        cand = [r for r in rows if r["bioproject"] in projects and r["gtdb_species"] == sp and r["status"] == "ok"
                and r["gyrA_matches_gtdb"] == "True" and r["exact_match"] != "" and r["atypical"] == "False"]
        controls += [r for r in cand if clean(r, False)][: args.controls]
    print(f"{len(targets)} disagreeing samples, {len(controls)} controls; excluded as contaminated: {excluded or 'none'}", flush=True)

    seqs = {lc: locus_sequences(lc, atypical, args.gyra_references) for lc in LOCI}
    out_rows: list[dict[str, str]] = []
    site_rows: list[dict[str, str]] = []
    hap_rows: list[dict[str, str]] = []
    verdicts: dict[str, str] = {}
    with tempfile.TemporaryDirectory(prefix="gyra_reads_", dir=args.out) as td:
        tmp_root = Path(td)
        for role, group in (("disagree", targets), ("control", controls)):
            for r in group:
                s, proj, gtdb, gcall = r["sample"], r["bioproject"], r["gtdb_species"], r["closest_species"]
                calls: dict[str, str] = {}
                for lc in LOCI:
                    tmp = tmp_root / lc
                    tmp.mkdir(exist_ok=True)
                    sites, counts, haps = check_sample(s, proj, lc, seqs[lc], tmp)
                    kt = pt = ot = 0
                    depths: list[int] = []
                    for p, kb, pb in sites:
                        k, pe = counts[(p, kb)], counts[(p, pb)]
                        tot = sum(v for (pp, _), v in counts.items() if pp == p)
                        depths.append(tot)
                        site_rows.append({"sample": s, "locus": lc, "ref_pos_1based": str(p + 1), "kansasii_base": kb,
                                          "persicum_base": pb, "n_kansasii": str(k), "n_persicum": str(pe), "n_other": str(tot - k - pe)})
                        if tot >= MIN_DEPTH:
                            kt, pt, ot = kt + k, pt + pe, ot + tot - k - pe
                    inf = kt + pt
                    fk = kt / inf if inf else None
                    fp = pt / inf if inf else None
                    good = [d for d in depths if d >= MIN_DEPTH]
                    calls[lc] = classify(fk, fp) if good else "no_data"
                    out_rows.append({"sample": s, "role": role, "bioproject": proj, "gtdb_species": gtdb, "gyrA_call_assembly": gcall,
                                     "locus": lc, "diagnostic_sites": str(len(sites)), "sites_with_depth": str(len(good)),
                                     "median_depth": str(sorted(depths)[len(depths) // 2] if depths else 0),
                                     "frac_kansasii": "" if fk is None else f"{fk:.4f}", "frac_persicum": "" if fp is None else f"{fp:.4f}",
                                     "other_bases": str(ot), "read_level_call": calls[lc]})
                    for hp, n in sorted(haps.items()):
                        hap_rows.append({"sample": s, "locus": lc, "read_pattern": hp, "n_reads": str(n)})
                verdicts[s] = verdict(gtdb, calls["gyrA"], calls)
                for o in out_rows:
                    if o["sample"] == s:
                        o["verdict"] = verdicts[s]
                print(f"{s} ({role}) GTDB {gtdb}, assembly gyrA {gcall}: {calls} -> {verdicts[s]}", flush=True)

    def write(name: str, data: list[dict[str, str]]) -> None:
        with (args.out / name).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t")
            w.writeheader()
            w.writerows(data)

    write("read_check.tsv", out_rows)
    write("diagnostic_sites.tsv", site_rows)
    write("read_haplotypes.tsv", hap_rows)
    lines = [f"thresholds: base quality >= {MIN_BASEQ}, MAPQ >= {MIN_MAPQ}, site depth >= {MIN_DEPTH}, pure >= {PURE}, mixed minor >= {MIXED}",
             f"excluded (CheckM contamination > {args.max_contamination} % or missing): {', '.join(excluded) or 'none'}",
             "verdict per sample (disagreeing samples first, then controls):"]
    for r in targets + controls:
        lines.append(f"  {r['sample']}: GTDB {r['gtdb_species']}, assembly gyrA {r['closest_species']} -> {verdicts[r['sample']]}")
    (args.out / "READ_CHECK_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
