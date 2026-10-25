#!/usr/bin/env python
"""gyrA references from the 72 public GTDB complex genomes, and a leave-one-out check of them.

1. The probe is the annotated gyrA CDS (GFF `gene=gyrA`) of M. persicum AFPC-000227, written to GYRA_REFERENCES.
2. The probe is laid over every genome (data/gtdb_genomes/.../mlsa-kansasii/<species>/<accession>.fna); the covered
   columns are that genome's gyrA, written as <accession>_gyrA.fasta with a `>acc Mycobacterium <species>` header.
3. Each reference is called against all the others (leave-one-out). Public genomes only: written to OUTPUT_PUBLIC.
"""
from __future__ import annotations

import sys

import pandas as pd
from Bio import SeqIO

from gyrA_refalign import GTDB_GENOMES, GYRA_REFERENCES, OUTPUT_PUBLIC, PROBE_ACCESSION, PROBE_FASTA, SPECIES
from gyrA_refalign.gyra import MIN_COVERAGE, call_gyra, extract_gyra, load_gyra_references


def probe_from_gff() -> str:
    species = next(d.name for d in GTDB_GENOMES.iterdir() if (d / f"{PROBE_ACCESSION}.fna").exists())
    fna, gff = GTDB_GENOMES / species / f"{PROBE_ACCESSION}.fna", GTDB_GENOMES / species / f"{PROBE_ACCESSION}.gff"
    contigs = {r.id: r.seq for r in SeqIO.parse(fna, "fasta")}
    best = None
    for line in open(gff, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or "gene=gyrA;" not in f[8] + ";":
            continue
        seq = contigs[f[0]][int(f[3]) - 1:int(f[4])]
        seq = seq.reverse_complement() if f[6] == "-" else seq
        if best is None or len(seq) > len(best):
            best = seq
    if best is None:
        sys.exit(f"no gyrA in {gff}")
    return str(best).upper()


def main() -> None:
    GYRA_REFERENCES.mkdir(parents=True, exist_ok=True)
    OUTPUT_PUBLIC.mkdir(parents=True, exist_ok=True)
    probe = probe_from_gff()
    print(f"probe {PROBE_ACCESSION}: {len(probe)} bp, starts {probe[:3]}, ends {probe[-3:]}")
    PROBE_FASTA.write_text(f">{PROBE_ACCESSION}_gyrA\n{probe}\n", encoding="utf-8")

    rows: list[dict[str, object]] = []
    for sdir in sorted(GTDB_GENOMES.iterdir()):
        for fna in sorted(sdir.glob("*.fna")):
            acc = fna.stem
            ext = extract_gyra(fna, probe)
            seq, cov = ext if ext else ("", 0.0)
            ok = cov >= MIN_COVERAGE
            if ok:
                (GYRA_REFERENCES / f"{acc}_gyrA.fasta").write_text(f">{acc} Mycobacterium {sdir.name}\n{seq}\n", encoding="utf-8")
            rows.append({"accession": acc, "gtdb_species": sdir.name, "gyrA_bp": len(seq), "gyrA_coverage": round(cov, 4),
                         "written": ok})
    df = pd.DataFrame(rows)
    print(f"{int(df['written'].sum())} of {len(df)} genomes give a gyrA reference")

    refs = load_gyra_references()
    loo: dict[str, dict[str, object]] = {}
    for r in refs:
        others = [o for o in refs if o.accession != r.accession]
        loo[r.accession] = call_gyra(r.accession, r.path, probe, others)
    res = pd.DataFrame.from_dict(loo, orient="index").rename_axis("accession").reset_index()
    out = df.merge(res, on="accession", how="left", suffixes=("", "_loo"))
    out["loo_correct"] = out["closest_species"] == out["gtdb_species"]
    out.to_csv(OUTPUT_PUBLIC / "reference_genomes.tsv", sep="\t", index=False)

    called = out[out["written"]]
    lines = [f"probe: {PROBE_ACCESSION} gyrA, {len(probe)} bp",
             f"genomes: {len(df)}; with gyrA (>= {MIN_COVERAGE:.0%} of probe columns): {len(called)}",
             "gyrA length per species: " + "; ".join(f"{s} {sorted(set(g['gyrA_bp']))}" for s, g in called.groupby("gtdb_species")),
             f"leave-one-out: {int(called['loo_correct'].sum())}/{len(called)} closest species = GTDB species",
             "leave-one-out status: " + str(called["status"].value_counts().to_dict())]
    for _, w in called[~called["loo_correct"]].iterrows():
        lines.append(f"  WRONG {w['accession']} ({w['gtdb_species']}) -> {w['closest_species']} ({w['status']})")
    for _, w in out[~out["written"]].iterrows():
        lines.append(f"  no gyrA: {w['accession']} ({w['gtdb_species']}), coverage {w['gyrA_coverage']}")
    kor = called[called["accession"].isin({"GCF_002705785.1", "GCF_002705825.1", "GCF_002705865.1"})]
    lines.append("Korean kansasii genomes (persicum-like hsp65) leave-one-out: " + ", ".join(kor["closest_species"]))
    (OUTPUT_PUBLIC / "REFERENCES_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
