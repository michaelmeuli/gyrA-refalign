#!/usr/bin/env python
"""gyrA species call of the public WGS samples (runs/public_PRJNA*) that GTDB places in the kansasii complex.

Public samples only (project isolates are in 02_call_project.py); the public reference genomes are in
OUTPUT_PUBLIC/reference_genomes.tsv from 01. Writes OUTPUT_PUBLIC/samples.tsv and SAMPLES_SUMMARY.txt.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from gyrA_refalign import OUTPUT_PUBLIC, PUBLIC_RUNS, RUNS, SPECIES
from gyrA_refalign.gyra import call_many


def gtdb_species(summary: str) -> str:
    row = pd.read_csv(summary, sep="\t", dtype=str).iloc[0]
    return str(row["classification"]).split(";")[-1].replace("s__Mycobacterium ", "").replace("s__", "")


def main() -> None:
    OUTPUT_PUBLIC.mkdir(parents=True, exist_ok=True)
    meta: dict[str, tuple[str, str]] = {}
    items: list[tuple[str, Path]] = []
    for run in PUBLIC_RUNS:
        for d in sorted((RUNS / run / "assembly" / "results").iterdir()):
            summ = list((d / "3_quality" / "GTDB").glob("*_gtdb_summary.tsv"))
            fasta = d / "1_unicycler" / f"{d.name}.fasta"
            if not summ or not fasta.exists():
                continue
            sp = gtdb_species(str(summ[0]))
            if sp in SPECIES:
                meta[d.name] = (run.removeprefix("public_"), sp)
                items.append((d.name, fasta))
    print(f"{len(items)} public samples in the complex", flush=True)
    calls = call_many(items, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    df = pd.DataFrame.from_dict(calls, orient="index").rename_axis("sample").reset_index()
    df.insert(1, "bioproject", df["sample"].map(lambda x: meta[x][0]))
    df.insert(2, "gtdb_species", df["sample"].map(lambda x: meta[x][1]))
    df["gyrA_matches_gtdb"] = df["closest_species"] == df["gtdb_species"]
    df.sort_values(["bioproject", "sample"]).to_csv(OUTPUT_PUBLIC / "samples.tsv", sep="\t", index=False)

    ok = df[df["status"] == "ok"]
    lines = [f"public samples in the complex (GTDB): {len(df)}", "GTDB species: " + str(df["gtdb_species"].value_counts().to_dict()),
             "gyrA status: " + str(df["status"].value_counts().to_dict()),
             f"gyrA call = GTDB species (status ok): {int(ok['gyrA_matches_gtdb'].sum())}/{len(ok)}",
             f"persicum-like gyrA: {int(df['persicum_like'].sum())} (GTDB persicum: {int((df['gtdb_species'] == 'persicum').sum())})"]
    for _, w in df[~df["gyrA_matches_gtdb"] & df["status"].isin(["ok", "ambiguous"])].iterrows():
        lines.append(f"  {w['sample']} ({w['bioproject']}): GTDB {w['gtdb_species']}, gyrA {w['closest_species']} "
                     f"({w['identity']}, margin {w['margin']}, {w['status']})")
    for _, w in df[~df["status"].isin(["ok"])].iterrows():
        lines.append(f"  {w['status']}: {w['sample']} ({w['bioproject']}) GTDB {w['gtdb_species']} coverage {w['gyrA_coverage']} identity {w['identity']}")
    k = df[df["gtdb_species"] == "kansasii"]
    lines.append(f"kansasii (GTDB): {len(k)} samples, {k['haplotype'].nunique()} gyrA haplotypes; "
                 f"exact match to a reference genome: {int((k['exact_match'] != '').sum())}")
    (OUTPUT_PUBLIC / "SAMPLES_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
