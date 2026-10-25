#!/usr/bin/env python
"""gyrA species call of the public WGS samples (runs/public_PRJNA*) that GTDB places in the kansasii complex.

Public samples only (project isolates are in 02_call_project.py); the public reference genomes are in
OUTPUT_PUBLIC/reference_genomes.tsv from 01. Writes OUTPUT_PUBLIC/samples.tsv and SAMPLES_SUMMARY.txt.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the gyrA_refalign package

from gyrA_refalign import OUTPUT_PUBLIC  # noqa: E402
from gyrA_refalign.gyra import call_many  # noqa: E402
from gyrA_refalign.kansasii_paths import GYRA_REFERENCES, PUBLIC_PROJECTS, SPECIES, public_assemblies  # noqa: E402


def gtdb_species(summary: str) -> str:
    row = pd.read_csv(summary, sep="\t", dtype=str).iloc[0]
    return str(row["classification"]).split(";")[-1].replace("s__Mycobacterium ", "").replace("s__", "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("projects", nargs="*", default=PUBLIC_PROJECTS, help="public BioProjects (default: the four of the project)")
    ap.add_argument("--gyra-references", type=Path, default=GYRA_REFERENCES, help="<accession>_gyrA.fasta of the 72 genomes (from 01)")
    ap.add_argument("--out", type=Path, default=OUTPUT_PUBLIC)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    meta: dict[str, tuple[str, str]] = {}
    items: list[tuple[str, Path]] = []
    for run in args.projects:
        for d in sorted(public_assemblies(run).iterdir()):
            summ = list((d / "3_quality" / "GTDB").glob("*_gtdb_summary.tsv"))
            fasta = d / "1_unicycler" / f"{d.name}.fasta"
            if not summ or not fasta.exists():
                continue
            sp = gtdb_species(str(summ[0]))
            if sp in SPECIES:
                meta[d.name] = (run, sp)
                items.append((d.name, fasta))
    print(f"{len(items)} public samples in the complex", flush=True)
    calls = call_many(items, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")), args.gyra_references)
    df = pd.DataFrame.from_dict(calls, orient="index").rename_axis("sample").reset_index()
    df.insert(1, "bioproject", df["sample"].map(lambda x: meta[x][0]))
    df.insert(2, "gtdb_species", df["sample"].map(lambda x: meta[x][1]))
    df["gyrA_matches_gtdb"] = df["closest_species"] == df["gtdb_species"]
    df.sort_values(["bioproject", "sample"]).to_csv(args.out / "samples.tsv", sep="\t", index=False)

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
    (args.out / "SAMPLES_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
