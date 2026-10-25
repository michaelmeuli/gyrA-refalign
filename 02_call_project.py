#!/usr/bin/env python
"""gyrA species call of every project isolate (Mkan329-*), from its own Unicycler assembly.

Compared with the WGS species (GTDB-Tk), the Sanger-based call and the hsp65 reference-alignment call.
Project isolates only: writes OUTPUT_PROJECT/isolates.tsv and SUMMARY.txt.
"""
from __future__ import annotations

import os

import pandas as pd

from gyrA_refalign import ASSEMBLIES, HSP65_ISOLATES, OUTPUT_PROJECT, REFERENCE_STRAINS, SCREENING_RESULTS
from gyrA_refalign.gyra import call_many


def main() -> None:
    OUTPUT_PROJECT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(SCREENING_RESULTS, dtype=str).fillna("")
    s = s[~s["PROBENNUMMER"].isin(REFERENCE_STRAINS) & (s["PROBENNUMMER"] != "")]
    items = [(pn, ASSEMBLIES / pn / "1_unicycler" / f"{pn}.fasta") for pn in s["PROBENNUMMER"]]
    items = [(pn, p) for pn, p in items if p.exists()]
    calls = call_many(items, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    df = pd.DataFrame.from_dict(calls, orient="index").rename_axis("probennummer").reset_index()
    meta = s.rename(columns={"PROBENNUMMER": "probennummer", "species": "wgs_species"})
    meta["contaminated"] = meta["contamination_flag"] != ""
    keep = ["probennummer", "wgs_species", "gtdb_ani", "contaminated", "species_sanger_mlsa", "species_sanger_refalign"]
    df = meta[keep].merge(df, on="probennummer", how="inner")
    if HSP65_ISOLATES.exists():
        h = pd.read_csv(HSP65_ISOLATES, sep="\t", dtype=str).fillna("")
        df = df.merge(h[["probennummer", "called_species", "verdict"]].rename(
            columns={"called_species": "hsp65_species", "verdict": "hsp65_verdict"}), on="probennummer", how="left")
    df["gyrA_matches_wgs"] = df["closest_species"] == df["wgs_species"]
    df.to_csv(OUTPUT_PROJECT / "isolates.tsv", sep="\t", index=False)

    clean = df[~df["contaminated"]]
    lines = [f"project isolates (reference strains 126-133 excluded): {len(df)}; uncontaminated: {len(clean)}",
             "gyrA status (all): " + str(df["status"].value_counts().to_dict()),
             "WGS species (uncontaminated): " + str(clean["wgs_species"].value_counts().to_dict()),
             "gyrA closest species (uncontaminated, status ok): "
             + str(clean[clean["status"] == "ok"]["closest_species"].value_counts().to_dict()),
             f"gyrA call = WGS species (uncontaminated, status ok): {int(clean[clean['status'] == 'ok']['gyrA_matches_wgs'].sum())}"
             f"/{int((clean['status'] == 'ok').sum())}",
             f"persicum-like gyrA (all): {int(df['persicum_like'].sum())}"]
    for _, w in df[df["persicum_like"] | (~df["gyrA_matches_wgs"] & (df["status"] == "ok"))].iterrows():
        lines.append(f"  {w['probennummer']}: WGS {w['wgs_species']}, gyrA {w['closest_species']} "
                     f"({w['identity']}, margin {w['margin']}), hsp65 {w.get('hsp65_species', '')}, contaminated {w['contaminated']}")
    non_ok = df[df["status"] != "ok"]
    for _, w in non_ok.iterrows():
        lines.append(f"  {w['status']}: {w['probennummer']} WGS {w['wgs_species']} gyrA {w['closest_species']} {w['identity']}")
    kans = clean[clean["wgs_species"] == "kansasii"]
    lines.append(f"kansasii gyrA haplotypes (uncontaminated): {kans['haplotype'].nunique()} distinct in {len(kans)}; "
                 f"exact match to a public genome: {int((kans['exact_match'] != '').sum())}")
    lines.append("hsp65 persicum-like but gyrA kansasii: " + (", ".join(
        df[(df.get("hsp65_species", "") == "persicum") & (df["closest_species"] == "kansasii")]["probennummer"]) or "none"))
    (OUTPUT_PROJECT / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
