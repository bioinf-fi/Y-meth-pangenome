#!/usr/bin/env python3

"""
Reads raw BED annotation files for chrY and classifies each region
(palindrome arm, spacer, inverted repeat, amplicon, …).
Rows belonging to excluded categories (PAR, CEN, DYZ, …) are dropped.
Output: cleaned BED per sample in the output directory.
"""

import os
import argparse
import pandas as pd


def is_palindrome_arm(name):
    """Return True if name matches the pattern P<digits> (e.g. P1, P8)."""
    if not name.startswith("P"):
        return False
    tail = name[1:]
    return tail.isdigit()


def classify_region(name):
    """Assign a region_type label based on the annotation name."""

    if is_palindrome_arm(name):
        return "palindrome_arm"

    if "spacer" in name:
        return "palindrome_spacer"

    # Sub-regions inside palindromes that also carry AZF locus names
    if name.startswith("P") and "AZF" in name:
        return "palindrome_subregion"

    if name.startswith("IR") or "-IR" in name:
        return "inverted_repeat"

    # Named colour amplicons used in Y-chromosome literature
    if name in [
        "blue", "green", "red", "yellow", "gray",
        "teal1", "teal2", "blue-plus"
    ]:
        return "amplicon"

    if name.startswith("AMPL"):
        return "amplicon"

    return "other"


def main(input_dir, output_dir):

    # Derive path to assemblies directory (sibling of input_dir)
    assemblies_dir = os.path.join(
        os.path.dirname(input_dir),
        "assemblies_unzipped"
    )

    # Collect sample names that have a chrY assembly
    valid_samples = set()
    for f in os.listdir(assemblies_dir):
        if f.endswith("_chrY.fa"):
            valid_samples.add(f.replace("_chrY.fa", ""))

    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(input_dir):

        if not fname.endswith(".bed"):
            continue

        sample = fname.replace("_new.bed", "")
        if sample not in valid_samples:
            continue

        path = os.path.join(input_dir, fname)
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None
        )

        if df.shape[1] < 6:
            print(f"[SKIP] {sample} – unexpected BED format")
            continue

        # Keep only the first 6 BED columns
        df = df.iloc[:, :6]
        df.columns = ["chrom", "start", "end", "name", "score", "strand"]
        df["name"] = df["name"].astype(str)
        df["strand"] = df["strand"].astype(str)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        # Drop random/unplaced contigs
        df = df[~df["chrom"].str.contains("random", na=False)]

        # Classify each annotation row
        df["region_type"] = df["name"].apply(classify_region)

        # Prefixes to exclude entirely (heterochromatin, PAR, satellites, …)
        EXCLUDE_PREFIXES = [
            "ERR", "PAR", "HET", "DYZ", "SAT",
            "CEN", "periCEN", "OTHER", "UNASSIGNED"
        ]

        def keep_row(name, region_type):
            if region_type == "other":
                return False
            for p in EXCLUDE_PREFIXES:
                if name.startswith(p):
                    return False
            return True

        df = df[df.apply(lambda r: keep_row(r["name"], r["region_type"]), axis=1)]

        if df.empty:
            print(f"[SKIP] {sample} – no valid regions after filtering")
            continue

        df["sample"] = sample
        df = df[[
            "chrom", "start", "end",
            "sample", "name", "region_type", "strand"
        ]]

        out = os.path.join(output_dir, f"{sample}.task02.cleaned.bed")
        df.to_csv(out, sep="\t", index=False)

        print(f"[OK] {sample}: {len(df)} regions → {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Clean and classify chrY BED annotations."
    )
    parser.add_argument("--input", required=True, help="Directory with raw BED annotations")
    parser.add_argument("--output", required=True, help="Output directory for cleaned BEDs")
    args = parser.parse_args()

    main(args.input, args.output)
