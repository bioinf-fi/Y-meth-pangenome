#!/usr/bin/env python3

"""
From cleaned annotation BEDs (output of clean_y_annotations.py),
keep only palindrome arms where BOTH arms (exactly 2 records per
palindrome name) are present. Outputs a 6-column BED per sample.
"""

import os
import argparse
import pandas as pd


def main(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(input_dir):
        if not fname.endswith(".task02.cleaned.bed"):
            continue

        sample = fname.replace(".task02.cleaned.bed", "")
        path = os.path.join(input_dir, fname)

        df = pd.read_csv(path, sep="\t")
        if df.empty:
            print(f"[SKIP] {sample} – empty file")
            continue

        # Keep only palindrome arm rows
        df = df[df["region_type"] == "palindrome_arm"].copy()
        if df.empty:
            print(f"[SKIP] {sample} – no palindrome arms")
            continue

        df["name"] = df["name"].astype(str)

        # Retain only palindromes that have exactly 2 arms (complete pairs)
        counts = df.groupby("name").size()
        valid_palindromes = counts[counts == 2].index

        df = df[df["name"].isin(valid_palindromes)]
        if df.empty:
            print(f"[SKIP] {sample} – no complete palindromes")
            continue

        # Build standard 6-column BED (score set to 0)
        bed_out = df[["chrom", "start", "end", "name", "strand"]].copy()
        bed_out["score"] = 0
        bed_out = bed_out[["chrom", "start", "end", "name", "score", "strand"]]

        out_path = os.path.join(output_dir, f"{sample}.bed")
        bed_out.to_csv(out_path, sep="\t", index=False, header=False)
        print(f"[OK] {sample}: {len(valid_palindromes)} palindromes → {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract complete palindrome arm pairs from cleaned BED files."
    )
    parser.add_argument("--input", required=True,
                        help="Directory with cleaned annotation BED files")
    parser.add_argument("--output", required=True,
                        help="Output directory for palindromes with both arms")
    args = parser.parse_args()
    main(args.input, args.output)
