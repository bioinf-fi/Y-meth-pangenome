#!/usr/bin/env python3

"""
Computes sliding-window percent identity for pairwise alignments
produced by EMBOSS Stretcher.
"""

import os
import argparse
from Bio import SeqIO
import pandas as pd


def sliding_window(seqA, seqB, window, step, min_valid):
    """
    Slide a window of `window` columns across an aligned sequence pair.
    Gap columns are skipped; windows with fewer than `min_valid` ungapped
    positions get percent_identity = None.
    Returns a list of dicts (one per window position).
    """
    rows = []
    aln_len = len(seqA)

    for start in range(0, aln_len - window + 1, step):
        end = start + window
        subA = seqA[start:end]
        subB = seqB[start:end]

        matches = 0
        valid = 0

        for a, b in zip(subA, subB):
            if a == "-" or b == "-":
                continue
            valid += 1
            if a == b:
                matches += 1

        if valid >= min_valid:
            pid = matches / valid * 100
        else:
            pid = None

        rows.append({
            "window_start": start,
            "window_end": end,
            "valid_positions": valid,
            "matches": matches,
            "percent_identity": pid
        })

    return rows


def process_alignment(fasta, sample, palindrome, args):
    """
    Parse a Stretcher output file (FASTA format with 2 records),
    run the sliding window, and tag each row with sample/palindrome metadata.
    """
    records = list(SeqIO.parse(fasta, "fasta"))
    if len(records) != 2:
        return []

    seqA = str(records[0].seq)
    seqB = str(records[1].seq)

    # Alignment sequences must be the same length
    if len(seqA) != len(seqB):
        return []

    rows = sliding_window(seqA, seqB, args.window, args.step, args.min_valid)

    for r in rows:
        r["sample"] = sample
        r["palindrome"] = palindrome
        r["alignment_length"] = len(seqA)

    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Sliding-window percent identity from Stretcher alignments."
    )
    parser.add_argument("--input", required=True,
                        help="stretcher_summaries directory (sample/palindrome.stretcher.txt)")
    parser.add_argument("--output", required=True,
                        help="Output TSV file")
    parser.add_argument("-w", "--window", type=int, default=200,
                        help="Window size in alignment columns (default: 200)")
    parser.add_argument("-s", "--step", type=int, default=50,
                        help="Step size between windows (default: 50)")
    parser.add_argument("--min-valid", type=int, default=10,
                        help="Minimum ungapped positions to report identity (default: 10)")

    args = parser.parse_args()

    all_rows = []

    # Traverse: input/<sample>/<palindrome>.stretcher.txt
    for sample in sorted(os.listdir(args.input)):
        sample_dir = os.path.join(args.input, sample)
        if not os.path.isdir(sample_dir):
            continue

        for fname in sorted(os.listdir(sample_dir)):
            if not fname.endswith(".stretcher.txt"):
                continue

            palindrome = fname.split(".")[0]
            fasta = os.path.join(sample_dir, fname)

            rows = process_alignment(fasta, sample, palindrome, args)
            all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    df = df.sort_values(["sample", "palindrome", "window_start"])

    df.to_csv(args.output, sep="\t", index=False)
    print(f"[OK] Written {len(df)} windows to {args.output}")


if __name__ == "__main__":
    main()
