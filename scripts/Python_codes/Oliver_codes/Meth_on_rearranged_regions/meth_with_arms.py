#!/usr/bin/env python3

"""
Intersects per-sample CpG methylation bedgraphs with cleaned chrY
annotation BEDs. Assigns arm labels (A/B) to palindrome arms and
propagates them to overlapping same-strand regions. Outputs a TSV
with mean and variance of methylation per annotated region.
"""

import os
import argparse
import pandas as pd


def intersect_sample(meth_path, bed_path, sample):
    """
    Load methylation bedgraph and annotation BED for one sample,
    assign armIDs to palindrome arms (A = proximal, B = distal by
    genomic coordinate), propagate armIDs to overlapping regions
    on the same strand, then intersect with methylation data.
    Returns a list of result dicts.
    """
    meth = pd.read_csv(
        meth_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "methyl"]
    )

    # Load cleaned annotation BED (has header from clean_y_annotations.py)
    ann = pd.read_csv(bed_path, sep="\t", header=0)
    ann["armID"] = "-"

    # --------------------------------------------------
    # 1) Assign armID to palindrome arms (A = first by start, B = second)
    # --------------------------------------------------
    pal_arms = ann[ann["region_type"] == "palindrome_arm"]

    for pname, grp in pal_arms.groupby("name"):
        grp = grp.sort_values("start")

        if len(grp) == 1:
            ann.loc[grp.index, "armID"] = f"{pname}A"
        else:
            ann.loc[grp.index[0], "armID"] = f"{pname}A"
            ann.loc[grp.index[1], "armID"] = f"{pname}B"

    # --------------------------------------------------
    # 2) Propagate armID to regions that overlap a palindrome arm
    #    on the SAME chromosome AND same strand
    # --------------------------------------------------
    for _, arm in ann[ann["region_type"] == "palindrome_arm"].iterrows():

        mask = (
            (ann["chrom"] == arm["chrom"]) &
            (ann["strand"] == arm["strand"]) &
            (ann["end"] > arm["start"]) &
            (ann["start"] < arm["end"]) &
            (ann["region_type"] != "palindrome_arm")   # don't overwrite arms themselves
        )

        ann.loc[mask, "armID"] = arm["armID"]

    # --------------------------------------------------
    # 3) Intersect each annotation region with methylation data
    # --------------------------------------------------
    results = []

    for _, r in ann.iterrows():

        # Find methylation positions that overlap this region
        sub = meth[
            (meth["end"] > r["start"]) &
            (meth["start"] < r["end"])
        ]

        if sub.empty:
            continue

        results.append({
            "sample": sample,
            "armID": r["armID"],
            "region_type": r["region_type"],
            "region_name": r["name"],
            "strand": r["strand"],
            "n_cpg": len(sub),
            "mean_methyl": sub["methyl"].mean(),
            "var_methyl": sub["methyl"].var()
        })

    return results


def main(ann_dir, meth_dir, output):

    all_results = []

    for fname in os.listdir(ann_dir):

        if not fname.endswith(".task02.cleaned.bed"):
            continue

        sample = fname.split(".")[0]

        # Expected methylation filename convention
        meth_file = os.path.join(meth_dir, f"{sample}_chrY_meth.bedgraph")

        if not os.path.exists(meth_file):
            print(f"[SKIP] No methylation file for {sample}")
            continue

        ann_path = os.path.join(ann_dir, fname)
        print(f"[RUN] {sample}")

        res = intersect_sample(meth_file, ann_path, sample)
        all_results.extend(res)

    df = pd.DataFrame(all_results)

    # Explicit column order for reproducibility
    df = df[[
        "sample", "armID", "region_type", "region_name",
        "strand", "n_cpg", "mean_methyl", "var_methyl",
    ]]

    df.to_csv(output, sep="\t", index=False)
    print(f"[OK] Written {len(df)} rows to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Intersect CpG methylation with chrY palindrome arm annotations."
    )
    parser.add_argument("--annotations", required=True, help="Directory with task02 cleaned BEDs")
    parser.add_argument("--methylation", required=True, help="Directory with methylation bedgraphs")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    main(args.annotations, args.methylation, args.output)
