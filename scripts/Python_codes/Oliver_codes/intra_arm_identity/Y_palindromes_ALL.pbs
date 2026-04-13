#!/bin/bash
#PBS -N Y_PAL_SINGLE
#PBS -q default
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=100gb
#PBS -l walltime=12:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT


# -------------------------------------------------------------
# INPUT CHECK
# Expects PAL_FASTA passed via: qsub -v PAL_FASTA=...
# -------------------------------------------------------------
if [ -z "$PAL_FASTA" ]; then
  echo "ERROR: provide input FASTA via:"
  echo "qsub -v PAL_FASTA=/path/to/sample_palindromes.fa Y_palindromes_ALL.pbs"
  exit 1
fi


# -------------------------------------------------------------
# PATHS – adjust BASE to your project root
# -------------------------------------------------------------
BASE="${PROJECT_BASE:-/path/to/project}"
OUTDIR="$BASE/03_palindrome_identity_manual"
SUMMARY_BASE="$OUTDIR/stretcher_summaries"

mkdir -p "$OUTDIR"
mkdir -p "$SUMMARY_BASE"


# -------------------------------------------------------------
# LOAD MODULES
# -------------------------------------------------------------
module add bedtools
module add mambaforge
module add seqtk
mamba activate emboss          # activate EMBOSS conda environment


# -------------------------------------------------------------
# COPY INPUT TO SCRATCH
# -------------------------------------------------------------
sample=$(basename "$PAL_FASTA" _palindromes.fa)
echo "[$(date)] Processing sample: $sample"

cp "$PAL_FASTA" "$SCRATCHDIR/palindromes.fa" || exit 2
cd "$SCRATCHDIR" || exit 3


# -------------------------------------------------------------
# OUTPUT TABLE HEADER
# -------------------------------------------------------------
echo -e "sample\tpalindrome\talignment_length\tpercent_identity" > results.tsv


# -------------------------------------------------------------
# MAIN LOOP – iterate over each unique palindrome in the FASTA
# For each palindrome: extract both arms, reverse-complement arm B,
# align with Stretcher (global Needleman-Wunsch), parse identity
# -------------------------------------------------------------
for p in $(grep "^>" palindromes.fa | cut -d: -f1 | sed 's/>//' | sort -u); do

  echo "[$(date)] $sample – Palindrome $p"

  # Extract all sequences belonging to this palindrome
  awk -v p=">$p::" '
  BEGIN {keep=0}
  /^>/ {keep = (index($0,p)==1)}
  keep {print}
  ' palindromes.fa > tmp.fa

  # Skip if we don't have exactly 2 arms
  if [ $(grep "^>" tmp.fa | wc -l) -ne 2 ]; then
    echo "WARNING: $sample $p does not have exactly 2 arms, skipping" >&2
    rm -f tmp.fa
    continue
  fi

  # Split into individual arm FASTA files
  awk '
  BEGIN {i=0}
  /^>/ {
    i++
    f = (i==1 ? "armA.fa" : "armB.fa")
    print > f
    next
  }
  { print > f }
  ' tmp.fa

  # Reverse-complement arm B (palindromes are inverted repeats)
  seqtk seq -r armB.fa > armB.rc.fa

  # Determine the shorter arm length for normalisation
  lenA=$(grep -v "^>" armA.fa | tr -d '\n' | wc -c)
  lenB=$(grep -v "^>" armB.rc.fa | tr -d '\n' | wc -c)

  if [ "$lenA" -lt "$lenB" ]; then
    shorter="$lenA"
  else
    shorter="$lenB"
  fi

  # Global pairwise alignment with EMBOSS Stretcher
  stretcher \
    -asequence armA.fa \
    -bsequence armB.rc.fa \
    -outfile aln.report \
    -auto

  # Save full alignment report for downstream sliding-window analysis
  mkdir -p "$SUMMARY_BASE/$sample"
  cp aln.report "$SUMMARY_BASE/$sample/${p}.stretcher.txt"

  # Parse identity from the report and normalise by the shorter arm length
  awk -v s="$sample" -v p="$p" -v short="$shorter" '
/^# Identity:/ {
  split($3, a, "/")
  matches = a[1]
}
END {
  if (short > 0) {
    norm_pid = (matches / short) * 100
    printf "%s\t%s\t%d\t%.2f\n", s, p, short, norm_pid
  }
}
' aln.report >> results.tsv

  rm -f armA.fa armB.fa armB.rc.fa aln.report tmp.fa
done

# -------------------------------------------------------------
# COPY RESULT BACK FROM SCRATCH
# -------------------------------------------------------------
cp results.tsv "$OUTDIR/${sample}.tsv" || export CLEAN_SCRATCH=false

echo "[$(date)] DONE sample $sample"
