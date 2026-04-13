#!/bin/bash

# Shell commands used for per-palindrome sequence composition analysis.
# Tools required: seqrequester, seqkit, dustmasker, RepeatMasker, trf
# All paths are relative to the working directory unless noted.

# =============================================================
# 1  Nucleotide composition with seqrequester (one file per sample)
# =============================================================
for f in fastapalindromes/*.fa; do
    name=$(basename "$f" .fa)
    seqrequester/build/bin/seqrequester summarize "$f" > "composition_${name}.txt"
done

# Aggregate into a single TSV with A/C/G/T counts and GC content
echo -e "sample\tA\tC\tG\tT\tGC" > nucleotide_composition.tsv

for f in composition_*.txt; do
    sample=$(basename "$f" .txt | sed 's/composition_//')

    A=$(awk '$3=="A"{print $2}' "$f")
    C=$(awk '$3=="C"{print $2}' "$f")
    G=$(awk '$3=="G"{print $2}' "$f")
    T=$(awk '$3=="T"{print $2}' "$f")

    # GC line follows "--GC--" header in seqrequester output
    GC=$(grep -A1 "\-\-GC\-\-" "$f" | tail -n1 | awk '{print $1}')

    echo -e "$sample\t$A\t$C\t$G\t$T\t$GC" >> nucleotide_composition.tsv
done

# Per-arm GC content using seqkit (optional, faster alternative)
seqkit fx2tab -n -g fastapalindromes/*.fa > arm_composition.tsv


# =============================================================
# 2  Low-complexity regions with DustMasker
#    Outputs: masked FASTA, then summary TSV (total bp, masked bp, density)
# =============================================================
mkdir -p dustmasker_results

for f in fastapalindromes/*.fa; do
    base=$(basename "$f" .fa)

    dustmasker \
        -in "$f" \
        -out "dustmasker_results/${base}_masked.fa" \
        -outfmt fasta \
        -hard_masking
done

# Parse masked (N) bases from each masked FASTA
for f in dustmasker_results/*_masked.fa; do
    sample=$(basename "$f" _palindromes_masked.fa)

    awk -v sample="$sample" '
    /^>/ { if (seq) { process(seq) }; header=$0; seq=""; next }
    { seq=seq$0 }
    END  { process(seq) }

    function process(s,   total, masked, i) {
        total  = length(s)
        masked = 0
        for (i=1; i<=total; i++) {
            if (substr(s,i,1) == "N") masked++
        }
        print sample "\t" header "\t" total "\t" masked "\t" masked/total
    }
    ' "$f"
done > low_complexity.tsv

# Simplify: extract sample, palindrome name, and masked fraction
awk '
{
    gsub(">", "", $2)
    split($2, a, "::")
    palindrome = a[1]
    print $1 "\t" palindrome "\t" $5
}
' low_complexity.tsv > low_complexity_clean.tsv


# =============================================================
# 3  Repeat annotation with RepeatMasker (species: human)
# =============================================================
mkdir -p ANALYSIS-per-sample/repeatmasker_results

for f in fastapalindromes/*_palindromes.fa; do
    sample=$(basename "$f" _palindromes.fa)

    RepeatMasker \
        -species human \
        -pa 4 \
        -gff \
        -dir ANALYSIS-per-sample/repeatmasker_results \
        "$f"
done


# =============================================================
# 4  Homopolymer runs (default minimum length = 10 bp)
#    Output columns: sample, sequence_name, start, end, base, run_length
# =============================================================
HOMOPOLYMER_LENGTH=10
mkdir -p ANALYSIS-per-arm/homopolymer_runs

for f in fastapalindromes/*_palindromes.fa; do
    sample=$(basename "$f" _palindromes.fa)

    # Linearise multi-line FASTA, then scan for runs of identical bases
    awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' "$f" | \
    awk -v len="${HOMOPOLYMER_LENGTH}" -v sample="$sample" -F '' '
    /^>/{
        gsub("^>", "", $0)
        seq = $0
        next
    }
    {
        base  = $1
        start = 1
        end   = 1

        for (i=2; i<=NF+1; i++) {
            if (base == $(i)) {
                end = i
            } else {
                if (end - start >= len - 1) {
                    print sample "\t" seq "\t" start "\t" end "\t" base "\t" (end-start+1)
                }
                start = i
                end   = i
                base  = $(i)
            }
        }
    }' > "ANALYSIS-per-arm/homopolymer_runs/${sample}.tsv"
done


# =============================================================
# 5  Microsatellite (tandem repeat) detection with TRF
#    Parameters: 2 7 7 80 10 50 500 (standard for short TRs)
#    Output columns: sample, sequence_name, total_bp_in_TRs, tr_count
# =============================================================
mkdir -p ANALYSIS-per-arm/trf_results

for f in fastapalindromes/*_palindromes.fa; do
    sample=$(basename "$f" _palindromes.fa)

    trf "$f" 2 7 7 80 10 50 500 -d -h

    mv ./*.dat "ANALYSIS-per-arm/trf_results/${sample}.dat"
done

# Summarise TRF .dat files: total bp covered by repeats and repeat count per sequence
for f in ANALYSIS-per-arm/trf_results/*.dat; do
    sample=$(basename "$f" .dat)

    awk -v sample="$sample" '
    BEGIN { OFS="\t" }

    /^Sequence:/ { seq = $2 }

    /^[0-9]/ {
        start = $1
        end   = $2
        len   = end - start + 1
        count[seq]++
        bp[seq] += len
    }

    END {
        for (s in bp) {
            print sample, s, bp[s], count[s]
        }
    }
    ' "$f"
done > ANALYSIS-per-arm/microsatellites.tsv
