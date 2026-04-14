# Methylation Variation Across 142 Human Y Chromosomes

This repository contains the code and analysis workflow used to investigate DNA methylation patterns across 142 human Y chromosomes, with a particular focus on palindromic regions, structural architecture, and sequence-derived features.

The project aims to characterize variation in methylation levels and methylation heterogeneity across different Y chromosome sequence contexts, and to relate these patterns to intra-palindromic sequence identity, local coverage profiles, and structural organization.

---

## Data used for the analysis

The analysis was done using:
- sample-specific assemblies, available [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/);
- aligned ONT data, available [...] ;
- sample-specific sequence classes annotation, available [...] .

---

## Overview of analyses

The analyses currently included in this repository are:

- **Methylation calling** using `modkit`
- **Methylation heterogeneity quantification** using `modkit`
- **Intra-palindromic arm identity analysis**
- **G-quadruplex (G4) annotation** using nbmst and G4hunter
- exploratory analysis of **methylation across Y chromosome architectural classes**
- **Palindromic structure characterization** using seqrequester

An additional exploratory analysis of **BAM-based coverage** using raaqa was performed but is not part of the core workflow.

---

## Repository structure

```text
.
├── README.md
├── scripts/
│   ├── Python_codes/
|   │   ├── Mariateresa_codes/
|   │   │   ├── Extract_FASTA_sequences.sh
|   │   │   ├── P2_alignment.sh
|   │   │   ├── Analyze_BAMs.sh
|   │   │   ├── filter_BED_file_seqname_column.sh
|   │   │   ├── Generate_cytosine_bed.sh
|   │   │   ├── methHet_modbamtools.sh
|   │   │   ├── modkit_entropy.sh
|   │   │   ├── modkit_extract.sh
|   │   │   ├── modkit_pileup.sh
|   │   │   ├── rename_BED_file_seqname_column.sh
|   │   ├── Oliver_codes/
|   │   │   ├── Meth_on_rearranged_regions/
|   │   │   |   ├── meth_with_arms.py
|   │   │   ├── Seqrequester_task/
|   │   │   |   ├── used_scripts.sh
|   │   │   ├── intra_arm_identity/
|   │   │   |   ├── clean_y_annotation.py
|   │   │   |   ├── extract_both_arms.py
|   │   │   |   ├── sliding_window_identity_batch.py
|   │   ├── README_MArie_Kratka_nonB_DNA.txt
│   ├── R_codes/
|   │   ├── Methylation_along_palindromic_arms.R
|   │   ├── Dinucleotide_composition_of_palindromes.R
|   │   ├── Intra_arm_identity.R
|   │   ├── Meth_AZFc_architectures.R
|   │   ├── Meth_DAZ_genes.R
|   │   ├── Meth_all_linear_circular_plots.R
|   │   ├── Meth_heterogeneity_across_gene_copies.R
|   │   ├── Methylation_across_palindromic_copies.R
|   │   ├── Methylation_across_sequence_classes.R
|   │   ├── Methylation_heterogeneity_across_sequence_classes.R
|   │   ├── meth_ampliconic_genes.R
|   │   ├── nonB_G4_analysis.R

