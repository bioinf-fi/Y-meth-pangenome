# Methylation Variation Across 143 Human Y Chromosomes

This repository contains the code and analysis workflow used to investigate DNA methylation patterns across 143 human Y chromosomes, with a particular focus on palindromic regions, structural architecture, and sequence-derived features.

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
- **G-quadruplex (G4) annotation**
- **BAM-based coverage visualization** using external alignment inspection tools (e.g. IGV, samtools)
- **Palindromic structure characterization**

An additional exploratory analysis of methylation across Y chromosome architectural classes was performed but is not part of the core workflow.

---

## Repository structure

```text
.
├── README.md
├── scripts/
│   ├── methylation_calling/
│   ├── methylation_heterogeneity/
│   ├── palindrome_identity/
│   ├── g4_annotation/
│   └── palindrome_structure/
├── metadata/
├── results/
│   ├── methylation/
│   ├── heterogeneity/
│   ├── identity/
│   ├── g4/
│   ├── coverage/
│   └── structure/
├── figures/
└── docs/
