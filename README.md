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

## R session and packages versions

> sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Biostrings_2.66.0     XVector_0.38.0        magick_2.9.1          ComplexHeatmap_2.14.0
 [5] magrittr_2.0.4        circlize_0.4.17       valr_0.9.1            patchwork_1.3.2      
 [9] ggridges_0.5.7        plyranges_1.18.0      rtracklayer_1.58.0    Gviz_1.42.1          
[13] data.table_1.18.2.1   corrplot_0.95         lubridate_1.9.5       forcats_1.0.0        
[17] readr_2.2.0           tidyverse_2.0.0       tidyr_1.3.2           stringr_1.6.0        
[21] dplyr_1.2.0           GenomicRanges_1.50.2  GenomeInfoDb_1.34.9   IRanges_2.32.0       
[25] S4Vectors_0.36.2      BiocGenerics_0.44.0   ggplot2_4.0.2         purrr_1.2.1          
[29] scales_1.4.0          RColorBrewer_1.1-3    zoo_1.8-11            ape_5.8-1            
[33] phylogram_2.1.0       ggnewscale_0.5.2      ggtreeExtra_1.21.0    ggtree_4.1.1.005     
[37] tibble_3.3.1         

loaded via a namespace (and not attached):
  [1] backports_1.5.0             Hmisc_5.2-5                 BiocFileCache_2.6.1        
  [4] systemfonts_1.3.2           lazyeval_0.2.2              BiocParallel_1.32.6        
  [7] digest_0.6.39               foreach_1.5.2               ensembldb_2.22.0           
 [10] yulab.utils_0.2.4           htmltools_0.5.9             checkmate_2.3.4            
 [13] memoise_2.0.1               BSgenome_1.66.3             doParallel_1.0.17          
 [16] cluster_2.1.4               tzdb_0.3.0                  matrixStats_1.5.0          
 [19] vroom_1.7.0                 timechange_0.4.0            prettyunits_1.2.0          
 [22] jpeg_0.1-11                 colorspace_2.1-2            blob_1.3.0                 
 [25] rappdirs_0.3.4              xfun_0.57                   crayon_1.5.3               
 [28] RCurl_1.98-1.18             jsonlite_2.0.0              iterators_1.0.14           
 [31] VariantAnnotation_1.44.1    glue_1.8.0                  gtable_0.3.6               
 [34] zlibbioc_1.44.0             GetoptLong_1.1.0            DelayedArray_0.24.0        
 [37] shape_1.4.6.1               fontquiver_0.2.1            DBI_1.3.0                  
 [40] Rcpp_1.1.1                  progress_1.2.3              htmlTable_2.4.3            
 [43] clue_0.3-67                 gridGraphics_0.5-1          tidytree_0.4.7             
 [46] foreign_0.8-83              bit_4.6.0                   Formula_1.2-5              
 [49] fontLiberation_0.1.0        htmlwidgets_1.6.4           httr_1.4.8                 
 [52] cpp11bigwig_0.1.3           pkgconfig_2.0.3             XML_3.99-0.23              
 [55] farver_2.1.2                deldir_2.0-4                nnet_7.3-18                
 [58] dbplyr_2.5.2                labeling_0.4.3              ggplotify_0.1.3            
 [61] tidyselect_1.2.1            rlang_1.1.7                 AnnotationDbi_1.60.2       
 [64] tools_4.2.2                 cachem_1.1.0                cli_3.6.5                  
 [67] generics_0.1.4              RSQLite_2.4.6               broom_1.0.3                
 [70] evaluate_1.0.5              fastmap_1.2.0               yaml_2.3.12                
 [73] knitr_1.51                  bit64_4.6.0-1               fs_1.6.7                   
 [76] KEGGREST_1.38.0             AnnotationFilter_1.22.0     nlme_3.1-160               
 [79] ggiraph_0.9.6               aplot_0.2.9                 xml2_1.5.2                 
 [82] biomaRt_2.54.1              compiler_4.2.2              rstudioapi_0.18.0          
 [85] filelock_1.0.3              curl_7.0.0                  png_0.1-9                  
 [88] treeio_1.35.0               stringi_1.8.7               GenomicFeatures_1.50.4     
 [91] gdtools_0.5.0               lattice_0.20-45             ProtGenerics_1.30.0        
 [94] Matrix_1.5-1                fontBitstreamVera_0.1.1     vctrs_0.7.1                
 [97] pillar_1.11.1               lifecycle_1.0.5             GlobalOptions_0.1.3        
[100] bitops_1.0-9                latticeExtra_0.6-31         R6_2.6.1                   
[103] BiocIO_1.8.0                gridExtra_2.3               codetools_0.2-18           
[106] dichromat_2.0-0.1           MASS_7.3-58.1               SummarizedExperiment_1.28.0
[109] rjson_0.2.23                withr_3.0.2                 GenomicAlignments_1.34.1   
[112] Rsamtools_2.14.0            GenomeInfoDbData_1.2.9      parallel_4.2.2             
[115] hms_1.1.4                   rpart_4.1.19                ggfun_0.2.0                
[118] rmarkdown_2.30              S7_0.2.1                    MatrixGenerics_1.10.0      
[121] biovizBase_1.46.0           Biobase_2.58.0              base64enc_0.1-6            
[124] interp_1.1-6                restfulr_0.0.16  

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

