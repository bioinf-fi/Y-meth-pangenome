### ------- || DAZ GENES || --------
#### ----- || Upload data || -----
filenames_genes = list.files(path = "Work/annotation/GeneAnnotationFiles/CombinedFiles/", full.names = T)
ldf_genes = lapply(filenames_genes,read.csv)
AZFc_inversions = read.csv("Work/annotation/AZFc_inversions.csv", sep = ",", stringsAsFactors = F, header = T)
AZFc_structures = read_tsv("Work/annotation/UniqueAZFc_Structures.tsv")
#### ----- || Organize data || -----
final_daz_df <- ldf_genes %>%
  # 1. Filter for DAZ AND force sample to character
  map(~ .x %>% 
        filter(grepl("DAZ", RM_Amplicon_gene)) %>% 
        mutate(sample = as.character(sample))) %>%
  # 2. Now list_rbind will work smoothly
  list_rbind(names_to = "source_list_name")
# Clean the single combined dataframe
# Clean the dataframes first
final_daz_df <- final_daz_df %>%
  filter(!grepl("random", contig, ignore.case = TRUE)) %>%
  mutate(contig = as.character(contig))
AZFc_inversions <- AZFc_inversions %>%
  filter(!grepl("random", ContigID, ignore.case = TRUE)) %>%
  mutate(ContigID = as.character(ContigID))
final_daz_df$contig <- "chrY"
AZFc_inversions$ContigID <- "chrY"
#### ----- || Methylation on DAZ genes based on their color block location || -----
# Ensure DAZ gene names are explicit (e.g., DAZ1, DAZ2)
final_daz_df = final_daz_df[which(final_daz_df$PangenomeGraph=="CALLED"),]
final_daz_df = final_daz_df[which(final_daz_df$sample!="HG04187"),]
# and clean any "random" contigs if not already done
daz_metadata <- final_daz_df %>%
  filter(!grepl("random", contig, ignore.case = TRUE)) %>%
  mutate(contig = "chrY") 
# Split by sample
genes_list <- split(daz_metadata, daz_metadata$sample)
p_list <- split(plot_df, plot_df$individual)
results_list <- map(common_samples, function(s) {
  # Get data for this sample
  g_df <- genes_list[[s]]
  p_df <- p_list[[s]]
  
  # --- SAFETY CHECK ---
  # Check if g_df is NULL, has 0 rows, or is missing the required columns
  required_cols <- c("start_merged", "end_merged", "contig")
  
  if (is.null(g_df) || nrow(g_df) == 0) {
    message(paste("Skipping sample:", s, "- No DAZ data found."))
    return(NULL)
  }
  
  if (!all(required_cols %in% colnames(g_df))) {
    message(paste("Skipping sample:", s, "- Missing required columns:", 
                  paste(setdiff(required_cols, colnames(g_df)), collapse=", ")))
    return(NULL)
  }
  # --------------------
  
  # 1. Create GRanges for genes
  g_gr <- makeGRangesFromDataFrame(g_df, keep.extra.columns = TRUE,
                                   seqnames.field = "contig", 
                                   start.field = "start_merged", 
                                   end.field = "end_merged")
  
  # 2. Create a MINIMAL p_gr to avoid metadata name conflicts
  p_gr <- GRanges(
    seqnames = p_df$seqnames, 
    ranges = IRanges(start = p_df$end, end = p_df$end)
  )
  
  # 3. Find overlaps
  hits <- findOverlaps(p_gr, g_gr, ignore.strand = TRUE)
  
  # 4. Extract metadata from the gene hits
  overlap_metadata <- data.frame(
    query_idx = queryHits(hits),
    gene_type = mcols(g_gr)$RM_Amplicon_gene[subjectHits(hits)],
    completeness = mcols(g_gr)$RM_Amplicon_geneCompleteness[subjectHits(hits)]
  )
  
  # 5. Join back to the ORIGINAL p_df using row index
  p_df_annotated <- p_df %>%
    mutate(idx = 1:n()) %>%
    inner_join(overlap_metadata, by = c("idx" = "query_idx")) %>%
    select(-idx)
  
  return(p_df_annotated)
})

# Filter out the NULLs (skipped samples) before binding
combined_daz_data <- bind_rows(results_list)

if(exists("blocks_order_v2")) {
  clean_levels <- unique(blocks_order_v2)
  combined_daz_data$blocks_order <- factor(combined_daz_data$blocks_order, levels = clean_levels)
}
ggplot(combined_daz_data, aes(x = factor(blocks_order), y = avg_meth, fill = palindrome)) +
  # Violin plots with median line
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE, draw_quantiles = 0.5) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.5) +
  # Facet by the DAZ gene type (DAZ1, DAZ2, etc.)
  facet_wrap(~ gene_type, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_fill_viridis_d(option = "mako") + # Professional palette for blocks
  labs(
    title = "DAZ Methylation Distribution per Genomic Block",
    subtitle = "Faceted by DAZ gene type; Median indicated by horizontal line",
    x = "Genomic Blocks (Ordered)",
    y = "Average Methylation",
    fill = "Palindrome"
  )
#### ----- || Match DAZ coordinates with inversions || -----
# Create GRanges
gr_daz <- makeGRangesFromDataFrame(final_daz_df, keep.extra.columns = TRUE,
                                   seqnames.field = "contig", 
                                   start.field = "start_merged", 
                                   end.field = "end_merged")
gr_inv <- makeGRangesFromDataFrame(AZFc_inversions, keep.extra.columns = TRUE,
                                   seqnames.field = "ContigID", 
                                   start.field = "Start.coord..proximal.repeat.", 
                                   end.field = "End.coord..distal.repeat.")
# 3. Find all overlaps
hits <- findOverlaps(gr_daz, gr_inv)
# Convert hits to a dataframe to filter by matching Sample IDs
overlap_df <- data.frame(
  daz_idx = queryHits(hits),
  inv_idx = subjectHits(hits),
  daz_sample = gr_daz$sample[queryHits(hits)],
  inv_sample = gr_inv$Sample.with.inversion[subjectHits(hits)]
)
# 4. Keep only hits where the sample IDs match exactly
true_matches <- overlap_df %>% 
  filter(as.character(daz_sample) == as.character(inv_sample))
# Add the inversion status back to your original dataframe
final_daz_df$is_inverted <- FALSE
final_daz_df$is_inverted[true_matches$daz_idx] <- TRUE
# Optional: Add the specific inversion type name
final_daz_df$inversion_type <- NA
final_daz_df$inversion_type[true_matches$daz_idx] <- gr_inv$Inversion[true_matches$inv_idx]
### ----- || Methylation on DAZ genes with no inversions vs. with inversions || -----
# Convert your final DAZ dataframe to GRanges
gr_daz_status <- makeGRangesFromDataFrame(final_daz_df, 
                                          keep.extra.columns = TRUE,
                                          seqnames.field = "contig",
                                          start.field = "start_merged",
                                          end.field = "end_merged")
# Ensure ldf is named by sample to match final_daz_df$sample
# names(ldf) should match the values in final_daz_df$sample
meth_results <- lapply(names(ldf), function(s_name) {
  
  # A. Get methylation data for this sample
  meth_df <- ldf[[s_name]]
  gr_meth <- makeGRangesFromDataFrame(meth_df, 
                                      keep.extra.columns = TRUE, 
                                      seqnames.field = "chr")
  
  # B. Get DAZ regions for this specific individual
  this_sample_daz <- gr_daz_status[gr_daz_status$sample == s_name]
  
  if(length(this_sample_daz) == 0) return(NULL)
  
  # C. Find overlaps
  hits <- findOverlaps(this_sample_daz, gr_meth)
  
  # --- FIX: Check if hits are empty before creating the dataframe ---
  if (length(hits) == 0) {
    return(NULL) 
  }
  
  # D. Calculate average 'meth'
  daz_meth_summary <- data.frame(
    sample = s_name,
    gene_name = this_sample_daz$gene_name[queryHits(hits)],
    is_inverted = this_sample_daz$is_inverted[queryHits(hits)],
    meth_val = gr_meth$meth[subjectHits(hits)]
  ) %>%
    group_by(sample, gene_name, is_inverted) %>%
    summarize(avg_meth = mean(meth_val, na.rm = TRUE), .groups = 'drop')
  
  return(daz_meth_summary)
})
# Combine into one final comparison table
comparison_df <- do.call(rbind, meth_results)
# plotting
# 1. Clean up potential NAs from the previous step
plot_data <- comparison_df[!is.na(comparison_df$avg_meth), ]
# 2. Create the plot using only ggplot2
ggplot(comparison_df, aes(x = is_inverted, y = avg_meth, fill = is_inverted)) +
  geom_violin(trim = FALSE, alpha = 0.6) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#404080"),
                    name = "DAZ Inverted?") +
  labs(title = "DAZ Gene Methylation by Inversion Status",
       subtitle = "Individual-specific overlaps; Mann-Whitney U test",
       x = "Inversion Present",
       y = "Average Methylation Level") +
  theme_minimal()
# statistical test
# T-test to see if there is a significant difference
t.test(avg_meth ~ is_inverted, data = comparison_df)
### ----- || Methylation on DAZ genes which show deletions vs. full genes || -----
samples_w_gr_del = c("HG005", "NA18952", "NA18983", "HG01099", "HG02135", "HG02129", "HG04228")
# 1. Convert both to GRanges
gr_meth <- makeGRangesFromDataFrame(combined_data, keep.extra.columns = TRUE, 
                                    seqnames.field = "chr")
gr_metadata <- makeGRangesFromDataFrame(final_daz_df, keep.extra.columns = TRUE,
                                        seqnames.field = "contig", 
                                        start.field = "start_merged", 
                                        end.field = "end_merged")
# 2. Find overlaps (handling potential coordinate differences)
hits <- findOverlaps(gr_meth, gr_metadata)
# 3. Create the merged dataframe using the hits
# We match by overlap AND by sample name
overlap_df <- data.frame(
  meth = mcols(gr_meth)$meth[queryHits(hits)],
  sample_meth = mcols(gr_meth)$sampleName[queryHits(hits)],
  sample_meta = mcols(gr_metadata)$sample[subjectHits(hits)],
  gene_type = mcols(gr_metadata)$RM_Amplicon_gene[subjectHits(hits)],
  completeness = mcols(gr_metadata)$RM_Amplicon_geneCompleteness[subjectHits(hits)]
) %>%
  filter(as.character(sample_meth) == as.character(sample_meta)) # Ensure same individual
# 1. Filter for DAZ1 specifically
daz1_data <- overlap_df %>%
  filter(gene_type == "DAZ1")
# 2. Create the Violin Plot
ggplot(daz1_data, aes(x = completeness, y = meth, fill = completeness)) +
  # Adds the violin shape
  geom_violin(alpha = 0.6, trim = FALSE) +
  # Adds a narrow boxplot inside to show medians/quartiles
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Optional: adds jittered points to show actual data distribution
  geom_jitter(width = 0.05, alpha = 0.2, size = 1) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Methylation Distribution: DAZ1 Truncated vs. Full",
    subtitle = "Based on Individual-specific overlaps",
    x = "Gene Completeness Status",
    y = "Methylation Level",
    fill = "Status"
  ) +
  theme(legend.position = "none") # Legend is redundant if X-axis is labeled
# 3. Optional: Quick Statistical Check for DAZ1
wilcox.test(meth ~ completeness, data = daz1_data)
### ----- || Methylation of DAZ genes on individuals with different numbers of copies || -----
# Standardize IDs to ensure a clean match
gene_copies$Sample.ID <- as.character(gene_copies$Sample.ID)
comparison_df$sample <- as.character(comparison_df$sample)
# Join the copy number information
# Note: We specifically select Sample.ID and DAZ (total copies)
comparison_with_copies <- comparison_df %>%
  left_join(gene_copies %>% select(Sample.ID, DAZ_count = DAZ), 
            by = c("sample" = "Sample.ID"))
# Convert DAZ_count to a factor for plotting
comparison_with_copies$DAZ_count <- as.factor(comparison_with_copies$DAZ_count)
ggplot(comparison_with_copies, aes(x = avg_meth, fill = DAZ_count)) +
  geom_density(alpha = 0.4) +
  # Using a scale that handles multiple categories well
  scale_fill_viridis_d(name = "DAZ Copies") +
  labs(title = "DAZ Methylation Distribution by Copy Number",
       x = "Average Methylation Level",
       y = "Density") +
  theme_minimal()
### ----- || Methylation on distinct copies of DAZ genes || -----
# 1. Rank the DAZ copies within each individual
daz_ranked <- final_daz_df %>%
  group_by(sample) %>%
  # Sort by start coordinate to rank them upstream -> downstream
  arrange(start_merged, .by_group = TRUE) %>%
  mutate(
    copy_rank = row_number(),
    copy_label = paste0("Copy ", copy_rank)
  ) %>%
  ungroup()
# Convert this ranked dataframe to GRanges
gr_daz_ranked <- makeGRangesFromDataFrame(daz_ranked, 
                                          keep.extra.columns = TRUE,
                                          seqnames.field = "contig",
                                          start.field = "start_merged",
                                          end.field = "end_merged")
# 2. Rank methylation by copy
ranked_meth_results <- lapply(names(ldf), function(s_name) {
  # A. Get individual methylation and DAZ data
  meth_df <- ldf[[s_name]]
  gr_meth <- makeGRangesFromDataFrame(meth_df, keep.extra.columns = TRUE, seqnames.field = "chr")
  
  this_sample_daz <- gr_daz_ranked[gr_daz_ranked$sample == s_name]
  if(length(this_sample_daz) == 0) return(NULL)
  
  # B. Find overlaps
  hits <- findOverlaps(this_sample_daz, gr_meth)
  if (length(hits) == 0) return(NULL)
  
  # C. Create summary including the rank
  data.frame(
    sample = s_name,
    copy_rank = this_sample_daz$copy_rank[queryHits(hits)],
    copy_label = this_sample_daz$copy_label[queryHits(hits)],
    is_inverted = this_sample_daz$is_inverted[queryHits(hits)],
    meth_val = gr_meth$meth[subjectHits(hits)]
  ) %>%
    group_by(sample, copy_rank, copy_label, is_inverted) %>%
    summarize(avg_meth = mean(meth_val, na.rm = TRUE), .groups = 'drop')
})
final_ranked_df <- do.call(rbind, ranked_meth_results)
# 3. plotting
ggplot(final_ranked_df, aes(x = factor(copy_rank), y = avg_meth)) +
  # 1. Boxplot to show distribution per rank
  geom_boxplot(fill = "#69b3a2", alpha = 0.5, outlier.shape = NA) +
  
  # 2. Simple jittered points (since we aren't dodging by inversion anymore)
  geom_jitter(width = 0.2, alpha = 0.4, color = "#4d8c7d") +
  
  # 3. Aesthetics
  labs(
    title = "DAZ Methylation by Genomic Position",
    subtitle = "Ranked 1 (most upstream) to N (most downstream)",
    x = "Copy Rank (Position along Chromosome)",
    y = "Average Methylation"
  ) +
  theme_minimal()
