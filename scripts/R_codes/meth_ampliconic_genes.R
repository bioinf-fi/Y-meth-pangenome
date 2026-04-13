## --------- || CORRELATION AMONG AMPLICONIC GENE COPIES BASED ON LOCATION || -----------
###  ----- || Create a table combining palindromic information and gene information || ---------
genes_df = read.csv("Work/annotation/AmplGene_coords_Pille_Hallast/HGSVC3_HPRC_CEPH-AmpliconicGenes_combined_wLiftOff_plusManualAnnotation.09302025_AllExonGenes_noDbGaP.csv", sep = ",", stringsAsFactors = F, header = T)
genes_df <- genes_df[!grepl("random", genes_df$contig), ]
genes_df$contig = rep("chrY", nrow(genes_df))
genes_df = genes_df[,-1]
genes_df = data.frame("seqnames" = genes_df$contig, "start" = genes_df$contigStart, "end" = genes_df$contigEnd, "strand" = genes_df$orientation, genes_df[,c(1:2,4:7,11:17)])
genes_df$strand[which(genes_df$strand=="C")] = "-"
# 1. Split into lists by sample
genes_list <- split(genes_df, genes_df$sampleName)
common_samples <- intersect(names(genes_list), names(ldf))
results_list <- map(common_samples, function(s) {
  
  g_df <- genes_list[[s]]
  p_df <- ldf[[s]]
  
  g_gr <- makeGRangesFromDataFrame(g_df, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
  p_gr <- makeGRangesFromDataFrame(p_df, keep.extra.columns = TRUE) 
  
  hits <- findOverlaps(p_gr, g_gr, ignore.strand = TRUE)
  
  overlap_metadata <- data.frame(
    query_idx = queryHits(hits),
    geneName = mcols(g_gr)$geneName[subjectHits(hits)],
    geneCompletenessInfo = mcols(g_gr)$geneCompletenessInfo[subjectHits(hits)],
    totalExons = mcols(g_gr)$totalExons[subjectHits(hits)]
  )
  
  overlap_summary <- overlap_metadata %>%
    group_by(query_idx) %>%
    summarise(
      gene_names = paste(unique(geneName), collapse = "; "),
      exon_info = paste(unique(geneCompletenessInfo), collapse = "; "),
      exon_count = paste(unique(totalExons), collapse = "; "),
      has_overlap = TRUE, # Renamed to avoid confusion with pre-existing column name
      .groups = "drop"
    )
  
  # FIX: Ensure the base column for overlap status is initialized as FALSE 
  # before merging, and merge using the new column name 'has_overlap'.
  p_df <- p_df %>%
    mutate(
      idx = 1:n(),
      sampleName = s,  
      # If the original raw_plot_df already has an 'overlaps_gene' column, 
      # we make sure we don't confuse the script by initializing a new temporary one.
      overlaps_gene_temp = FALSE 
    ) %>%
    left_join(overlap_summary, by = c("idx" = "query_idx")) %>%
    # Use coalesce to pick TRUE from has_overlap if present, otherwise keep the default FALSE
    mutate(
      overlaps_gene = coalesce(has_overlap, overlaps_gene_temp)
    ) %>%
    select(-idx, -has_overlap, -overlaps_gene_temp) # Clean up temporary columns
  
  return(p_df)
})
combined_data <- bind_rows(results_list)
combined_data = combined_data[!is.na(combined_data$gene_names),]
### ------- || Create plot to show methylation changes on ampliconic genes based on their color block location || --------
# 1. Prepare the lists
# Using final_daz_df which we cleaned earlier
genes_list <- split(genes_df, genes_df$sampleName)
# Using your existing plot_df (binned methylation data)
colnames(plot_df)[2] = "seqnames"
pal_list <- split(plot_df, plot_df$individual) 
common_samples <- intersect(names(genes_list), names(pal_list))
# 2. Map and Annotate
results_list <- map(common_samples, function(s) {
  g_df <- genes_list[[s]]
  p_df <- pal_list[[s]]
  
  # Ensure DAZ coordinates are GRanges
  g_gr <- makeGRangesFromDataFrame(g_df, keep.extra.columns = TRUE,
                                   seqnames.field = "seqnames",
                                   start.field = "start",
                                   end.field = "end")
  
  # Ensure binned data has chr, start, end for GRanges
  p_gr <- makeGRangesFromDataFrame(p_df, keep.extra.columns = TRUE) 
  
  # Find overlaps between bins and DAZ genes
  hits <- findOverlaps(p_gr, g_gr, ignore.strand = TRUE)
  
  overlap_metadata <- data.frame(
    query_idx = queryHits(hits),
    geneName = mcols(g_gr)$geneName[subjectHits(hits)],
    geneCompleteness = mcols(g_gr)$geneCompletenessInfo[subjectHits(hits)]
  )
  
  # Join gene names to the binned data
  p_df_annotated <- p_df %>%
    mutate(idx = 1:n()) %>%
    left_join(overlap_metadata, by = c("idx" = "query_idx")) %>%
    filter(!is.na(gene)) %>% 
    select(-idx)
  
  return(p_df_annotated)
})

combined_gene_data <- bind_rows(results_list)

# 3. Factorize Blocks for Ordering
if(exists("blocks_order_v2")) {
  clean_levels <- unique(blocks_order_v2)
  combined_gene_data <- combined_gene_data %>%
    mutate(blocks_order = factor(blocks_order, levels = clean_levels))
}

# 4. Final Plotting
ggplot(combined_gene_data, aes(x = blocks_order, y = avg_meth, fill = palindrome)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE, draw_quantiles = 0.5) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.5) +
  # Facet by the specific gene type
  facet_wrap(~ geneName, ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Gene Methylation Distribution per Genomic Block",
    subtitle = "Faceted by Gene Name; Medians marked by horizontal lines",
    x = "Genomic Blocks (Ordered)",
    y = "Average Methylation",
    fill = "Palindrome"
  )
###  ----- || View ampliconic genes location || ---------
# Ensure names are consistent within the annotation list
ldf_ann_flat <- ldf_ann %>%
  imap_dfr(~ {
    .x %>%
      # Ensure seqnames is used instead of chr
      rename(seqnames = any_of("chr")) %>%
      # Add the individual ID from the list name (.y)
      mutate(individual_id = .y) %>%
      # Keep relevant columns for the map
      select(individual_id, seqnames, start, end, sequence_class, any_of("class"))
  })
# Create a clean mapping from (individual, class) -> start coordinate
palindrome_start_map <- ldf_ann_flat %>%
  group_by(individual_id, sequence_class) %>%
  summarise(
    palindrome_start = min(start, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(sampleName = individual_id) # Match the plot code names
# Join genes_df with the map to get the start of the parent palindrome
colnames(combined_gene_data)[c(6,9)] = c("sequence_class", "sampleName")
genes_normalized_df <- combined_gene_data %>%
  rename(sampleName = any_of(c("Sample.ID", "sampleName"))) %>% # Standardize ID name
  left_join(palindrome_start_map, by = c("sampleName", "sequence_class")) %>%
  mutate(
    # Normalize by subtracting the palindrome start
    normalized_start = start - palindrome_start,
    normalized_end = end - palindrome_start
    # normalized length is simply: normalized_end - normalized_start
  )
ggplot(genes_normalized_df, aes(color = sampleName)) +
  # Map normalized coordinates to the X axis
  geom_linerange(aes(xmin = normalized_start, xmax = normalized_end, y = geneName),
                 position = position_dodge(width = 0.7), 
                 alpha = 0.8,
                 size = 2) + 
  theme_minimal() +
  labs(
    title = "Normalized Gene Locations within Palindromes",
    x = "Relative Position (bp from Palindrome Start)",
    y = "Gene Name",
    color = "Individual"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    # Rotate Y axis labels to match your request for "lateral view"
    axis.text.y = element_text(angle = 0, hjust = 1) 
  ) +
  # Arrange into a single column of plots (one plot per gene)
  facet_wrap(~ geneName, ncol = 1, scales = "free_y")



ggplot(genes_df, aes(x = geneName, color = sampleName)) +
  # Use ymin and ymax to define the start and end of the range
  geom_linerange(aes(ymin = start, ymax = end), position = position_dodge(width = 0.5), alpha = 1) + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    title = "Palindrome Location Conservation and Length (2026)",
    x = "Gene Name(s)",
    y = "Genomic Position" # Y-axis now represents the full range
  ) 
# Removed: scale_x_discrete(limits = paste0("P", 1:8))
###  ----- || View ampliconic genes location in individuals with different copies - VCY e BPY2 || ---------
gene_copies = read.table("Work/annotation/AmplGene_coords_Pille_Hallast/GeneCopies_HGSVC3_assemblies.csv", sep = ",", stringsAsFactors = F, header = T)
rownames(gene_copies) = gene_copies$Sample.ID
gene_copies = gene_copies[,c(1,4:17)]
# We need to pivot this copy number data so each row is a gene-sample pair
gene_copies_long <- gene_copies %>%
  pivot_longer(cols = c(BPY2, VCY), names_to = "geneName", values_to = "copy_number")
gene_copies_long = gene_copies_long[,c(1,14:15)]
colnames(gene_copies_long)[1] = "sampleName"
# === Merge Location Data (genes_df) with Copy Number Data (gene_copies_long) ===
# The 'genes_df' must have a column named 'Sample.ID' (or 'sampleName') and 'geneName'
merged_df <- left_join(genes_df, gene_copies_long, by = c("sampleName", "geneName"))
# === Filter and Plotting ===
# Filter for relevant genes and copy number range [1, max_copies]
max_copies <- max(merged_df$copy_number, na.rm = TRUE)
filtered_data <- merged_df %>%
  filter(geneName %in% c("BPY2", "VCY"), 
         copy_number >= 1, 
         copy_number <= max_copies)
# Assume 'merged_df' contains columns: Sample.ID, geneName, start, end, copy_number
# reate the summarized gene range table
gene_ranges_summary <- filtered_data %>%
  group_by(sampleName, geneName, copy_number) %>%
  summarise(
    gene_start = min(start, na.rm = TRUE),
    gene_end = max(end, na.rm = TRUE),
    .groups = "drop" # Ensures the output is a standard dataframe
  )
print(head(gene_ranges_summary))
ggplot(gene_ranges_summary, aes(x = factor(copy_number), color = sampleName)) +
  geom_linerange(aes(ymin = gene_start, ymax = gene_end), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.7,
                 linewidth = 1.5) + 
  theme(legend.position = "none") +
  labs(
    title = "Full Gene Ranges by Copy Number for BPY2 and VCY",
    x = "Gene Copy Number (Count)",
    y = "Genomic Position (bp)",
    color = "Individual Sample"
  ) +
  # Separate plots for each gene name (BPY2/VCY)
  facet_grid(. ~ geneName) 
###  ----- || View methylation on distinct copies - VCY e BPY2 || ---------
# Instead of your summarise() command, use this:
gene_ranges_individual <- filtered_data %>%
  # Group by sample and gene to identify copies within each person
  group_by(sampleName, geneName) %>%
  # Sort by starting position so 1st is the most upstream
  arrange(start, .by_group = TRUE) %>%
  # Assign the 'copy' rank (1, 2, 3...)
  mutate(
    copy_rank = row_number(),
    copy_label = paste0("Copy ", copy_rank)
  ) %>%
  # Rename columns to match your anchor logic if necessary
  rename(gene_start = start, gene_end = end) %>%
  ungroup()
colnames(combined_data)[6] = "geneName"
# Merge ranked anchor with CpG data
ranked_meth_summary <- gene_ranges_individual %>%
  inner_join(combined_data, 
             by = join_by(sampleName, geneName, 
                          gene_start <= start, gene_end >= end)) %>%
  group_by(sampleName, geneName, copy_rank, copy_label)
ggplot(ranked_meth_summary, aes(x = factor(copy_rank), y = meth, fill = geneName)) +
  # Use boxplots to show distribution across individuals for each copy rank
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  # Overlay jittered points to see individual data
  facet_wrap(~geneName, scales = "free_x") +
  labs(
    title = "Methylation by Gene Copy Occurrence",
    subtitle = "Ranked from 1st (most upstream) to last (most downstream)",
    x = "Copy Occurrence (Position Rank)",
    y = "Average Methylation",
    fill = "Gene"
  ) +
  theme_minimal()
###  ----- || Create a table to see areas covered vs. not covered by genes || ---------
colnames(all_individuals_data)[1] = "seqnames"
pal_list <- split(all_individuals_data, all_individuals_data$individual_id)
common_samples <- intersect(names(genes_list), names(pal_list))
# 2. Iterate through each sample and compute overlaps
#### Split into lists by sample
results_list <- map(common_samples, function(s) {
  
  g_df <- genes_list[[s]]
  p_df <- pal_list[[s]]
  
  # SAFETY CHECK: Skip if the palindrome data is empty for this sample
  if (is.null(p_df) || nrow(p_df) == 0) return(NULL)
  
  # Convert to GRanges - EXPLICITLY define columns to prevent ".find_start_end_cols" error
  p_gr <- makeGRangesFromDataFrame(p_df, 
                                   seqnames.field = "seqnames", 
                                   start.field = "start", 
                                   end.field = "end", 
                                   keep.extra.columns = TRUE)
  
  # Check if genes exist for this sample; if not, return p_df with defaults
  if (is.null(g_df) || nrow(g_df) == 0) {
    p_df$overlaps_gene <- FALSE
    p_df$geneName <- NA_character_
    return(p_df)
  }
  
  g_gr <- makeGRangesFromDataFrame(g_df, 
                                   seqnames.field = "seqnames", 
                                   start.field = "start", 
                                   end.field = "end", 
                                   keep.extra.columns = TRUE)
  
  # Find Overlaps
  hits <- findOverlaps(p_gr, g_gr, ignore.strand = TRUE)
  
  # Initialize columns
  p_df$overlaps_gene <- FALSE
  p_df$geneName <- NA_character_
  
  if (length(hits) > 0) {
    p_df$overlaps_gene[unique(queryHits(hits))] <- TRUE
    
    # Map and collapse gene names
    overlap_map <- data.frame(
      p_idx = queryHits(hits),
      g_name = as.character(g_df$geneName[subjectHits(hits)])
    ) %>%
      group_by(p_idx) %>%
      summarise(all_genes = paste(unique(g_name), collapse = ", "), .groups = "drop")
    
    p_df$geneName[overlap_map$p_idx] <- overlap_map$all_genes
  }
  
  return(p_df)
})

# Use compact() to remove any NULL entries (samples we skipped)
raw_plot_df <- bind_rows(purrr::compact(results_list))
# 4. Final Plotting Preparation
raw_plot_df$overlap_group <- factor(
  raw_plot_df$overlaps_gene,
  levels = c(FALSE, TRUE),
  labels = c("No Gene Overlap", "Covered by Gene")
)
# 3. Create a violin plot using ggplot2
ggplot(plot_df, aes(x = overlap_group, y = avg_meth, fill = overlap_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) + # trim=FALSE extends violins to extreme data points
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) + # Optional: add a mini boxplot inside
  labs(
    title = "Methylation Levels: Gene-Covered vs. Uncovered Palindromic Regions (Violin Plot)",
    x = "Genomic Region Status",
    y = "Average Methylation (%)",
    fill = "Genomic Region Status"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme(legend.position = "none") # Hide legend if x-axis labels are clear enough
wilcox.test(plot_df$avg_meth[which(plot_df$overlaps_gene==TRUE)],plot_df$avg_meth[which(plot_df$overlaps_gene==FALSE)])
### ----------- || Barplot containing number of genes on each palindrome (Works as is) || ----------
# Create a data frame
# Counts how many times each 'palindrome' appears for each 'individual'
occurrence_table <- combined_data %>%
  count(palindrome, geneName)
occurrence_table
palindrome_data <- data.frame(
  Palindrome = c("P7", "P6", "P8", "P5", "P4", "P3", "P2", "P1"),
  Gene_Count = c(0, 0, 1, 3, 1, 2, 2, 5)
)
palindrome_data$Palindrome <- factor(palindrome_data$Palindrome, levels = unique(palindrome_data$Palindrome)) 
ggplot(palindrome_data, aes(x = Palindrome, y = Gene_Count)) +
  geom_col(fill = "steelblue") + 
  labs(title = "Number of Genes per Palindrome",
       x = "Palindrome ID",
       y = "Number of Genes") +
  theme_minimal()
### ----------- || Density plot to see how methylation differs with gene density || ----------
# This table links palindrome IDs (P1-P8) to the known number of genes
genes_P <- data.frame(
  palindrome = c("P8", "P5", "P4", "P3", "P2", "P1"), 
  defined_gene_count = c(1, 3, 1, 2, 2, 5)
)
# --- 2. Join the Mapping to plot_df and Classify ---
plot_df_classified <- plot_df %>%
  dplyr::left_join(genes_P, by = "palindrome") %>%
  na.omit() # Removes spacers or NAs that don't match
# Create the categories for plotting (1 gene, 3 genes, >5 genes)
plot_df_classified$PalindromeClass <- cut(
  plot_df_classified$defined_gene_count,
  # Breaks define the edges. To get exactly 1, 3-4, and 5+:
  breaks = c(0, 1, 4, Inf), 
  labels = c("1 gene", "2-3 genes", "5 genes"),
  right = TRUE
)
# --- 3. Generate the density plot using the classified plot_df ---
ggplot(plot_df_classified, aes(x = avg_meth, color = PalindromeClass, fill = PalindromeClass)) +
  geom_density(alpha = 0.5) + # Smooths the distribution of methylation values
  labs(
    title = "Smoothed Density Plot of Methylation by Palindrome Gene Count",
    x = "Average Methylation per Bin (%)",
    y = "Density"
  ) +
  theme_minimal() +
  # Manually order the legend/colors if needed
  scale_color_discrete(limits = c("1 gene", "2-3 genes", "5 genes")) +
  scale_fill_discrete(limits = c("1 gene", "2-3 genes", "5 genes"))
### ------- || Create boxplot to show methylation changes between complete and truncated exons || --------
# --- ASSUMPTION: 'plot_df' is already generated in your environment ---
# plot_df columns: bin_id, chr, start, end, avg_meth, individual, palindrome, arm
# 1. Filter plot_df to isolate relevant annotations (e.g., specific palindromes/arms)
# We assume 'full exons' and 'truncated exons' relate to some specific 'palindrome' or 'arm' type.
# Based on the image and previous context, let's compare 'left' vs 'right' arms as an example:
# Filter for the arms you want to compare
filtered_plot_df <- plot_df %>%
  filter(arm %in% c("left", "right")) # Adjust this filter based on what defines "full" vs "truncated"
# 2. Reshape the data from wide-by-individual to long-by-arm for plotting
methylation_long <- filtered_plot_df %>%
  rename(Exon_Type = arm) %>% # Rename 'arm' to 'Exon_Type' for consistency with the previous plot labels
  mutate(Exon_Type = recode(Exon_Type,
                            "left" = "Full Exons", # Assuming 'left' maps to 'full' for your analysis
                            "right" = "Truncated Exons")) %>% # Assuming 'right' maps to 'truncated'
  # Aggregate the binned data back to a per-sample, per-type average
  group_by(individual, Exon_Type) %>%
  summarise(Mean_Methylation = mean(avg_meth, na.rm = TRUE), .groups = "drop")
# 3. Statistical Analysis
# Perform a paired t-test to compare the means across all samples
t_test_result <- t.test(Mean_Methylation ~ Exon_Type, 
                        data = methylation_long, 
                        paired = TRUE)
print(t_test_result)
# 4. Generate the Boxplot Visualization
ggplot(methylation_long, aes(x = Exon_Type, y = Mean_Methylation, fill = Exon_Type)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.6, width = 0.1) +
  labs(
    title = "Methylation Levels: Full vs. Truncated Exons (from plot_df)",
    x = "Exon Type",
    y = "Mean Methylation Value",
    fill = "Exon Type"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired")
### ------- || Create boxplot to show methylation changes across gene copies || --------
# --- STEP 1: Spatial Overlap (remains similar but converted to DT) ---
methylation_gr <- GRanges(
  seqnames = plot_df$chr, ranges = IRanges(start = plot_df$start, end = plot_df$end),
  individual = plot_df$individual, avg_meth = plot_df$avg_meth
)
genes_gr <- GRanges(
  seqnames = combined_data$chr, ranges = IRanges(start = combined_data$start, end = combined_data$end),
  geneName = combined_data$geneName, sampleName = combined_data$sampleName
)
overlaps_spatial <- findOverlaps(methylation_gr, genes_gr)
# Create the linked dataframe and immediately convert to data.table
merged_dt <- data.frame(
  individual = genes_gr[subjectHits(overlaps_spatial)]$sampleName,
  avg_meth = methylation_gr[queryHits(overlaps_spatial)]$avg_meth,
  geneName = genes_gr[subjectHits(overlaps_spatial)]$geneName
) %>% distinct() %>% setDT()
# --- STEP 2: The Crash-Fix (Pre-Collapse combined_data) ---
# Create a unique reference table with only essential columns.
# This prevents the "Many-to-Many" row explosion.
combined_ref_dt <- as.data.table(combined_data)[
  !is.na(totalExons) & !is.na(geneCompleteness), 
  .(totalExons = first(totalExons), geneCompleteness = first(geneCompleteness)), 
  by = .(sampleName, geneName)
]
# --- STEP 3: Memory-Efficient Join & Calculation ---
# Join the overlapping hits with the unique reference table
final_analysis_dt <- combined_ref_dt[
  merged_dt, 
  on = .(sampleName = individual, geneName), 
  nomatch = 0
]
# Calculate copy number
gene_copies_per_sample <- final_analysis_dt[
  geneCompleteness == "Full_Exon_Copy", 
  .(n_occurrences = .N), 
  by = .(sampleName, geneName, totalExons)
][, gene_copy_number := n_occurrences / totalExons]
# --- STEP 4: Cleanup ---
rm(merged_dt, combined_ref_dt)
gc() # Manually trigger garbage collection to free RAM
# 4. Use merge() for data.table - more robust against scoping errors
# This performs an inner join (nomatch=0 equivalent)
plot_ready_dt <- merge(
  gene_copies_per_sample, 
  final_analysis_dt, 
  by = c("sampleName", "geneName")
)
# Optional: Ensure column names match for plotting if you used 'individual' earlier
# setnames(plot_ready_dt, "sampleName", "individual") 
# 5. Visualize with 2026 settings for large genomic data
ggplot(plot_ready_dt, aes(x = gene_copy_number, y = avg_meth)) +
  geom_point(alpha = 0.1, size = 0.5) + # Lower alpha and smaller size for density
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE) + 
  labs(
    title = "Methylation Across Gene Copies (Integrated Data)",
    x = "Gene Copy Number (n_occurrences / totalExons)",
    y = "Average Methylation (avg_meth)"
  ) +
  theme_minimal()
### ------- || Create plot to show methylation changes on ampliconic genes in relation to their location on the palindrome || --------
# Assumes 'plot_df' is available from previous steps and has 'arm', 'palindrome', 'start', 'end', 'avg_meth'
# step 1. Define the function to assign the "position" to each gene
assign_position <- function(df) {
  dt <- as.data.table(df)
  
  dt[, `:=`(
    min_c = min(start, na.rm = TRUE),
    max_c = max(end, na.rm = TRUE)
  ), by = .(sampleName, palindrome)]
  
  dt[, range_l := max_c - min_c]
  dt[, `:=`(d_thr = min_c + 0.25 * range_l, c_thr = max_c - 0.25 * range_l)]
  
  dt[, position := "middle"]
  dt[(arm == "left" & start <= d_thr) | (arm == "right" & end >= c_thr), position := "distal"]
  dt[(arm == "left" & end >= c_thr) | (arm == "right" & start <= d_thr), position := "center"]
  
  # Remove temp columns
  dt[, c("min_c", "max_c", "range_l", "d_thr", "c_thr") := NULL]
  
  return(dt)
}
# Apply the function (Assuming your combined_data is loaded)
combined_data_with_position <- assign_position(combined_data)
# step 2. merge it back to plot_df
combine_methylation_and_genes <- function(plot_df, combined_data_with_position) {
  # 1. Ensure both are data.tables
  dt_meth <- as.data.table(plot_df)
  dt_gene <- as.data.table(combined_data_with_position)
  # 2. Handle the "individual" vs "sampleName" naming early
  if ("individual" %in% names(dt_meth)) {
    setnames(dt_meth, "individual", "sampleName")
  }
  # 3. Create GRanges (standard for genomic overlaps)
  gr_meth <- makeGRangesFromDataFrame(dt_meth, keep.extra.columns = TRUE)
  gr_gene <- makeGRangesFromDataFrame(dt_gene, keep.extra.columns = TRUE)
  # 4. Find the intersections
  hits <- findOverlaps(gr_meth, gr_gene)
  # 5. Build the result using explicit indexing
  res <- data.table(
    sample_meth = dt_meth$sampleName[queryHits(hits)],
    sample_gene = dt_gene$sampleName[subjectHits(hits)],
    avg_meth    = dt_meth$avg_meth[queryHits(hits)],
    geneName    = dt_gene$geneName[subjectHits(hits)],
    position    = dt_gene$position[subjectHits(hits)],
    palindrome  = dt_gene$palindrome[subjectHits(hits)] # Added this column
  )
  # 6. Filter for matching samples and cleanup
  final_result <- res[sample_meth == sample_gene, .(
    sampleName = sample_meth, 
    avg_meth, 
    geneName, 
    position,
    palindrome # Added this column
  )]
  # 7. Fast deduplication
  return(unique(final_result))
}
# Apply the combining function (Assuming combined_data_with_position is ready)
combined_df_for_plot <- combine_methylation_and_genes(plot_df, combined_data_with_position)
# step 3. plotting the density plot
plot_methylation_density <- function(combined_df_for_plot) {
  ggplot(combined_df_for_plot, aes(x = avg_meth, fill = position)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Methylation Density by Gene Position (Center, Middle, Distal)",
      x = "Average Methylation (%)",
      fill = "Position"
    ) +
    theme_minimal() +
    # Optional: adjust x-axis if your methylation is 0-100%
    xlim(0, 100)
}
# Apply the plotting function (Assuming combined_df_for_plot is ready)
density_plot <- plot_methylation_density(combined_df_for_plot)
print(density_plot)
# step 4. plot as a geom_smooth on each palindromic region
# Assuming you have run the function above and stored the result here:
# final_combined_data <- combine_methylation_and_genes(plot_df_data, combined_data_with_position_data)
create_methylation_plot <- function(data_to_plot) {
  ggplot(data_to_plot, aes(x = position, y = avg_meth, color = palindrome, group = palindrome)) +
    geom_smooth(method = "loess", se = TRUE) + # Uses a LOESS smooth by default
    labs(
      title = "Smoothed Methylation Levels by Position and Palindrome",
      x = "Genomic Position Category",
      y = "Average Methylation (%)",
      color = "Palindrome ID"
    ) +
    theme_minimal() +
    # This ensures "distal", "center", "middle" appear in a logical order on the X-axis
    scale_x_discrete(limits = c("distal", "middle", "center")) 
}
# How to generate the plot:
my_plot <- create_methylation_plot(combined_df_for_plot)
print(my_plot)
