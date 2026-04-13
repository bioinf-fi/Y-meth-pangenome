## --------- || UPDATED: METHYLATION ACROSS PALINDROME COPIES || -----------
### ----- || Methylation on distinct copies || -----
### 1. Prepare Reference Table (from your Part 1 data_arms)
# Use dplyr:: explicitly to avoid conflicts with plyr
expected_counts <- data_arms %>%
  dplyr::group_by(palindrome) %>%
  dplyr::summarise(expected_n = dplyr::n(), .groups = 'drop') %>%
  dplyr::rename(palindrome_id = palindrome)
### 2. Calculate Normalized Occurrences (using combined_ann_ldf)
palindrome_counts <- purrr::map_dfr(names(combined_ann_ldf), function(indiv_name) {
  df <- combined_ann_ldf[[indiv_name]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    dplyr::group_by(parent_palindrome) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate(individual = indiv_name)
})

normalized_palindromes <- palindrome_counts %>%
  dplyr::left_join(expected_counts, by = c("parent_palindrome" = "palindrome_id")) %>%
  dplyr::mutate(normalized_count = count / expected_n)

### 3. Assign Copy Numbers to the Annotation List
combined_ann_with_copy <- lapply(names(combined_ann_ldf), function(indiv_name) {
  df <- combined_ann_ldf[[indiv_name]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    dplyr::left_join(expected_counts, by = c("parent_palindrome" = "palindrome_id")) %>%
    dplyr::group_by(parent_palindrome, color_block, arm) %>%
    dplyr::mutate(occurrence_id = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(copy = dplyr::case_when(
      occurrence_id <= 1 ~ 1, 
      TRUE ~ occurrence_id     
    )) %>%
    dplyr::select(-expected_n, -occurrence_id) %>%
    dplyr::mutate(individual = indiv_name)
})
names(combined_ann_with_copy) <- names(combined_ann_ldf)
### 4. Merge Raw Methylation with updated Annotation (FIXED FOR SEQLEVELS)
merge_raw_methylation_updated <- function(meth_df, ann_df) {
  # 1. Standardize column names
  colnames(meth_df)[1:4] <- c("chr", "start", "end", "meth")
  
  # 2. FORCE CHROMOSOME CONSISTENCY
  # Strip "chr" and add it back to ensure both are "chrY", "chr1", etc.
  meth_df$chr <- paste0("chr", gsub("chr", "", meth_df$chr))
  ann_df$chr  <- paste0("chr", gsub("chr", "", ann_df$chr))
  
  # 3. Create GRanges
  meth_gr <- GenomicRanges::GRanges(
    seqnames = meth_df$chr, 
    ranges = IRanges::IRanges(meth_df$start, meth_df$end), 
    meth = meth_df$meth
  )
  
  ann_gr <- GenomicRanges::GRanges(
    seqnames = ann_df$chr, 
    ranges = IRanges::IRanges(ann_df$start, ann_df$end),
    palindrome = ann_df$parent_palindrome,
    arm = ann_df$arm,
    copy = ann_df$copy
  )
  
  # 4. Find overlaps
  hits <- suppressWarnings(GenomicRanges::findOverlaps(ann_gr, meth_gr))
  
  if (length(hits) > 0) {
    res <- as.data.frame(hits)
    
    # Extract coordinates and values from the methylation GRanges (subject)
    res$chr   <- as.character(GenomicRanges::seqnames(meth_gr)[res$subjectHits])
    res$start <- GenomicRanges::start(meth_gr)[res$subjectHits]
    res$end   <- GenomicRanges::end(meth_gr)[res$subjectHits]
    res$meth  <- GenomicRanges::mcols(meth_gr)$meth[res$subjectHits]
    
    # Extract metadata from the annotation GRanges (query)
    res$palindrome <- GenomicRanges::mcols(ann_gr)$palindrome[res$queryHits]
    res$arm        <- GenomicRanges::mcols(ann_gr)$arm[res$queryHits]
    res$copy       <- GenomicRanges::mcols(ann_gr)$copy[res$queryHits]
    
    # Return the dataframe with the new coordinate columns included
    return(res %>% dplyr::select(chr, start, end, meth, palindrome, arm, copy))
  } else {
    return(data.frame())
  }
  
}
### 5. Execute Merge
raw_data_list <- list()
for (indiv in names(ldf)) {
  if (indiv %in% names(combined_ann_with_copy) && !is.null(combined_ann_with_copy[[indiv]])) {
    raw_data_list[[indiv]] <- merge_raw_methylation_updated(ldf[[indiv]], combined_ann_with_copy[[indiv]]) %>%
      dplyr::mutate(individual = indiv)
  }
}
raw_plot_df <- dplyr::bind_rows(raw_data_list) 
### 5. Visualization
#### 1) Boxplot of methylation on distinct copies of P2
##### Filter the raw data for P2 and ensure 'copy' is a discrete factor
p2_methylation_data <- raw_plot_df %>%
  filter(palindrome == "P2") %>%
  mutate(copy = as.factor(copy)) %>%
  # Optionally remove NA methylation values if any exist
  filter(!is.na(meth))
##### Generate the Boxplot
ggplot(p2_methylation_data, aes(x = copy, y = meth, fill = copy)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    title = "Methylation Distribution for P2 by Copy Number",
    x = "P2 Copy Number",
    y = "Methylation Level (0 to 1)"
  ) +
  theme_classic() +
  theme(legend.position = "none")
# Optional: Run statistical tests as originally intended
if (n_distinct(p2_methylation_data$copy) > 1) {
  print("Wilcoxon test results comparing copies:")
  # Note: These tests will only run if both copies exist in the filtered data
  # tryCatch handles potential errors if a copy group is empty
  tryCatch({
    w12 <- wilcox.test(meth ~ copy, data = filter(p2_methylation_data, copy %in% c(1, 2)))
    print(paste("P-value Copy 1 vs 2:", w12$p.value))
    w13 <- wilcox.test(meth ~ copy, data = filter(p2_methylation_data, copy %in% c(1, 3)))
    print(paste("P-value Copy 1 vs 3:", w13$p.value))
  }, error = function(e) {
    print("Could not perform all Wilcoxon tests (likely missing data for one copy number)")
  })
}
##### Perform Kruskal-Wallis test (non-parametric ANOVA)
kruskal_p2 <- kruskal.test(meth ~ copy, data = p2_methylation_data)
print("Kruskal-Wallis Test Results for P2 Copies:")
print(kruskal_p2)
### ----- || Methylation on distinct P2 copies across left and right arms || -----
# Assuming 'plot_df' is the data frame with binned methylation for all individuals
# And 'combined_ann_with_copy' is the list of data frames with correct 'copy' numbers
# Function to perform the spatial merge for a single individual (FIXED COLUMNS AND CBIND)
merge_bins_with_annotations <- function(bins_df, ann_df) {
  if (is.null(bins_df) || is.null(ann_df) || nrow(bins_df) == 0 || nrow(ann_df) == 0) {
    return(data.frame())
  }
  
  # Remove potentially conflicting columns from bins_df before conversion to avoid 'arm.1' issue
  bins_df_clean <- bins_df %>% dplyr::select(-any_of(c("copy", "arm", "palindrome")))
  
  # Prepare GRanges objects for spatial overlap
  bins_gr <- makeGRangesFromDataFrame(bins_df_clean, keep.extra.columns = TRUE)
  
  # Ensure annotation dataframe has necessary columns for GRanges conversion
  ann_df$copy <- ann_df$copy # Make sure 'copy' is present
  ann_df$arm <- ann_df$arm   # Make sure 'arm' is present
  ann_df$parent_palindrome <- ann_df$parent_palindrome # Make sure 'parent_palindrome' is present
  ann_gr <- makeGRangesFromDataFrame(ann_df, keep.extra.columns = TRUE)
  
  # Find overlaps (each bin should overlap exactly one annotation region)
  hits <- suppressWarnings(findOverlaps(bins_gr, ann_gr, type = "within"))
  
  if (length(hits) > 0) {
    # Extract annotations (copy, arm, parent_palindrome) that match the bins
    matched_ann <- as.data.frame(mcols(ann_gr[subjectHits(hits)])) %>% 
      dplyr::select(copy, arm, parent_palindrome)
    
    # Attach these annotations back to the original bins data using base::cbind
    results_df <- as.data.frame(bins_gr[queryHits(hits)]) %>%
      dplyr::select(-width, -strand) %>%
      dplyr::rename(chr = seqnames) %>%
      base::cbind(matched_ann)
    
    return(results_df)
  } else {
    return(data.frame())
  }
}
# Iterate through individuals in plot_df and merge them with the list data
final_plot_data_list <- lapply(unique(plot_df$individual), function(indiv_name) {
  # Filter plot_df for the current individual's bins
  indiv_bins <- plot_df %>% dplyr::filter(individual == indiv_name)
  # Get the corresponding annotation data from the list
  indiv_ann <- combined_ann_with_copy[[indiv_name]]
  
  if (nrow(indiv_bins) > 0 && !is.null(indiv_ann)) {
    message("Merging annotations for individual: ", indiv_name)
    merged_indiv_df <- merge_bins_with_annotations(indiv_bins, indiv_ann)
    # Add the individual column back as it might have been lost in GRanges conversion
    if(nrow(merged_indiv_df) > 0) {
      merged_indiv_df$individual <- indiv_name
    }
    return(merged_indiv_df)
  } else {
    message("Skipping individual: ", indiv_name, " due to missing data.")
    return(NULL)
  }
})
# Combine all individuals back into a single, final data frame
final_plot_df <- dplyr::bind_rows(final_plot_data_list[!sapply(final_plot_data_list, is.null)]) %>%
  # Rename 'parent_palindrome' back to 'palindrome' for consistency
  dplyr::rename(palindrome = parent_palindrome) %>%
  # Filter only for left and right arms as requested
  dplyr::filter(arm %in% c("left", "right")) %>%
  # Ensure palindrome column is a factor for correct plotting order
  dplyr::mutate(palindrome = factor(palindrome, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")))
# Assuming 'final_plot_df' is available from the previous step and libraries are loaded
# Data cleaning and preparation for plotting
# The final_plot_df is already mostly clean, just add factors and filter for completeness
plot_df_clean <- final_plot_df %>%
  filter(!is.na(avg_meth), !is.na(bin_id), arm %in% c("left", "right")) %>%
  mutate(
    palindrome = factor(palindrome, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")),
    arm = factor(arm, levels = c("left", "right")),
    copy = factor(copy) # Treat copy number as a factor for coloring
  ) %>%
  # Group by all relevant factors to ensure enough data points per trend line
  group_by(palindrome, arm, individual, copy) %>%
  filter(n() >= 10) %>% # Ensure at least 10 bins are present for a smooth line
  ungroup()
# Plot 1: Smoothed Trends Across All Palindromes (as originally requested)
# This plot uses color to distinguish arms
plot_smooth_all_palindromes <- ggplot(plot_df_clean, aes(x = bin_id, y = avg_meth, color = arm, fill = arm)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2, se = TRUE) +
  facet_wrap(~palindrome, ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  labs(
    title = "Methylation Trends Across All Palindrome Arms (Left vs. Right)",
    subtitle = "Aggregated trends using GAM smoothing (95% CI)",
    x = "Bin (Proximal 1 -> Distal 10)", 
    y = "Mean Methylation (%)",
    color = "Arm", fill = "Arm"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
print(plot_smooth_all_palindromes)
# Plot 2: Detailed comparison of P2 Copies (focusing on your specific request)
# This plot uses color to distinguish P2 Copy 1 vs Copy 2
plot_p2_copies <- plot_df_clean %>%
  filter(palindrome == "P2") %>%
  ggplot(aes(x = bin_id, y = avg_meth, color = interaction(copy, arm), group = interaction(copy, arm))) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  labs(
    title = "Comparison of Methylation in P2 Copies (Copy 1 vs. Copy 2)",
    subtitle = "Trends across 20 bins (Left and Right arms combined for smoothness)",
    x = "Bin Number",
    y = "Mean Methylation (%)",
    color = "P2 Copy"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(plot_p2_copies)
### ----- || Example of methylation on P2 1st vs. 2nd occurrences and on P1 || -----
# Allow Gviz to use non-standard chromosome names (as a fallback/safety measure)
options(ucscChromosomeNames = FALSE) 
# Load Gviz and necessary libraries
# Assuming 'final_plot_df' is the data frame from previous steps
# 1. Filter the data for P1 and P2 regions of interest, specifically for individual HG00096 on chrY
indiv_example <- "HG00096" 
plot_data_p1_p2 <- raw_plot_df %>%
  filter(individual == indiv_example, 
         palindrome %in% c("P1", "P2"), 
         copy %in% c(1, 2),
         chr == "chrY") %>% 
  arrange(start)
if(nrow(plot_data_p1_p2) == 0) {
  stop(paste("No P1 or P2 data found on chrY in final_plot_df for individual", indiv_example))
}
# 2. Create separate DataTracks for each region
track_p1 <- DataTrack(
  range = plot_data_p1_p2 %>% filter(palindrome == "P1") %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE),
  name = "P1 Methylation",
  type = "a", col = "blue", ylim = c(0, 100)
)
track_p2_c1 <- DataTrack(
  range = plot_data_p1_p2 %>% filter(palindrome == "P2", copy == 1) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE),
  name = "P2 Copy 1",
  type = "a", col = "red", ylim = c(0, 100)
)
track_p2_c2 <- DataTrack(
  range = plot_data_p1_p2 %>% filter(palindrome == "P2", copy == 2) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE),
  name = "P2 Copy 2",
  type = "a", col = "darkgreen", ylim = c(0, 100)
)
# 3. Overlay the DataTracks into a single panel
overlay_track <- OverlayTrack(
  trackList = list(track_p1, track_p2_c1, track_p2_c2),
  name = "Methylation Trends"
)
# 4. Create an AnnotationTrack to show the structure (optional, can be combined with overlay)
annotation_track <- AnnotationTrack(
  range = plot_data_p1_p2 %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE),
  name = "Regions",
  stacking = "dense", showID = TRUE, col = "black",
  feature = paste0(plot_data_p1_p2$palindrome, "_C", plot_data_p1_p2$copy) # Use feature mapping for IDs
)
# 5. Create a GenomeAxisTrack
axis_track <- GenomeAxisTrack()
# 6. Plot the tracks together
plotTracks(list(axis_track, annotation_track, overlay_track), 
           from = min(plot_data_p1_p2$start), 
           to = max(plot_data_p1_p2$end),
           title.width = 1.5,
           main = paste("Methylation in P1 and P2 Regions for Individual", indiv_example, "on chrY"))
### ----- || Barplot with n.observed/expected occurrences || -----
# Assumes 'data_arms' and 'combined_ann_ldf' are available from previous steps
# 1. Calculate Expected counts *per palindrome* from the data_arms reference
expected_ref_pal <- data_arms %>%
  dplyr::group_by(palindrome) %>%
  dplyr::summarise(total_exp_in_pal = dplyr::n(), .groups = 'drop') %>%
  dplyr::rename(parent_group = palindrome)
# 2. Count observed occurrences *per palindrome* from the overlap results (combined_ann_ldf)
observed_counts_pal <- purrr::map_dfr(names(combined_ann_ldf), function(indiv_name) {
  df <- combined_ann_ldf[[indiv_name]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    dplyr::group_by(parent_palindrome) %>%
    dplyr::summarise(total_obs_in_pal = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate(individual = indiv_name) %>%
    dplyr::rename(parent_group = parent_palindrome)
})
# 3. Calculate the Ratio and prepare the data frame for plotting
pal_ratio_df <- observed_counts_pal %>%
  dplyr::left_join(expected_ref_pal, by = "parent_group") %>%
  dplyr::mutate(ratio = total_obs_in_pal / total_exp_in_pal) %>%
  # Ensure P1-P8 order
  dplyr::mutate(parent_group = factor(parent_group, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")))
## Occurrences in P2 are 3
pal_ratio_df$total_exp_in_pal[which(pal_ratio_df$parent_group=="P2")] = 3
pal_ratio_df$ratio = pal_ratio_df$total_obs_in_pal/pal_ratio_df$total_exp_in_pal
# 4. Generate the Barplot with Jittered points
ggplot(pal_ratio_df, aes(x = parent_group, y = ratio, fill = parent_group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.7, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, na.rm = TRUE) +
  geom_jitter(aes(color = parent_group), width = 0.2, size = 1, alpha = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  theme_minimal() +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Palindrome Integrity: Observed / Expected Ratio (Grouped by Palindrome)",
    x = "Palindrome Group (P1-P8)",
    y = "Ratio (Total Observed Blocks / Total Expected Blocks)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
### ----- || P2 Methylation by Copy Number Boxplot || -----
# Filter the raw data for P2 and ensure 'copy' is a discrete factor
p2_methylation_data <- raw_plot_df %>%
  dplyr::filter(palindrome == "P2") %>%
  dplyr::mutate(copy = as.factor(copy)) %>%
  dplyr::filter(!is.na(meth)) # Remove NA methylation values
# Generate the Boxplot
ggplot(p2_methylation_data, aes(x = copy, y = meth, fill = copy)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    title = "Methylation Distribution for P2 by Copy Number",
    x = "P2 Copy Number",
    y = "Methylation Level (0 to 1)"
  ) +
  theme_classic() +
  theme(legend.position = "none")
