## ------ || INTRA-ARM IDENTITY ANALYSIS || ------
### ----- || Upload data || -----
identity_all = read_tsv("Work/Data/all_palindrome_identity.tsv")
identity_window = read_tsv("Work/Data/04_sliding_window_identity_P2.tsv")
identity_window = na.omit(identity_window)
# Rename coordinate columns
colnames(identity_window)[1:2] <- c("start", "end")
# Add chr column
identity_window <- data.frame(
  seqnames = "chrY",
  identity_window,
  stringsAsFactors = FALSE
)
identity_window$start <- as.numeric(identity_window$start)
identity_window$end   <- as.numeric(identity_window$end)
### ----- || Intra-arm sequence identity for P1-P8 - Global View || -----
ggplot(identity_window, aes(x = palindrome, y = percent_identity, color = palindrome)) +
  # Jitter spreads points horizontally so overlaps are visible
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  # Add a crossbar for the median
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  theme_minimal() +
  labs(
    title = "Intra-arm Sequence Identity by Palindrome (Jitter Plot)",
    x = "Palindrome Category",
    y = "Percent Identity (%)"
  )
### ----- || Intra-arm sequence identity for P1-P8 along the arms || -----
ggplot(identity_window, aes(x = end, y = percent_identity, group = palindrome)) +
  # Color points by individual
  geom_jitter(aes(color = sample), alpha = 0.5, size = 1.2, width = 0.2, height = 0) +
  # One column layout with shared axes
  facet_wrap(~palindrome, scales = "free", ncol = 1) + 
  # Remove the legend for individuals
  guides(color = "none") +
  # Styling
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "gray80", fill = NA)
  ) +
  labs(
    title = "Methylation vs Identity by Palindrome",
    subtitle = "Colored by Individual | Single Column Layout",
    x = "Sequence Identity (%)",
    y = "Average Methylation"
  )
####OR
# --- 1. Data Preparation & Mirroring ---
palindromes <- sort(unique(identity_window$palindrome))
samples     <- sort(unique(identity_window$sample))
# Empty container
binned_data <- data.frame()
for (p in palindromes) {
  for (s in samples) {
    
    sub_data <- identity_window[
      identity_window$palindrome == p &
        identity_window$sample == s, ]
    
    if (nrow(sub_data) > 0) {
      
      # Order by genomic coordinate
      sub_data <- sub_data[order(sub_data$end), ]
      
      # 30 equal-count bins (more robust than cut())
      sub_data$bin <- ceiling(
        rank(sub_data$end, ties.method = "first") /
          (nrow(sub_data) / 30)
      )
      
      sub_data$bin[sub_data$bin > 30] <- 30
      
      # Aggregate identity per bin
      agg <- aggregate(percent_identity ~ bin,
                       data = sub_data,
                       FUN = mean)
      
      agg$palindrome <- p
      agg$sample     <- s
      
      binned_data <- rbind(binned_data, agg)
    }
  }
}
par(mfrow = c(length(palindromes), 2),
    mar = c(1.5, 1.5, 1.5, 1),
    oma = c(4, 4, 2, 1),
    bg = "black")
sample_ids <- unique(binned_data$sample)
sample_cols <- setNames(
  rainbow(length(sample_ids)),
  sample_ids
)
for (p in palindromes) {
  
  p_data <- binned_data[binned_data$palindrome == p, ]
  
  y_range <- range(p_data$percent_identity, na.rm = TRUE)
  
  # ---------------- LEFT ARM ----------------
  plot(NA,
       xlim = c(1, 30),
       ylim = y_range,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       main = "",
       col = "white")
  
  grid(col = "gray30")
  
  for (s in sample_ids) {
    s_data <- p_data[p_data$sample == s, ]
    lines(s_data$bin,
          s_data$percent_identity,
          col = sample_cols[s],
          lwd = 1.5)
  }
  
  # ---------------- RIGHT ARM (mirrored) ----------------
  plot(NA,
       xlim = c(30, 1),  # reversed axis
       ylim = y_range,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       main = "",
       col = "white")
  
  grid(col = "gray30")
  
  for (s in sample_ids) {
    s_data <- p_data[p_data$sample == s, ]
    lines(s_data$bin,
          s_data$percent_identity,
          col = sample_cols[s],
          lwd = 1.5)
  }
}
mtext("Position (Bins)", side = 1, outer = TRUE, line = 2, col = "white")
mtext("Sequence Identity (%)", side = 2, outer = TRUE, line = 2, col = "white")
#OR
library(dplyr)
library(ggplot2)

# 1. Prepare mirrored data
identity_mirrored <- identity_window %>%
  group_by(palindrome) %>%
  mutate(bin_num = as.numeric(cut_interval(end, n = 50))) %>%
  ungroup()

identity_mirrored <- bind_rows(
  identity_mirrored %>%
    mutate(arm = "Arm A", display_bin = bin_num),
  identity_mirrored %>%
    mutate(arm = "Arm B", display_bin = 51 - bin_num)
) %>%
  mutate(
    palindrome = factor(
      palindrome,
      levels = sort(unique(palindrome), decreasing = TRUE)
    )
  )

# 2. Plotting mirrored boxplots
ggplot(identity_mirrored,
       aes(x = factor(display_bin), y = percent_identity, fill = palindrome)) +
  geom_boxplot(outlier.size = 0.1, alpha = 0.8, lwd = 0.1) +
  facet_grid(palindrome ~ arm, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = grid::unit(0.05, "lines"),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_viridis_d(option = "mako") +
  labs(
    title = "Mirrored Intra-arm Sequence Identity (50 Bins per Arm)",
    subtitle = "Left (Arm A): Start -> End | Right (Arm B): End -> Start",
    x = "Relative Position (Mirrored at center)",
    y = "Sequence Identity (%)"
  )
### ----- || Methylation heterogeneity and sequence identity || -----
#### merge all files in ldf_het list
ldf_het_all <- bind_rows(ldf_het, .id = "Individual_ID")
#### since there is no P2, add it only with "red","P2-spacer","red" occurrences
#### and compute met. het. mean
#### then filter only for P1-P8
ldf_het_all_P <- ldf_het_all %>%
  # 1. Name "red","P2-spacer","red" schemes as P2
  mutate(V4 = ifelse(V4 %in% c("red", "P2-spacer","red"), 
                                 "P2", 
                                 V4)) %>%
  # 2. Group per sequence class and individual
  group_by(Individual_ID, V4) %>%
  # 3. Measure heterogeneity mean (es. 'heterogeneity_value')
  summarise(V10 = mean(V10, na.rm = TRUE), 
            .groups = "drop") %>%
  # 4. FIlter for P1-P8
  filter(V4 %in% paste0("P", 1:8))
#### merge identity infos with heterogeneity
#### 1. Prepare data for joining: rename columns for uniformity
# Rename in ldf_het_all_P: V4 to sequence_class, V10 to heterogeneity
ldf_het_all_P_clean <- ldf_het_all_P %>%
  rename(sequence_class = V4, 
         heterogeneity = V10) %>%
  # Clean up white spaces as done previously
  mutate(sequence_class = trimws(sequence_class))
#### Rename in identity_all: sample to Individual_ID, palindrome to sequence_class
identity_all <- identity_all %>%
rename(Individual_ID = sample,
       sequence_class = palindrome) %>%
  mutate(sequence_class = trimws(sequence_class))
#### 2. Join the two data frames
# Use inner_join to keep only rows that have a match in both datasets
merged_data <- inner_join(ldf_het_all_P_clean, identity_all, 
                          by = c("Individual_ID", "sequence_class"))
#### 3. Create the requested plot
# Assuming 'merged_data' is clean and contains 'sequence_class', 
# 'percent_identity', and 'heterogeneity'
# Ensure the palindromes are treated as factors in the correct P1-P8 order
merged_data_ordered <- merged_data %>%
  mutate(sequence_class = factor(sequence_class, levels = paste0("P", 1:8)))
# --- Plot 1: Lateral Boxplot for Identities (Left Panel) ---
plot_identity_lateral <- ggplot(merged_data_ordered, 
                                aes(x = percent_identity, y = sequence_class)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  labs(title = "Percent Identity per Palindrome",
       x = "Percent Identity (%)",
       y = NULL) + 
  theme_minimal() + 
  theme(axis.text.y = element_text(hjust = 0))
# --- Plot 2: Ridge Density Plot for Heterogeneity (Right Panel) ---
plot_het_density <- ggplot(merged_data_ordered, 
                           aes(x = heterogeneity, y = sequence_class)) +
  geom_density_ridges_gradient(aes(fill = after_stat(x)), scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Het. Value") + 
  labs(title = "Heterogeneity Distribution",
       x = "Heterogeneity (met.het)",
       y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
# --- Combine Plots using patchwork ---
combined_plot_lateral <- plot_identity_lateral | plot_het_density
print(combined_plot_lateral)
### ----- || Methylation and sequence identity x palindrome || -----
# 1. Rebin methylation data
plot_df_rebinned <- plot_df %>%
  dplyr::group_by(individual, palindrome, arm) %>%
  dplyr::mutate(
    bin_id = as.integer(cut(start, breaks = 10, labels = 1:10))
  ) %>%
  dplyr::group_by(individual, palindrome, arm, bin_id) %>%
  dplyr::summarise(
    avg_meth = mean(avg_meth, na.rm = TRUE),
    .groups = "drop"
  )
# 2. Bin identity_window
binned_identity <- identity_window %>%
  dplyr::group_by(sample, palindrome) %>%
  dplyr::mutate(
    base_bin = as.integer(cut(start, breaks = 10, labels = 1:10))
  ) %>%
  dplyr::group_by(sample, palindrome, base_bin) %>%
  dplyr::summarise(
    avg_identity = mean(percent_identity, na.rm = TRUE),
    .groups = "drop"
  )
# 3. Create mirrored mappings
identity_left <- binned_identity %>%
  dplyr::mutate(arm = "left", bin_id = base_bin)
identity_right <- binned_identity %>%
  dplyr::mutate(arm = "right", bin_id = 11 - base_bin)
# 4. Combine and rename
identity_mirrored <- dplyr::bind_rows(identity_left, identity_right) %>%
  dplyr::rename(individual = sample)
# 5. Final match and plot
final_plot_data <- plot_df_rebinned %>%
  dplyr::inner_join(
    identity_mirrored,
    by = c("individual", "palindrome", "arm", "bin_id")
  )
ggplot(final_plot_data, aes(x = avg_identity, y = avg_meth, group = palindrome)) +
  # Color points by individual
  geom_jitter(aes(color = individual), alpha = 0.5, size = 1.2, width = 0.2, height = 0) +
  # One column layout with shared axes
  facet_wrap(~palindrome, scales = "fixed", ncol = 1) + 
  # Remove the legend for individuals
  guides(color = "none") +
  # Styling
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "gray80", fill = NA)
  ) +
  labs(
    title = "Methylation vs Identity by Palindrome",
    subtitle = "Colored by Individual | Single Column Layout",
    x = "Sequence Identity (%)",
    y = "Average Methylation"
  )
### ----- || Methylation and sequence identity x palindrome - <90, 90-99, 100% regions || -----
raw_plot_df_filtered <- raw_plot_df %>%
  filter(class %in% unique(data_arms$palindrome))
# Convert filtered data frame to a list by individual_id
raw_plot_df_list <- split(raw_plot_df_filtered, raw_plot_df_filtered$individual_id)
# Convert identity_window data frame to a list by the 'sample' column
identity_window_list <- split(identity_window, identity_window$sample)
# Find the intersection of names (individual IDs)
common_ids <- intersect(names(raw_plot_df_list), names(identity_window_list))
# Filter both lists to keep only the common IDs
raw_plot_df_list_aligned <- raw_plot_df_list[common_ids]
identity_window_list_aligned <- identity_window_list[common_ids]
# Also ensure the annotation list has corresponding names for accurate matching in Step 3
# Rename the 'chr' column to 'seqnames' in every data frame within the list
# Safely rename 'chr' to 'seqnames' only if 'chr' exists
ldf_ann <- purrr::map(ldf_ann, ~ {
  if ("chr" %in% names(.x)) rename(.x, seqnames = chr) else .x
})

ldf_ann_aligned <- ldf_ann[common_ids]
calculate_relative_pos_fast <- function(df_list, ann_list) {
  map2(df_list, ann_list, function(df, ann) {
    
    # Standardize column names for valr
    df_valr <- df %>% rename(chrom = seqnames)
    ann_valr <- ann %>% rename(chrom = seqnames)
    
    # Perform intersection
    res <- bed_intersect(df_valr, ann_valr, suffix = c("", ".ann")) %>%
      mutate(relative_pos = (end - start.ann) / (end.ann - start.ann))
    
    # Return to original column naming style for the next steps
    res %>% 
      rename(seqnames = chrom) %>%
      select(any_of(c(names(df), "relative_pos")))
  })
}
# Apply the function to the aligned lists
raw_plot_df_rel_list <- calculate_relative_pos_fast(raw_plot_df_list_aligned, ldf_ann_aligned)
# repeat for identity_window
calculate_relative_pos_identity_only <- function(df_list) {
  map(df_list, function(df) {
    
    # 1. Ensure column name consistency
    if ("sample" %in% names(df)) {
      df <- df %>% rename(individual_id = sample)
    }
    
    # 2. Find the max end value (the end of the palindrome in this reset scale)
    # We use max(df$end) because 0 is the start.
    max_coord <- max(df$end, na.rm = TRUE)
    
    # 3. Calculate relative position
    df <- df %>%
      mutate(
        relative_pos = end / max_coord
      )
    
    return(df)
  })
}
# Apply this to the identity list (no ann_list needed here)
identity_window_rel_list <- calculate_relative_pos_identity_only(identity_window_list_aligned)
# Define the breaks (0, 0.05, 0.10 ... 1.0)
bins <- seq(0, 1, length.out = 21)
create_bins <- function(df_list) {
  map(df_list, function(df) {
    # Ensure individual_id is consistent
    if ("sample" %in% names(df)) df <- rename(df, individual_id = sample)
    
    df %>%
      mutate(
        # cut() assigns each relative_pos to a bin (e.g., "(0.05, 0.1]")
        bin = cut(relative_pos, breaks = bins, include.lowest = TRUE, labels = FALSE)
      ) %>%
      filter(!is.na(bin)) # Remove rows that fell outside [0,1]
  })
}
# Apply binning to both lists
meth_binned_list <- create_bins(raw_plot_df_rel_list)
iden_binned_list <- create_bins(identity_window_rel_list)
# Convert lists to large data tables
dt_meth <- as.data.table(bind_rows(meth_binned_list))
dt_iden <- as.data.table(bind_rows(iden_binned_list))
colnames(dt_iden)[8] = "class"
# Merge based on your requested keys: individual, palindrome (class), and bin
# We use 'class' here assuming it contains the palindrome ID
# 1. Collapse identity data so there is only ONE row per (individual, class, bin)
# We take the mean of percent_identity for each bin
dt_iden_avg <- dt_iden[, .(
  percent_identity = mean(percent_identity, na.rm = TRUE)
), by = .(individual_id, class, bin)]

# 2. Now perform the merge
# This will be much smaller and faster because each methylation row 
# now matches exactly ONE identity value
final_combined_df <- merge(
  dt_meth, 
  dt_iden_avg, 
  by = c("individual_id", "class", "bin"),
  suffixes = c(".meth", ".identity"),
  all.x = TRUE # Keep all methylation data even if no identity window matches
)
# Prepare categories
plot_df <- final_combined_df %>%
  dplyr::mutate(
    identity_group = dplyr::case_when(
      percent_identity == 100 ~ "100%",
      percent_identity >= 90 & percent_identity < 100 ~ "90-99%",
      percent_identity < 90 ~ "<90%",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::mutate(
    identity_group = factor(
      identity_group,
      levels = c("<90%", "90-99%", "100%")
    )
  )
# Plot: 1 plot per Palindrome (class), stacked in one column
ggplot(plot_df, aes(x = meth, fill = identity_group)) +
  geom_density(alpha = 0.4) +
  # 'class' is the palindrome ID column
  facet_grid(class ~ ., scales = "free_y") + 
  theme_minimal() +
  labs(
    title = "Methylation Density by Identity Bins",
    x = "Methylation Value",
    y = "Density",
    fill = "Identity Tier"
  ) +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "top"
  )
### ----- || Barplot with sequence identity per palindrome || -----
# 1. Categorize identity values
binned_counts <- final_plot_data %>%
  dplyr::mutate(
    identity_group = dplyr::case_when(
      avg_identity == 100 ~ "100%",
      avg_identity >= 90 & avg_identity < 100 ~ "90-99%",
      avg_identity < 90 ~ "< 90%",
      TRUE ~ NA_character_
    )
  ) %>%
  # Ensure the categories stay in the correct order on the plot
  dplyr::mutate(
    identity_group = factor(
      identity_group,
      levels = c("< 90%", "90-99%", "100%")
    )
  )
# 2. Create the barplot
ggplot(binned_counts, aes(x = palindrome, fill = identity_group)) +
  # CHANGE: position = "fill" normalizes bars to 100%
  geom_bar(position = "fill", color = "black") +
  
  # Format Y axis as 0-100%
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = c("< 90%" = "#e41a1c", 
                               "90-99%" = "#377eb8", 
                               "100%" = "#4daf4a")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Relative Proportion of Identity Classes",
    x = "Palindrome",
    y = "Percentage (%)",
    fill = "Identity Group"
  )
### ----- || Intra-arm sequence identity on regions covered by genes vs. devoid of genes || -----
# Filter raw_plot_df for sequence_class P1-P8
# Assuming the column name in the dataframe is 'class' based on your image
raw_plot_df_filtered <- raw_plot_df %>%
  filter(class %in% unique(data_arms$sequence_class))
# Convert filtered data frame to a list by individual_id
raw_plot_df_list <- split(raw_plot_df_filtered, raw_plot_df_filtered$individual_id)
# Convert identity_window data frame to a list by the 'sample' column
identity_window_list <- split(identity_window, identity_window$sample)
# Find the intersection of names (individual IDs)
common_ids <- intersect(names(raw_plot_df_list), names(identity_window_list))
# Filter both lists to keep only the common IDs
raw_plot_df_list_aligned <- raw_plot_df_list[common_ids]
identity_window_list_aligned <- identity_window_list[common_ids]
# Also ensure the annotation list has corresponding names for accurate matching in Step 3
# Rename the 'chr' column to 'seqnames' in every data frame within the list
# Safely rename 'chr' to 'seqnames' only if 'chr' exists
ldf_ann <- purrr::map(ldf_ann, ~ {
  if ("chr" %in% names(.x)) rename(.x, seqnames = chr) else .x
})

ldf_ann_aligned <- ldf_ann[common_ids]
calculate_relative_pos_fast <- function(df_list, ann_list) {
  map2(df_list, ann_list, function(df, ann) {
    
    # Standardize column names for valr
    df_valr <- df %>% rename(chrom = seqnames)
    ann_valr <- ann %>% rename(chrom = seqnames)
    
    # Perform intersection
    res <- bed_intersect(df_valr, ann_valr, suffix = c("", ".ann")) %>%
      mutate(relative_pos = (end - start.ann) / (end.ann - start.ann))
    
    # Return to original column naming style for the next steps
    res %>% 
      rename(seqnames = chrom) %>%
      select(any_of(c(names(df), "relative_pos")))
  })
}
# Apply the function to the aligned lists
raw_plot_df_rel_list <- calculate_relative_pos_fast(raw_plot_df_list_aligned, ldf_ann_aligned)
# repeat for identity_window
calculate_relative_pos_identity_only <- function(df_list) {
  map(df_list, function(df) {
    
    # 1. Ensure column name consistency
    if ("sample" %in% names(df)) {
      df <- df %>% rename(individual_id = sample)
    }
    
    # 2. Find the max end value (the end of the palindrome in this reset scale)
    # We use max(df$end) because 0 is the start.
    max_coord <- max(df$end, na.rm = TRUE)
    
    # 3. Calculate relative position
    df <- df %>%
      mutate(
        relative_pos = end / max_coord
      )
    
    return(df)
  })
}
# Apply this to the identity list (no ann_list needed here)
identity_window_rel_list <- calculate_relative_pos_identity_only(identity_window_list_aligned)
# Define the breaks (0, 0.05, 0.10 ... 1.0)
bins <- seq(0, 1, length.out = 21)
create_bins <- function(df_list) {
  map(df_list, function(df) {
    # Ensure individual_id is consistent
    if ("sample" %in% names(df)) df <- rename(df, individual_id = sample)
    
    df %>%
      mutate(
        # cut() assigns each relative_pos to a bin (e.g., "(0.05, 0.1]")
        bin = cut(relative_pos, breaks = bins, include.lowest = TRUE, labels = FALSE)
      ) %>%
      filter(!is.na(bin)) # Remove rows that fell outside [0,1]
  })
}
# Apply binning to both lists
meth_binned_list <- create_bins(raw_plot_df_rel_list)
iden_binned_list <- create_bins(identity_window_rel_list)
# Convert lists to large data tables
dt_meth <- as.data.table(bind_rows(meth_binned_list))
dt_iden <- as.data.table(bind_rows(iden_binned_list))
colnames(dt_iden)[8] = "class"
# Merge based on your requested keys: individual, palindrome (class), and bin
# We use 'class' here assuming it contains the palindrome ID
# 1. Collapse identity data so there is only ONE row per (individual, class, bin)
# We take the mean of percent_identity for each bin
dt_iden_avg <- dt_iden[, .(
  percent_identity = mean(percent_identity, na.rm = TRUE)
), by = .(individual_id, class, bin)]

# 2. Now perform the merge
# This will be much smaller and faster because each methylation row 
# now matches exactly ONE identity value
final_combined_df <- merge(
  dt_meth, 
  dt_iden_avg, 
  by = c("individual_id", "class", "bin"),
  suffixes = c(".meth", ".identity"),
  all.x = TRUE # Keep all methylation data even if no identity window matches
)
# 3. Create a violin plot using ggplot2
ggplot(final_combined_df, aes(x = overlap_group, y = percent_identity, fill = overlap_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) + # trim=FALSE extends violins to extreme data points
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) + # Optional: add a mini boxplot inside
  labs(
    title = "Methylation Levels: Gene-Covered vs. Uncovered Palindromic Regions (Violin Plot)",
    x = "Genomic Region Status",
    y = "Average Methylation (%)",
    fill = "Genomic Region Status"
  ) +
  theme_minimal() + ylim(99.8,100) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme(legend.position = "none") # Hide legend if x-axis labels are clear enough
wilcox.test(final_combined_df$percent_identity[which(final_combined_df$overlaps_gene==TRUE)],final_combined_df$percent_identity[which(final_combined_df$overlaps_gene==FALSE)])
### ------- || Create plot to show % identity on ampliconic genes based on their color block location || --------
# Ensure 'class' exists in data_arms for the join
# If it's called 'sequence_class' in data_arms, we rename it
if ("sequence_class" %in% names(data_arms)) {
  data_arms <- data_arms %>% rename(class = sequence_class)
}
# Merge the ordering and color information
plot_df_ordered <- final_combined_df %>%
  left_join(data_arms %>% 
              select(class, blocks_order), 
            by = "class") %>%
  # Reorder the 'class' factor based on the color_block_order
  mutate(class = reorder(class, blocks_order))
ggplot(plot_df_ordered, aes(x = blocks_order, y = percent_identity, group = geneName)) +
  # Area plot showing the identity profile
  geom_area(fill = "steelblue", alpha = 0.6) +
  # Facet by the ordered class (1 per row)
  facet_grid(class ~ ., scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Sequence Identity Profile across Palindromes",
    subtitle = "Ordered by color_block_order",
    x = "Relative Position (Bins 1-20)",
    y = "Percent Identity"
  ) +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.grid.minor = element_blank()
  )
