## ------ || METHYLATION ON DIFFERENT AZFc ARCHITECTURES || ------
### ------ || Methylation of P1, P2, and P3 in individuals with different architectures || ------
### Upload data
meth_AZFc_arch_all = read_tsv("Work/Data/methylation_region_summary_pal_with_AZFc.tsv")
AZFc_arch = read.delim("Work/Data/UniqueAZFc_Structures.tsv", header = F, stringsAsFactors = F)
### 1. Identify and create the P2 combined rows
# Assuming your dataframe is sorted by genomic position
meth_AZFc_arch_all_updated <- meth_AZFc_arch_all %>%
  group_by(sample, architecture) %>%
  mutate(
    # Logic: If current is P2-spacer, and neighbors are "red" (or the appropriate amplicon)
    is_p2_scheme = region_name == "P2-spacer" & 
      lag(region_name) == "red" & 
      lead(region_name) == "red"
  )

### Extract the P2 combined data
p2_combined <- meth_AZFc_arch_all_updated %>%
  filter(is_p2_scheme | lag(is_p2_scheme) | lead(is_p2_scheme)) %>%
  group_by(sample, architecture, group_id = cumsum(is_p2_scheme | (lag(is_p2_scheme) & !is_p2_scheme))) %>%
  summarise(
    region_name = "P2",
    region_type = "combined_palindrome",
    # Weighted average of methylation based on CpG count
    mean_methyl = sum(mean_methyl * n_cpg) / sum(n_cpg),
    n_cpg = sum(n_cpg),
    .groups = "drop"
  )

### 2. Combine back, Filter, and Plot
final_plot_df <- bind_rows(meth_AZFc_arch_all, p2_combined) %>%
  # Filter for your specific targets
  filter(region_name %in% c("P1", "P2", "P3"))

### 3. Create the Violin Plot
ggplot(final_plot_df, aes(x = region_name, y = mean_methyl, fill = region_name)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  facet_wrap(~ architecture) +
  theme_minimal() +
  labs(
    title = "Methylation Values for P1, P2, and P3",
    subtitle = "Aggregated across all architectures",
    x = "Region",
    y = "Mean Methylation (%)",
    fill = "Region"
  )
### ------ || AZFc methylation in individuals with different architectures || ------
# 1. Identify "Normal" samples (those in plot_df but not in the unique list)
normal_data <- plot_df %>%
  filter(!(individual %in% unique(meth_AZFc_arch_all$sample))) %>%
  mutate(architecture = "Normal") %>%
  # Rename columns to match the 'meth_AZFc_arch_all' structure
  select(architecture, mean_methyl = avg_meth) 
# 2. Extract and simplify the unique architecture data
unique_data_to_bind <- meth_AZFc_arch_all %>%
  select(architecture, mean_methyl)
# 3. Combine both and set the order (Normal first)
plot_combined_arch <- bind_rows(normal_data, unique_data_to_bind) %>%
  mutate(architecture = factor(architecture, 
                               levels = c("Normal", sort(unique(unique_data_to_bind$architecture)))))
# 4. Create the updated plot
ggplot(plot_combined_arch, aes(x = architecture, y = mean_methyl, fill = architecture)) +
  geom_violin(alpha = 0.6, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.5) +
  # Using a distinct color for Normal vs the rest if you like
  scale_fill_manual(values = c("Normal" = "grey70", 
                               setNames(scales::hue_pal()(length(unique(unique_data_to_bind$architecture))), 
                                        sort(unique(unique_data_to_bind$architecture))))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Methylation Distribution: Normal vs. Unique Architectures",
    subtitle = "Global comparison across all genomic structures",
    x = "Architecture Type",
    y = "Mean Methylation (%)"
  )
### ------ || N. of individuals with different architectures || ------
# 1. Summarize the data to count unique samples per architecture
architecture_counts <- meth_AZFc_arch_all %>%
  group_by(architecture) %>%
  summarise(n_individuals = n_distinct(sample))
# 2. Create the bar plot
ggplot(architecture_counts, aes(x = architecture, y = n_individuals, fill = architecture)) +
  geom_col() +
  # Add text labels on top of bars for clarity
  geom_text(aes(label = n_individuals), vjust = -0.5, size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Number of Individuals per Genomic Architecture",
    x = "Architecture",
    y = "Count of Unique Individuals"
  )
### ------ || Methylation on AZFc in individuals with Architecture:6 || ------
# --- 1. Filter individuals with Architecture:6 ---
# Clean IDs: remove _chrY to match ldf/ldf_ann names
arch6_ids <- AZFc_arch %>%
  filter(V1 == "Architecture:6") %>%
  mutate(sample_clean = str_remove(V2, "_chrY")) %>%
  pull(sample_clean)
ldf_filt <- ldf[names(ldf) %in% arch6_ids]
ldf_ann_filt <- ldf_ann[names(ldf_ann) %in% arch6_ids]
# --- 2. Rename blocks by occurrence (Sample-specific) ---
# Define the colors/names we are looking for in sequence_class
# The official Architecture:6 order provided
arch6_template <- c("blue", "yellow", "green", "red", "red", "green", "yellow", "blue", 
                          "green", "red", "red", "gray", "blue", "teal", "P3-spacer", "teal", "blue", "gray")
process_ann <- function(df) {
  df %>%
    # Filter for blocks that are part of the AZFc sequence
    filter(sequence_class %in% arch6_template) %>%
    arrange(as.numeric(start)) %>%
    # Generate labels like blue1, blue2 based on their appearance in this sample
    mutate(color_blocks = paste0(sequence_class, 
                                 ave(as.character(sequence_class), 
                                     as.character(sequence_class), 
                                     FUN = seq_along))) %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    select(start, end, color_blocks)
}
ldf_ann_renamed <- lapply(ldf_ann_filt, process_ann)
# --- 3. Overlap and 5-bin division ---
# Function to join methylation with blocks and create bins
bin_meth_gr <- function(m_df, a_df) {
  if (is.null(a_df) || nrow(a_df) == 0) return(NULL)
  
  # 1. Create GenomicRanges (using your specific column names)
  # We use the 'start' column from ldf as the position
  m_gr <- GRanges("chrY", IRanges(start = as.numeric(m_df$start), end = as.numeric(m_df$end)))
  a_gr <- GRanges("chrY", IRanges(start = a_df$start, end = a_df$end))
  mcols(a_gr)$block_id <- a_df$color_blocks
  
  # 2. Fast Overlap
  hits <- findOverlaps(m_gr, a_gr)
  
  # 3. Build result table
  res <- data.frame(
    pos = as.numeric(m_df$start[queryHits(hits)]),
    meth = as.numeric(m_df$meth[queryHits(hits)]),
    block = mcols(a_gr)$block_id[subjectHits(hits)],
    # Get the block start/end to calculate relative position for binning
    b_start = start(a_gr)[subjectHits(hits)],
    b_end = end(a_gr)[subjectHits(hits)]
  )
  
  # 4. Create 5 bins per block
  res %>%
    group_by(block) %>%
    # Calculate bin based on relative position within the block
    mutate(bin = cut(pos, 
                     breaks = seq(unique(b_start), unique(b_end), length.out = 6), 
                     labels = FALSE, include.lowest = TRUE)) %>%
    group_by(block, bin) %>%
    summarise(mean_meth = mean(meth, na.rm = TRUE), .groups = "drop")
}

# Apply to all samples
binned_list <- mapply(bin_meth_gr, ldf_filt, ldf_ann_renamed, SIMPLIFY = FALSE)
plot_df_arch6 <- bind_rows(binned_list, .id = "individual")
# --- 4. Mapping Methylation (Smooth plot on one row) ---
# Ensure the blocks stay in the genomic order on the X-axis
ordered_levels <- paste0(arch6_template, 
                         ave(arch6_template, arch6_template, FUN = seq_along))
# This tells ggplot: "Follow THIS order, not alphabetical order"
plot_df_arch6$block <- factor(plot_df_arch6$block, levels = ordered_levels)
ggplot(plot_df_arch6, aes(x = bin, y = mean_meth, group = individual, color = individual)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.5) +
  # nrow = 1 ensures they stay in a single horizontal line
  facet_wrap(~block, nrow = 1, scales = "free_x") + 
  theme_minimal() +
  labs(title = "AZFc Architecture:6 Methylation (Genomic Order)", 
       x = "Bins (5 per block)", 
       y = "Mean Methylation") +
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.1, "lines"), # Keeps blocks physically close
        strip.text = element_text(size = 6, angle = 45), # Tilts names to fit
        legend.position = "none") # Optional: remove if there are too many samples
### ------ || Methylation on AZFc in individuals with unique vs. normal architectures || ------
# Get the list of individuals with unique architectures from the top table
unique_list <- unique(meth_AZFc_arch_all$sample)
# Add a column to plot_df to identify unique vs normal individuals
plot_df <- plot_df %>%
  mutate(arch_status = if_else(individual %in% unique_list, "Unique", "Normal"))
# Filter for AZFc palindromes P1, P2, and P3
azfc_data <- plot_df %>%
  filter(palindrome %in% c("P1", "P2", "P3"))
# Plot the comparison
# Plot the comparison merging P1, P2, and P3
ggplot(azfc_data, aes(x = arch_status, y = avg_meth, fill = arch_status)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.7) +
labs(title = "Total AZFc Methylation (P1-P3): Unique vs. Normal",
       x = "Architecture Type",
       y = "Average Methylation (%)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = "none")
### ------ || Methylation on AZFc vs. other sequences in individuals with unique and normal architectures || ------
#### Analyze specific individual HG03065
hg03065_data <- plot_df %>%
  filter(individual == "HG03065") %>%
  mutate(region_group = if_else(palindrome %in% c("P1", "P2", "P3"), 
                                "AZFc (P1-P3)", 
                                "Other Y Regions"))
# Visualization for HG03065
ggplot(hg03065_data, aes(x = region_group, y = avg_meth, fill = region_group)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, color = "black") +
  labs(title = "Individual HG03065: AZFc vs. Other Y Regions",
       x = "Region Category",
       y = "Average Methylation (%)") +
  theme_classic() +
  guides(fill = "none")
#### 1. Categorize all individuals (Unique IDs vs Normal)
### 1. Categorize individuals as "Normal" vs. "Unique (Combined)"
grouped_arch_data <- plot_df %>%
  # Join architecture metadata
  left_join(meth_AZFc_arch_all %>% select(sample, architecture) %>% distinct(), 
            by = c("individual" = "sample")) %>%
  mutate(
    # Simple logic: If architecture is NA, they are Normal. Otherwise, Unique.
    arch_type = if_else(is.na(architecture), "Normal", "Unique"),
    region_group = if_else(palindrome %in% c("P1", "P2", "P3"), 
                           "AZFc (P1-P3)", 
                           "Other Y Regions")
  ) %>%
  mutate(arch_type = factor(arch_type, levels = c("Normal", "Unique")))

### 2. Plot: Normal vs. Combined Unique (2 Panels)
ggplot(grouped_arch_data, aes(x = region_group, y = avg_meth, fill = region_group)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA) +
  facet_wrap(~ arch_type, nrow = 1) + 
  labs(
    title = "Methylation Comparison: Normal vs. Combined Unique Architectures",
    x = "Region Category",
    y = "Average Methylation (%)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

### 3. Redo Statistical Test
grouped_stats <- grouped_arch_data %>%
  group_by(arch_type) %>%
  filter(n_distinct(region_group) == 2) %>%
  summarise(
    p_value = wilcox.test(avg_meth ~ region_group)$p.value,
    statistic = wilcox.test(avg_meth ~ region_group)$statistic,
    .groups = "drop"
  )

print(grouped_stats)
