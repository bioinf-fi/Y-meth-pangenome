## ------- || METHYLATION SCORE IN RIGHT VS. LEFT ARM + CENTER VS. DISTAL REGION || ---------
### ----- || 0. Rename input list columns globally (NO NEED TO RERUN THIS!) || -----
### Define the precise names the function requires
target_ann_names <- c("chr", "start", "end", "sequence_class", "strand", "V1", "thickStart", "thickEnd", "color")
target_meth_names <- c("chr", "start", "end", "meth", "score_2") 
### Function to rename columns by index for ldf (methylation data)
rename_ldf_globally <- function(df) {
  if (ncol(df) >= 4) {
    names(df)[1:4] <- c("chr", "start", "end", "meth") 
  }
  return(df)
}
### Function to rename columns by index for ldf_ann (annotation data)
rename_ann_globally <- function(df) {
  if (ncol(df) >= 4) {
    names(df)[1:4] <- c("chr", "start", "end", "sequence_class")
  }
  if (ncol(df) == 9) {
    names(df) <- target_ann_names
  }
  return(df)
}
ldf <- lapply(ldf, rename_ldf_globally)
ldf_ann <- lapply(ldf_ann, rename_ann_globally)
print("Columns renamed globally.")
### ----- || 0.2: Check P3, P2, P1 ranges in the sequence_class data - NO NEED TO RERUN THIS || -----
# Modify ldf_ann to embed new P2 container rows
# Define the specific sequences for each arm
sequences <- list(
  P3_Left = c("blue", "teal1", "P3-spacer"),
  P3_Right = c("P3-spacer", "blue", "teal2"),
  P2 = c("red", "P2-spacer", "red"),
  P1_Left = c("gray", "blue", "yellow", "green", "red"),
  P1_Right = c("red", "green", "yellow", "blue", "gray")
)

# Function to detect patterns and inject container rows
inject_arm_containers <- function(df) {
  if (ncol(df) >= 4) names(df)[1:4] <- c("chr", "start", "end", "sequence_class")
  
  new_rows_list <- list()
  
  for (arm_name in names(sequences)) {
    pattern <- sequences[[arm_name]]
    pattern_len <- length(pattern)
    
    # Find all starting indices of the exact sequence pattern
    starts <- which(sapply(1:(nrow(df) - pattern_len + 1), function(i) {
      all(df$sequence_class[i:(i + pattern_len - 1)] == pattern)
    }))
    
    if (length(starts) > 0) {
      arm_type <- gsub("_.*", "", arm_name) # Extract P1, P2, P3
      
      rows_to_add <- map_df(starts, function(i) {
        data.frame(
          chr = df$chr[i],
          start = df$start[i],
          end = df$end[i + pattern_len - 1],
          sequence_class = arm_type, # Label as P1, P2, or P3
          # Fill other columns with NAs to ensure compatibility when binding rows
          stringsAsFactors = FALSE
        ) %>%
          bind_cols(df[i, setdiff(names(df), names(.))]) # Ensure all columns exist
      })
      new_rows_list[[arm_name]] <- rows_to_add
    }
  }
  
  # Combine all new container rows with the original data
  bind_rows(c(list(df), new_rows_list)) %>%
    arrange(chr, start)
}
# Apply this sophisticated injection logic to your list
ldf_ann_2 <- lapply(ldf_ann, inject_arm_containers)
ldf_ann = ldf_ann_2
rm(ldf_ann_2)
### ----- || 1. Global Setup and Metadata Definitions || -----
# Original definitions are correct up to here
sequence_class_order <- c("P8", "P8", "P7", "P7", "P6", "P6", "P5", "P5", "P4", "P4", 
                          "blue", "teal1", "P3-spacer", "teal2", "blue", "green",
                          "red", "P2-spacer", "red", "gray", "blue", "yellow", 
                          "green", "red", "red", "green", "yellow", "blue", "gray")
blocks_order <- c("P8_left", "P8_right", "P7_left", "P7_right", "P6_left", "P6_right", "P5_left", "P5_right", "P4_left", "P4_right", 
                          "blue1", "teal1", "P3-spacer", "teal2", "blue2", "green1",
                          "red1", "P2-spacer", "red2", "gray1", "blue3", "yellow1", 
                          "green2", "red3", "red4", "green3", "yellow2", "blue4", "gray2")
arms_vector <- c("left", "right", "left", "right", "left", "right", "left", "right", "left", 
                 "right", "left", "left", "spacer", "right", "right","spacer", "left", 
                 "spacer", "right", "left", "left", "left", "left", "left", "right", 
                 "right", "right", "right", "right")
palindrome_vector <- c("P8", "P8", "P7", "P7", "P6", "P6", "P5", "P5", "P4", "P4",  
                       "P3", "P3", "P3","P3", "P3", "spacer", "P2", "P2", 
                       "P2", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1")
NUM_BINS <- 10

# The issue is that the reference 'data_arms' table needs to know about
# the *second* occurrence of the entire P2 block later in the genome.
# We need to manually add the second full set of P2 descriptors here:
sequence_class_order_v2 <- c(sequence_class_order, 
                             "red", "P2-spacer", "red") # The second P2 block (rows 6, 7, 8 in your image)
blocks_order_v2 <- c(blocks_order, 
                             "red3", "P2-spacer.2", "red4") # The second P2 block (rows 6, 7, 8 in your image)
arms_vector_v2 <- c(arms_vector, 
                    "left", "spacer", "right") # They are NA in your image output, so we define them as such here
palindrome_vector_v2 <- c(palindrome_vector, 
                         "P2", "P2", "P2") # They belong to P2
data_arms <- data.frame(
  sequence_class = sequence_class_order_v2, 
  arm = arms_vector_v2, 
  palindrome = palindrome_vector_v2,
  blocks_order = blocks_order_v2
) %>%
  # Prepare data_arms for joining by calculating the occurrence ID
  group_by(palindrome, sequence_class) %>%
  mutate(occ_in_p = row_number()) %>%
  ungroup()
### ----- || 2. Create a list containing ONLY the palindrome containers (P1-P8) || -----
# Assumes standardization/P2 injection ran previously, using 'sequence_class' column
P_ldf <- lapply(ldf_ann, function(df) {
  df %>% filter(grepl("^P[1-8]$|^spacer$", sequence_class, ignore.case = TRUE))
})
### ----- || 3. Create a list containing ONLY the color blocks || -----
color_ldf <- lapply(ldf_ann, function(df) { # Use ldf_ann_2 here too
  df %>% filter(sequence_class %in% sequence_class_order)
})
### ----- || 4. Create a combined annotation list by overlapping color blocks with palindromes || -----
combined_ann_ldf <- lapply(names(P_ldf), function(indiv) {
  
  p_df <- P_ldf[[indiv]]
  c_df <- color_ldf[[indiv]]
  
  if(nrow(p_df) == 0 || nrow(c_df) == 0) return(NULL)
  
  # Ensure c_df has the metadata needed for mapping
  p_gr <- GRanges(seqnames = p_df$chr, ranges = IRanges(p_df$start, p_df$end), p_id = p_df$sequence_class)
  c_gr <- GRanges(seqnames = c_df$chr, ranges = IRanges(c_df$start, c_df$end), 
                  c_id = c_df$sequence_class)
  
  hits <- findOverlaps(c_gr, p_gr, type = "any")
  
  mapping <- data.frame(
    chr = as.character(seqnames(c_gr[queryHits(hits)])), 
    start = start(c_gr[queryHits(hits)]), 
    end = end(c_gr[queryHits(hits)]),
    color_block = mcols(c_gr)$c_id[queryHits(hits)], 
    parent_palindrome = mcols(p_gr)$p_id[subjectHits(hits)]
  ) 
  
  df_filtered_mapping <- mapping %>%
    dplyr::filter(
      !color_block %in% c("spacer", "green") | 
        parent_palindrome == color_block
    )
  
  df_joined <- df_filtered_mapping %>%
    dplyr::arrange(chr, start) %>% 
    dplyr::group_by(parent_palindrome, color_block) %>%
    dplyr::mutate(occ_in_p = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    # The join now pulls 'blocks_order' and 'arm' from your data_arms reference table
    dplyr::left_join(data_arms, 
                     by = c("parent_palindrome" = "palindrome", 
                            "color_block" = "sequence_class", 
                            "occ_in_p"))
  
  df_clean <- df_joined %>%
    dplyr::filter(!is.na(arm), !is.na(parent_palindrome)) %>%
    # Ensure blocks_order is kept and positioned clearly
    dplyr::select(chr, start, end, parent_palindrome, color_block, blocks_order, arm, everything())
  
  return(df_clean)
})
names(combined_ann_ldf) <- names(P_ldf)
# AFTER RUNNING THIS CODE, run the following to see the problem rows:
# View(combined_ann_ldf[[1]] %>% filter(parent_palindrome == "P2" | is.na(arm)))
### ----- || 5. Binning Function (Updated to use the unified annotation list) || -----
bin_and_process <- function(meth_df, ann_df_joined, n_bins) {
  names(meth_df)[1:4] <- c("chr", "start", "end", "meth")
  meth_gr <- GRanges(seqnames = meth_df$chr, ranges = IRanges(start = meth_df$start, end = meth_df$end), meth = meth_df$meth)
  results <- list()
  
  for (j in 1:nrow(ann_df_joined)) {
    region <- ann_df_joined[j, ]
    
    # Check for blocks_order as well
    if (is.na(region$parent_palindrome) || is.na(region$arm)) {
      next
    }
    
    bin_width <- (region$end - region$start) / n_bins
    b_starts <- seq(region$start, region$end - bin_width, length.out = n_bins)
    
    # ADDED: Include blocks_order in the bins metadata
    bins_df <- data.frame(
      chr = region$chr, start = b_starts, end = b_starts + bin_width, bin_id = 1:n_bins,
      parent_palindrome = rep(region$parent_palindrome, n_bins),
      arm = rep(region$arm, n_bins),
      blocks_order = rep(region$blocks_order, n_bins) # <--- Carry this over
    )
    
    bins_gr <- makeGRangesFromDataFrame(bins_df, keep.extra.columns = TRUE)
    hits <- findOverlaps(bins_gr, meth_gr)
    
    if (length(hits) > 0) {
      ov_df <- as.data.frame(hits)
      ov_df$meth <- mcols(meth_gr)$meth[ov_df$subjectHits]
      
      # Select blocks_order here so it's preserved after the join
      bin_info <- as.data.frame(bins_gr) %>% 
        dplyr::select(seqnames, start, end, bin_id, parent_palindrome, arm, blocks_order)
      
      names(bin_info)[names(bin_info) == "seqnames"] <- "chr"
      bin_info$queryHits <- as.integer(rownames(bin_info))
      
      ov_df <- dplyr::left_join(ov_df, bin_info, by = "queryHits")
      
      bin_summary <- ov_df %>%
        dplyr::group_by(bin_id, chr, start, end) %>%
        dplyr::summarise(avg_meth = mean(meth, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(
          palindrome = region$parent_palindrome, 
          arm = region$arm,
          blocks_order = region$blocks_order # <--- Ensure it is in the summary
        )
      
      results[[j]] <- bin_summary
    }
  }
  return(dplyr::bind_rows(results))
}
### ----- || 6. Execution Loop (Simplified to use the prepared combined_ann_final) || -----
final_data_list <- list()
for (indiv in names(ldf)) {
  if (indiv %in% names(combined_ann_ldf) && !is.null(combined_ann_ldf[[indiv]])) {
    message("Processing ", indiv)
    m_df <- ldf[[indiv]]
    colnames(m_df)[1:4] <- c("chr", "start", "end", "meth")
    
    # Process the data
    processed_df <- bin_and_process(m_df, combined_ann_ldf[[indiv]], NUM_BINS) %>%
      dplyr::mutate(individual = indiv)
    
    # CRITICAL FIX: Check if processed data is valid before storing
    if (!is.null(processed_df) && nrow(processed_df) > 0) {
      final_data_list[[indiv]] <- processed_df
      n_pal <- dplyr::n_distinct(final_data_list[[indiv]]$palindrome)
      message(paste0("Processed ", indiv, ": Found ", n_pal, " palindromes."))
    } else {
      # If empty, set the list entry to NULL, so bind_rows ignores it
      message(paste0("Processed ", indiv, ": Found 0 palindromes or empty data returned. Setting entry to NULL."))
      final_data_list[[indiv]] <- NULL 
    }
    
  } else {
    message("Skipping ", indiv, " due to missing or empty annotation data.")
  }
}

# CRITICAL FIX: Filter out NULL entries from the list BEFORE binding rows
plot_df <- dplyr::bind_rows(final_data_list[!sapply(final_data_list, is.null)]) %>% 
  filter(arm %in% c("left", "right")) %>%
  mutate(palindrome = factor(palindrome, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")))
# ----- || 7. Visualization || -----
# 1. Rigorous data cleaning to remove rows that cause issues
plot_df_clean <- plot_df %>%
  filter(!is.na(avg_meth), !is.na(bin_id), arm %in% c("left", "right")) %>%
  mutate(
    palindrome = factor(palindrome, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")),
    arm = factor(arm, levels = c("left", "right"))
  ) %>%
  group_by(palindrome, arm, individual) %>%
  filter(n() >= 10) %>%
  ungroup()
# 2. Plot 1: Smoothed Trends (GAM)
plot_smooth <- ggplot(plot_df_clean, aes(x = bin_id, y = avg_meth, color = arm, fill = arm)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2) +
  facet_wrap(~palindrome, ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  labs(
    title = "Methylation Trends Across Palindrome Arms",
    subtitle = "Aggregated trends using GAM smoothing (95% CI)",
    x = "Bin (Proximal 1 -> Distal 10)", 
    y = "Mean Methylation (%)",
    color = "Arm", fill = "Arm"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
print(plot_smooth)
# 3. Plot 2: Unified Boxplot
plot_unified_boxplot <- ggplot(plot_df_clean, aes(x = palindrome, y = avg_meth, fill = arm)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_point(aes(color = arm), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), 
             alpha = 0.2, size = 0.5) +
  labs(
    title = "Methylation Profile: Left vs Right Arms",
    subtitle = "Distribution across all detected Palindromic units",
    x = "Palindrome",
    y = "Average Methylation per Bin (%)",
    fill = "Arm Type", color = "Arm Type"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("left" = "#377eb8", "right" = "#e41a1c")) +
  scale_color_manual(values = c("left" = "#377eb8", "right" = "#e41a1c")) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold")
  )
print(plot_unified_boxplot)
