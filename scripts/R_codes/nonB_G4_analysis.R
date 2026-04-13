## ------ || NONB-DNA ANALYSIS || ------
### ----- || Upload data || -----
# 1. NBMST
all_paths <- list.files(path = "Work/annotation/nbmst", 
                        recursive = TRUE, 
                        full.names = TRUE, 
                        pattern = "\\.tsv$")
# Extract Names (Splitting at the last underscore)
filenames_clean <- gsub("\\.tsv$", "", basename(all_paths))
metadata <- strsplit(filenames_clean, "_(?=[^_]+$)", perl = TRUE)
sample_names <- sapply(metadata, `[`, 1)
dna_types <- sapply(metadata, `[`, 2)
# Read files (Fixing the row.names error)
all_data <- lapply(all_paths, function(x) {
  read.delim(x, header = TRUE, row.names = NULL, sep = "\t")
})
# Group by DNA type and name by Sample
# We'll use a loop to make the naming crystal clear
unique_types <- unique(dna_types)
final_output <- list()
for (type in unique_types) {
  # Find indices belonging to this DNA type
  idx <- which(dna_types == type)
  # Extract those dataframes and name them
  type_list <- all_data[idx]
  names(type_list) <- sample_names[idx]
  # Store in the master list
  final_output[[type]] <- type_list
}
# 1. g4hunter
# 1. Locate the files
# This looks into "g4hunter/" and finds only the "Merged.txt" files in any subfolder
g4h_paths <- list.files(path = "Work/annotation/g4hunter/", 
                        recursive = TRUE, 
                        full.names = TRUE, 
                        pattern = "-Merged\\.txt$")
# 2. Extract Sample Names
# This removes the "-Merged.txt" suffix to leave "SampleX_chrY"
g4h_filenames <- basename(g4h_paths)
g4h_sample_names <- gsub("-Merged\\.txt$", "", g4h_filenames)
# 3. Read files, skipping the first row (>Sample_Chr)
g4h_data_list <- lapply(g4h_paths, function(x) {
  # skip = 1 ignores the metadata line at the top
  read.delim(x, header = TRUE, row.names = NULL, sep = "\t", skip = 1)
})
# Assign the sample names to the list elements
names(g4h_data_list) <- g4h_sample_names
# 4. Final Structure
# Wrapping it in a list called "G4" to match your previous format
final_output_g4hunter <- list("G4" = g4h_data_list)
### ----- || Reorganize data || -----
# NBMST
# 1. & 2. Reformat each file to have one data frame per nonB-DNA type 
# and reorganize columns to "seqnames", "start", "end", "strand"
reorganized_data <- lapply(names(final_output), function(dna_type) {
  # Extract the list of dataframes for this specific DNA type
  dna_list <- final_output[[dna_type]]
  # Combine all individuals into one dataframe and add a 'sample' column
  combined_df <- bind_rows(dna_list, .id = "sample") %>%
    # Rename and select columns as requested
    # Note: 'seqnames' is set to "chrY" for all entries
    mutate(seqnames = "chrY") %>%
    select(
      seqnames, 
      start = Start, 
      end = Stop, 
      strand = Strand, 
      sample, 
      everything() # Keep other columns like Score, Repeat, etc., for analysis
    )
  
  return(combined_df)
})
# Name the top-level list by the DNA type
names(reorganized_data) <- names(final_output)
#g4hunter
# 1. Clean the list: Force EVERY column to be numeric
cleaned_list <- lapply(final_output_g4hunter$G4, function(df) {
  df %>%
    # mutate(across(...)) converts all columns to character then to numeric
    # this safely handles the ">ID" rows by turning them into NAs
    mutate(across(everything(), ~as.numeric(as.character(.x)))) %>%
    # Now we remove those NA rows (using Start as the reference)
    filter(!is.na(Start))
})
# 2. Merge all 141 dataframes
g4hunter_df <- bind_rows(cleaned_list, .id = "Sample_ID")
# 3. Check the result
str(g4hunter_df)
# Split Sample_ID into 'Sample' and 'chr'
# We'll put 'chr' at the start and 'Sample' at the very end as requested
g4hunter_df <- g4hunter_df %>%
  separate(Sample_ID, into = c("Sample", "chr"), sep = "_") %>%
  select(chr, Start, End, Length, Score, Sample) # Reorder columns
head(g4hunter_df)
#overlap nbmst and g4hunter results
# 1. Get unique sample IDs (ensuring they match names)
all_samples <- intersect(unique(g4hunter_df$Sample), unique(gq_annotated$sample))
# 2. Loop through each sample and find overlaps
overlap_list <- lapply(all_samples, function(s) {
  
  # Subset data for this specific sample
  sub_g4h <- g4hunter_df[g4hunter_clean$Sample == s, ]
  sub_nbmst <- gq_annotated[gq_annotated$sample == s, ]
  
  # Create GRanges
  gr_g4h <- GRanges(seqnames = sub_g4h$chr, ranges = IRanges(sub_g4h$Start, sub_g4h$End))
  gr_nbmst <- GRanges(seqnames = sub_nbmst$chr, ranges = IRanges(sub_nbmst$Start, sub_nbmst$End))
  
  # Find Overlaps
  hits <- findOverlaps(gr_g4h, gr_nbmst)
  
  if(length(hits) > 0) {
    # Combine rows side-by-side for this sample
    res <- cbind(sub_g4h[queryHits(hits), ], sub_nbmst[subjectHits(hits), ])
    return(res)
  } else {
    return(NULL)
  }
})
# 3. Bind all sample results back into one dataframe
# We use do.call(rbind, ...) or dplyr::bind_rows
final_overlaps <- do.call(rbind, overlap_list)
# 4. Check results
print(paste("Total confirmed overlaps across all samples:", nrow(final_overlaps)))
head(final_overlaps)

### ----- || NonB-DNA location on the Y || -----
# 1. Ensure the data is merged (using the object from the previous step)
all_obs_merged <- rbindlist(reorganized_data, idcol = "dna_type")
# 2. Plotting with samples side-by-side
ggplot(all_obs_merged, aes(x = dna_type, color = sample)) +
  # Use position_dodge to prevent samples from overlapping each other
  geom_linerange(
    aes(ymin = start, ymax = end), 
    position = position_dodge(width = 0.8), 
    linewidth = 1, 
    alpha = 0.8
  ) + 
  coord_flip() + 
  theme_minimal() +
  # Use a large palette if you have many samples (e.g., rainbow or viridis)
  scale_color_discrete() + 
  theme(
    legend.position = "none", # Removes the legend as requested
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "nonB-DNA Location Conservation across Individuals",
    subtitle = "Color-coded by Sample | Grouped by Motif Type",
    y = "Genomic Position on chrY (bp)"
  )
### ----- || Convert lists to master data.tables with a 'sample' column || -----
# Convert lists to master data.tables
classes_master <- rbindlist(ldf_ann, idcol = "sample")
meth_master    <- rbindlist(ldf, idcol = "sample")
nonb_master    <- rbindlist(reorganized_data, idcol = "dna_type") 

# Clean sample names: Strip "_chrY" to ensure IDs match across lists
nonb_master[, sample    := sub("_chrY$", "", sample)]
classes_master[, sample := sub("_chrY$", "", sample)]
meth_master[, sample    := sub("_chrY$", "", sample)]

# Helper function to remove NAs and ensure numeric coordinates
clean_genomic_dt <- function(dt) {
  dt <- dt[!is.na(start) & !is.na(end)]
  dt[, `:=`(start = as.numeric(start), end = as.numeric(end))]
  return(dt)
}
nonb_master    <- clean_genomic_dt(nonb_master)
classes_master <- clean_genomic_dt(classes_master)
meth_master    <- clean_genomic_dt(meth_master)
# Set keys for fast genomic overlaps (Essential for foverlaps)
setkey(classes_master, sample, start, end)
setkey(nonb_master, sample, start, end)
setkey(meth_master, sample, start, end)
### ----- || Analysis: % nonB-DNA on sequence classes || -----
# Join nonB motifs with their sequence class
# Find which sequence class (PAR, XTR, etc.) each nonB-DNA motif belongs to
nonb_with_class <- foverlaps(nonb_master, classes_master, 
                             by.x = c("sample", "start", "end"), type = "any")
nonb_with_class[is.na(sequence_class), sequence_class := "Other"]
# Calculate percentages per DNA type
class_summary <- nonb_with_class[, .(count = .N), by = .(sequence_class, dna_type)]
class_summary[, percentage := (count / sum(count)) * 100, by = dna_type]
ggplot(class_summary, aes(x = sequence_class, y = percentage, fill = dna_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "nonB-DNA Distribution across Sequence Classes", 
    x = "Sequence Class", 
    y = "% of Motif Type", 
    fill = "nonB Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
### ----- || Number of G4s on different sequence classes || -----
seq_classes <- c("PAR1", "XDR", "XTR", "IR", "SAT", "CEN","DYZ17","P", "DYZ19", "HET", "PAR2")
seq_classes_v2 <- c("PAR1", "XDR", "XTR", "IR", "SAT", "CEN","DYZ17","P1","P2","P3","P4","P5","P6","P7","P8", "DYZ19", "HET", "PAR2")
# 1. Clean the G4 dataframe sample column
gq_annotated <- reorganized_data[["GQ"]] %>%
  mutate(sample = sub("_chrY$", "", sample))
# --- 1. AGGREGATE DATA GLOBALLY ---
# We group by sequence class across ALL samples
g4_global_summary <- gq_annotated %>%
  mutate(category = case_when(
    grepl("PAR1", sequence_class)  ~ "PAR1",
    grepl("XDR", sequence_class)   ~ "XDR",
    grepl("XTR", sequence_class)   ~ "XTR",
    grepl("IR", sequence_class)    ~ "IR",
    grepl("SAT", sequence_class)   ~ "SAT",
    grepl("CEN", sequence_class)   ~ "CEN",
    grepl("DYZ17", sequence_class) ~ "DYZ17",
    grepl("AMPL", sequence_class)  ~ "AMPL", # Adjusted 'P' to 'AMPL' based on your previous code
    grepl("DYZ19", sequence_class) ~ "DYZ19",
    grepl("HET", sequence_class)   ~ "HET",
    grepl("PAR2", sequence_class)  ~ "PAR2",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("PAR1", "XDR", "XTR", "IR", "SAT", "CEN",
                                                "DYZ17", "AMPL", "DYZ19", "HET", "PAR2")))

# --- 2. DEFINE THE MATCHING PALETTE ---
class_colors <- c(
  PAR1  = "#f4a261",
  XDR   = "#e76f51",
  XTR   = "#ffbe0b",
  IR    = "#90be6d",
  SAT   = "#adb5bd",
  CEN   = "#6c757d",
  DYZ17 = "#457b9d",
  AMPL  = "#c77dff",
  DYZ19 = "#4cc9f0",
  HET   = "#9d4edd",
  PAR2  = "#f28482"
)

# --- 3. CREATE CIRCULAR BARPLOT ---
ggplot(g4_global_summary, aes(x = category, y = total_count, fill = category)) +
  # Use geom_col for the bars
  geom_col(width = 1, color = "white", size = 0.3) + 
  # This transforms the plot into a circle
  coord_polar(start = 0) + 
  scale_fill_manual(values = class_colors) +
  labs(
    title = "Global G4 Distribution by Sequence Class",
    subtitle = "Total counts across all individuals",
    x = NULL,
    y = "Total G4 Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "none" # Labels are on the axis, so legend is optional
  )
### ----- || G4s distribution on sequence classes across samples || -----
# Ensure this matches your tree labels and metadata
seq_classes <- c("PAR1", "XDR", "XTR", "IR", "SAT", "CEN","DYZ17","AMPL", "DYZ19", "HET", "PAR2")
gq_dt <- as.data.table(g4hunter_df)
colnames(gq_dt) <- c("chr","start",'end','length','score','sample')
gq_dt[, sample := sub("_chrY$", "", sample)]
ann_dt <- rbindlist(ldf_ann, idcol = "sample")
setnames(ann_dt, "chr", "seqnames")
# Clean NAs for foverlaps
gq_dt_clean <- gq_dt[!is.na(start) & !is.na(end)]
ann_dt_clean <- ann_dt[!is.na(start) & !is.na(end)]
# Set keys for the overlap (sample, start, end)
setkey(gq_dt_clean, sample, start, end)
setkey(ann_dt_clean, sample, start, end)
# --- 2. DATA.TABLE OVERLAP (G4s vs Annotations) ---
# This finds which G4 falls into which sequence class
gq_annotated <- foverlaps(gq_dt_clean, ann_dt_clean, 
                          by.x = c("sample", "start", "end"), 
                          by.y = c("sample", "start", "end"), 
                          type = "within", 
                          nomatch = NULL)
# --- 3. CATEGORIZE AND SUMMARIZE ---
pattern <- paste(seq_classes, collapse = "|")
gq_grouped_final <- gq_annotated[grepl(pattern, sequence_class)] %>%
  .[, .(count = .N), by = .(sample, sequence_class)] %>%
  # Clean up sequence_class names to match your levels
  .[, category := tstrsplit(sequence_class, "_", keep = 1)] %>% 
  .[category %in% seq_classes] %>%
  .[, .(count = sum(count)), by = .(sample, category)] %>%
  # Calculate percentages per sample
  .[, percentage := (count / sum(count)) * 100, by = sample]
# --- 4. FORMAT FOR GGTREE ---
gq_grouped_final <- as_tibble(gq_grouped_final) %>%
  mutate(
    sample = factor(sample, levels = my_tree$tip.label),
    category = factor(category, levels = seq_classes)
  )
# --- 5. BUILD CIRCULAR PLOT ---
p <- ggtree(my_tree, layout = "circular", branch.length = "none", size = 0.2) +
  geom_fruit(
    data = gq_grouped_final,
    geom = geom_col,
    mapping = aes(y = sample, x = percentage, fill = category),
    orientation = "y",
    offset = 0.05,
    pwidth = 2, # Increase this to make the bars wider/more visible
    color = NA  # Remove outlines for a smoother look like your image
  ) +
  scale_fill_manual(values = class_colors, name = "Sequence class") +
  theme(legend.position = "right")

print(p)
### POLAR VIEW
# 1. Aggregate counts globally across all samples
gq_polar_data <- gq_grouped_final %>%
  group_by(category) %>%
  summarise(total_count = sum(count)) %>%
  mutate(category = factor(category, levels = seq_classes))

# 2. Build the Rose Chart
p_rose <- ggplot(gq_polar_data, aes(x = category, y = total_count, fill = category)) +
  # Draw the bars
  geom_col(width = 1, color = "white", size = 0.2) + 
  # Transform to polar coordinates to make it circular
  coord_polar(start = 0) +
  # Apply your specific class colors
  scale_fill_manual(values = class_colors) +
  # Clean up the theme to match the reference image
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    legend.position = "none" # Image has labels on the plot, so legend is often redundant
  ) +
  # Add the labels around the circle
  labs(title = "G4 Distribution by Sequence Class")

print(p_rose)
### ----- || Barplot of most represented G4 sequences || -----
# --- 1. CALCULATE TOP 10 FREQUENCIES ---
# Identify the top 10 sequences globally
top_10_global <- gq_annotated[, .(count = .N), by = Sequence][order(-count)][1:10]
# Calculate percentage relative to these 10 sequences
top_10_global[, percentage := (count / sum(count)) * 100]
# Set factor levels for sorting in the legend/plot
top_10_global[, Sequence := factor(Sequence, levels = Sequence)]
# --- 2. DEFINE VISUALLY MATCHING PALETTE ---
# These are the specific hex codes from your provided image/code
top10_colors <- c(
  "#f4a261", # Orange (PAR1)
           "#e76f51", # Coral (XDR)
           "#ffbe0b", # Yellow (XTR)
           "#90be6d", # Green (IR)
           "#457b9d", # Steel Blue (DYZ17)
           "#c77dff", # Lavender (AMPL)
           "#4cc9f0", # Sky Blue (DYZ19)
           "#9d4edd", # Deep Purple (HET)
           "#f28482", # Salmon/Pink (PAR2)
           "#6c757d"  # Dark Gray (CEN)
)
# --- 3. CREATE THE STACKED BARPLOT ---
ggplot(top_10_global, aes(x = "Overall Distribution", y = percentage, fill = Sequence)) +
  geom_col(width = 0.6, color = "white", size = 0.2) +
  scale_fill_manual(values = top10_colors) +
  labs(
    title = "Top 10 Most Frequent G4 Sequences",
    subtitle = "Global distribution across all samples",
    y = "Percentage of Top 10 (%)",
    x = NULL,
    fill = "Sequence Pattern"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(family = "mono", size = 9), # Mono font for DNA sequences
    legend.title = element_text(face = "bold")
  )
### ----- || G4s density along different Y chromosomes || -----
# 1. Prepare G4 data as a list of Data Tables
colnames(g4hunter_df) = c("chr","start",'end','length','score','sample')
gq_dt <- as.data.table(g4hunter_df)
g4_list_per_sample <- split(gq_dt, by = "sample")
# 2. Process into a 5000-bin DENSITY matrix
n_bins <- 5000
common_samples <- intersect(names(g4_list_per_sample), names(ldf_ann))
common_samples <- intersect(common_samples, phylog_annotation$sample)
g4_density_list <- map(common_samples, function(s) {
  sample_g4 <- g4_list_per_sample[[s]]
  ann_df <- ldf_ann[[s]]
  
  if (is.null(sample_g4) || nrow(sample_g4) == 0 || is.null(ann_df) || nrow(ann_df) == 0) {
    return(rep(0, n_bins))
  }
  
  # Normalize to sample-specific Y length
  max_len <- max(ann_df$end, na.rm = TRUE)
  
  # Calculate relative positions for G4 starts
  rel_pos <- sample_g4$start / max_len
  
  # Assign each G4 to a bin (1 to 5000)
  bin_idx <- cut(rel_pos, breaks = seq(0, 1, length.out = n_bins + 1), 
                 labels = FALSE, include.lowest = TRUE)
  
  # Count G4s per bin (Density)
  # We use tabulate for speed over 5000 bins
  bin_counts <- tabulate(bin_idx, nbins = n_bins)
  
  return(bin_counts)
})
# 3. Create Matrix and Sync with Tree
g4_matrix <- do.call(rbind, g4_density_list)
rownames(g4_matrix) <- common_samples
common_ids <- intersect(rownames(g4_matrix), my_tree$tip.label)
g4_matrix <- g4_matrix[common_ids, ]
my_tree_pruned <- keep.tip(my_tree, common_ids)
dend <- as.dendrogram(my_tree_pruned)
# 4. Define Density Gradient (White to Violet)
# Using 95th percentile to handle outliers so the heatmap isn't too faint
max_val <- quantile(g4_matrix, 0.95, na.rm = TRUE)
if(max_val == 0) max_val <- max(g4_matrix) 
col_density = colorRamp2(c(0, max_val), c("#FFFFFF", "#6A0DAD"))
# 5. Build Heatmap
h_g4 = Heatmap(g4_matrix, 
               name = "G4 Count", 
               cluster_rows = dend,           
               row_dend_side = "left",        
               row_dend_width = unit(4, "cm"), 
               cluster_columns = FALSE,       
               col = col_density,             
               show_row_names = FALSE,        
               show_column_names = FALSE,
               use_raster = TRUE,             
               raster_quality = 50,
               column_title = "Y-Chromosome G4 Density (Normalized Bins)"
)
draw(h_g4)
### ----- || G4s stability along different Y chromosomes || -----
# 1. Synchronize IDs across all datasets
# Ensure the tree, G4 data, and annotations all share the same samples
common_samples <- intersect(names(g4_list_per_sample), names(ldf_ann))
# 2. Process into a 5000-bin STABILITY matrix
n_bins <- 5000
stability_list <- lapply(common_samples, function(s) {
  sample_g4 <- g4_list_per_sample[[s]]
  ann_df    <- ldf_ann[[s]]
  
  if (is.null(sample_g4) || nrow(sample_g4) == 0 || is.null(ann_df)) {
    return(rep(NA, n_bins))
  }
  
  # DYNAMIC CALCULATION: Get max length from sample-specific annotation
  max_len <- max(ann_df$end, na.rm = TRUE)
  
  # Normalize positions (0 to 1) to account for varying Y lengths
  rel_pos <- as.numeric(sample_g4$start) / max_len
  
  # Binning
  bin_idx <- cut(rel_pos, 
                 breaks = seq(0, 1, length.out = n_bins + 1), 
                 labels = FALSE, include.lowest = TRUE)
  
  # Calculate MEAN stability score per bin
  bin_data <- data.table(bin_id = bin_idx, score = as.numeric(sample_g4$score))
  bin_means <- bin_data[, .(m = mean(score, na.rm = TRUE)), by = bin_id]
  
  # Fill vector with NA so empty bins don't look like '0' stability
  full_vec <- rep(NA, n_bins)
  full_vec[bin_means$bin_id] <- bin_means$m
  return(full_vec)
})

# 3. Final Matrix and Tree Integration
stability_matrix <- do.call(rbind, stability_list)
rownames(stability_matrix) <- common_samples

# Prune tree to match the final sample set
pruned_tree <- keep.tip(my_tree, common_samples)
dend_pruned <- as.dendrogram(pruned_tree)

# 4. Prepare Geographic Metadata and Color Map
row_regions_filtered <- ordered_metadata$Main_geographic_prevalence[match(common_samples, ordered_metadata$ID)]
row_regions_filtered[is.na(row_regions_filtered)] <- "Other"

# Generate color map dynamically to prevent "object not found" error
unique_regions <- unique(row_regions_filtered)
region_colors <- colorRampPalette(brewer.pal(min(length(unique_regions), 12), "Set3"))(length(unique_regions))
region_color_map <- setNames(region_colors, unique_regions)

# 5. Define Heatmap Components
left_ann = rowAnnotation(
  Region = row_regions_filtered,
  col = list(Region = region_color_map),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.5, "cm")
)

# Colors: Blue = Less Stable, White = Neutral, Red = Highly Stable
col_stability = colorRamp2(c(-1, 0, 1), c("#0077b6", "#FFFFFF", "#e63946"))

# 6. Final Heatmap Plot
h_stability = Heatmap(stability_matrix, 
                      name = "Avg Stability", 
                      
                      # Tree/Row Settings
                      cluster_rows = dend_pruned,      
                      row_dend_side = "left",        
                      row_dend_width = unit(4, "cm"), 
                      left_annotation = left_ann,
                      
                      # Column/Genomic Settings
                      cluster_columns = FALSE,       
                      col = col_stability,           
                      na_col = "#f7f7f7", # Light grey for missing data
                      
                      # Performance/Aesthetics
                      use_raster = TRUE,             
                      raster_quality = 50,
                      border = TRUE,
                      column_title = "Y-Chromosome G4 Stability Profile"
)

# Draw the plot
draw(h_stability)
### ----- || number of G4s on PAR1 and PAR2 || -----
# --- 1. EXTRACT AND BIN PAR REGIONS ---
par_distribution <- lapply(unique(gq_annotated$sample), function(indiv) {
  
  # 1. Get and clean sample data
  sample_data <- gq_annotated %>%
    filter(sample == indiv, sequence_class %in% c("PAR1", "PAR2")) %>%
    filter(!is.na(start), !is.na(end))
  
  if (nrow(sample_data) == 0) return(NULL)
  
  # 2. Summarize boundaries
  par_ann <- sample_data %>%
    group_by(sequence_class) %>%
    summarise(start = min(start), end = max(end), .groups = "drop")
  
  # 3. CORRECTION 1: Check if PAR2 is contained within PAR1
  if (all(c("PAR1", "PAR2") %in% par_ann$sequence_class)) {
    p1 <- par_ann %>% filter(sequence_class == "PAR1")
    p2 <- par_ann %>% filter(sequence_class == "PAR2")
    
    # If PAR2 start and end are inside PAR1, discard PAR2
    if (p2$start >= p1$start && p2$end <= p1$end) {
      par_ann <- par_ann %>% filter(sequence_class == "PAR1")
    }
  }
  
  # 4. Get G4 data for this individual
  g4_sub <- g4_df_clean %>% filter(sample == indiv)
  if (nrow(g4_sub) == 0) return(NULL)
  
  g4_gr <- GRanges(seqnames = "chrY", ranges = IRanges(g4_sub$start, g4_sub$end))
  
  bin_results <- list()
  NUM_BINS <- 20
  
  for(i in 1:nrow(par_ann)) {
    reg <- par_ann[i,]
    bin_width <- (reg$end - reg$start) / NUM_BINS
    b_starts <- seq(reg$start, by = bin_width, length.out = NUM_BINS)
    
    bins_df <- data.frame(
      bin_id = 1:NUM_BINS,
      start = b_starts,
      end = b_starts + bin_width,
      region = reg$sequence_class
    )
    
    bins_gr <- GRanges(seqnames = "chrY", 
                       ranges = IRanges(round(bins_df$start), round(bins_df$end)))
    
    # Count G4s
    bins_df$g4_count <- countOverlaps(bins_gr, g4_gr)
    
    # CORRECTION 2: Orientation (Bin 1 = Telomere)
    # PAR1: Start is telomere. Bin 1 is already at the start.
    # PAR2: End is telomere. We reverse bin_id so Bin 1 is the highest coordinate (the telomere).
    if(reg$sequence_class == "PAR2") {
      bins_df$bin_id <- rev(bins_df$bin_id)
    }
    
    bin_results[[i]] <- bins_df
  }
  
  return(bind_rows(bin_results) %>% mutate(sample = indiv))
}) %>% bind_rows()

# --- 2. PLOT PAR1 vs PAR2 SMOOTHED ---
ggplot(par_distribution, aes(x = bin_id, y = g4_count, color = region, fill = region)) +
  
  # 1. BOTTOM LAYER: Bars (Mean G4 count across all individuals)
  stat_summary(
    fun = mean, 
    geom = "col", 
    alpha = 0.4,       # Light transparency so it's a background "bar"
    width = 0.8, 
    color = NA         # Remove bar outlines for a cleaner look
  ) +
  
  # 2. MIDDLE LAYER: Individual points (highly transparent)
  geom_jitter(
    height = 0, 
    width = 0.2, 
    alpha = 0.1, 
    size = 0.5,
    show.legend = FALSE
  ) +
  
  # 3. TOP LAYER: Smoothed trend line (drawn over the bars)
  geom_smooth(
    method = "loess", 
    span = 0.4, 
    alpha = 0.2,       # Shaded confidence interval
    linewidth = 1.2    # Thicker line to stand out
  ) +
  
  # Faceting and Aesthetics
  facet_wrap(~region) +
  scale_color_manual(values = c("PAR1" = "#f4a261", "PAR2" = "#f28482")) +
  scale_fill_manual(values = c("PAR1" = "#f4a261", "PAR2" = "#f28482")) +
  
  labs(
    title = "G4 Distribution: PAR1 vs PAR2",
    subtitle = "Bars = Mean Count | Line = Smoothed Trend (Loess)",
    x = "Relative Position (Telomere → MSY)",
    y = "G4 Count per Bin"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank() # Clean up vertical lines
  )
### ----- || G4s stability on PAR1 and PAR2 || -----
# Convert to data.table for speed
colnames(g4hunter_df) = c("chr","start",'end','length','score','sample')
gq_dt <- as.data.table(g4hunter_df)
# Create the list per individual (Sample column)
g4_list_per_sample <- split(gq_dt, by = "sample")
y_lengths <- sapply(ldf_ann, function(x) max(x$end, na.rm = TRUE))
par_stability_final <- lapply(names(g4_list_per_sample), function(indiv) {
  # 1. Access the specific list element
  # [[indiv]] extracts the actual data frame or GRanges
  ann_content <- ldf_ann[[indiv]]
  if (is.null(ann_content)) return(NULL)
  
  # 2. Convert to a standard data.table
  sample_ann_dt <- as.data.table(ann_content)
  
  # 3. Identify the class column (Checking common names used in your previous code)
  # Change "sequence_class" if the column in your list has a different name
  potential_cols <- c("sequence_class", "type", "category", "name")
  class_col <- intersect(potential_cols, names(sample_ann_dt))[1]
  
  if (is.na(class_col)) return(NULL)
  
  # 4. Filter for PAR regions
  # Using standard subsetting to avoid variable name confusion
  par_regions <- sample_ann_dt[sample_ann_dt[[class_col]] %in% c("PAR1", "PAR2"), ]
  
  if (nrow(par_regions) == 0) return(NULL)
  
  sample_g4s <- g4_list_per_sample[[indiv]]
  NUM_PAR_BINS <- 20
  bin_results <- list()
  
  for(i in 1:nrow(par_regions)) {
    # Force coordinates to numeric to avoid "atomic type" errors
    c_start <- as.numeric(par_regions$start[i])
    c_end   <- as.numeric(par_regions$end[i])
    c_label <- as.character(par_regions[[class_col]][i])
    
    # Subset G4s strictly within the PAR boundaries
    # Using $start ensures we look at the column, not the function start()
    reg_g4s <- sample_g4s[sample_g4s$start >= c_start & sample_g4s$end <= c_end, ]
    
    if (nrow(reg_g4s) > 0) {
      bin_edges <- seq(c_start, c_end, length.out = NUM_PAR_BINS + 1)
      # Assign bins manually to avoid variable scoping issues
      current_bins <- cut(as.numeric(reg_g4s$start), breaks = bin_edges, labels = FALSE, include.lowest = TRUE)
      reg_g4s[, bin_id := current_bins]
      
      bin_avg <- reg_g4s[, .(avg_stability = mean(score, na.rm = TRUE)), by = bin_id]
    } else {
      bin_avg <- data.table(bin_id = integer(), avg_stability = numeric())
    }
    
    # Ensure a full set of 20 bins for the visualization
    skeleton <- data.table(bin_id = 1:NUM_PAR_BINS, region = c_label)
    bin_avg <- merge(skeleton, bin_avg, by = "bin_id", all.x = TRUE)
    bin_avg[is.na(avg_stability), avg_stability := 0]
    
    # Flip orientation for PAR2 (Telomere to MSY)
    if (c_label == "PAR2") bin_avg[, bin_id := rev(bin_id)]
    
    bin_results[[i]] <- bin_avg
  }
  
  return(rbindlist(bin_results)[, sample := indiv])
}) %>% rbindlist()
# --- 2. PLOT PAR1 vs PAR2 G4 STABILITY ---
ggplot(par_stability_final, aes(x = bin_id, y = avg_stability, color = region, fill = region)) +
  
  # 1. BOTTOM LAYER: Bars (Mean Stability across all individuals per bin)
  stat_summary(
    fun = mean, 
    geom = "col", 
    alpha = 0.4,       
    width = 0.8, 
    color = NA         
  ) +
  
  # 2. MIDDLE LAYER: Individual points (Stability per bin per sample)
  # Using geom_jitter helps visualize the spread of individual sample stabilities
  geom_jitter(
    height = 0, 
    width = 0.2, 
    alpha = 0.1, 
    size = 0.5,
    show.legend = FALSE
  ) +
  
  # 3. TOP LAYER: Smoothed trend line
  geom_smooth(
    method = "loess", 
    span = 0.4, 
    alpha = 0.2,       
    linewidth = 1.2    
  ) +
  
  # Faceting and Aesthetics
  facet_wrap(~region) +
  scale_color_manual(values = c("PAR1" = "#f4a261", "PAR2" = "#f28482")) +
  scale_fill_manual(values = c("PAR1" = "#f4a261", "PAR2" = "#f28482")) +
  
  labs(
    title = "G4 Stability Profile: PAR1 vs PAR2",
    subtitle = "Bars = Mean Stability | Points = Individual Samples | Line = Loess Trend",
    x = "Relative Position (Telomere → MSY)",
    y = "Average G4 Stability (G4Hunter Score)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

### ----- || number of G4s across palindromes || -----
# --- 0. DATA.TABLE OVERLAP (GQs vs Annotations) ---
gq_dt <- as.data.table(g4hunter_df)
colnames(gq_dt) <- c("chr","start",'end','length','score','sample')
gq_dt[, sample := sub("_chrY$", "", sample)]
ann_dt <- rbindlist(ldf_ann, idcol = "sample")
setnames(ann_dt, "chr", "seqnames")
# Clean NAs for foverlaps
gq_dt_clean <- gq_dt[!is.na(start) & !is.na(end)]
ann_dt_clean <- ann_dt[!is.na(start) & !is.na(end)]
setkey(gq_dt_clean, sample, start, end)
setkey(ann_dt_clean, sample, start, end)
gq_annotated <- foverlaps(gq_dt_clean, ann_dt_clean, 
                          by.x = c("sample", "start", "end"), 
                          by.y = c("sample", "start", "end"), 
                          type = "within", 
                          nomatch = NULL)
# --- CATEGORIZE BY SEQ_CLASSES ---
pattern <- paste(seq_classes_v2, collapse = "|")
gq_stacked_data <- gq_annotated[grepl(pattern, sequence_class), .(count = .N), by = .(sample, sequence_class)]
gq_stacked_data[, percentage := (count / sum(count)) * 100, by = sample]
gq_stacked_final <- as_tibble(gq_stacked_data) %>% rename(category = sequence_class)
# --- 1. FILTER FOR PALINDROMES P1-P8 ONLY ---
gq_annotated_filtered = gq_annotated[which(gq_annotated$score>1),]
g4_palindromes_summary <- gq_annotated_filtered %>%
  mutate(category = case_when(
    grepl("P1", sequence_class) ~ "P1",
    grepl("P2", sequence_class) ~ "P2",
    grepl("P3", sequence_class) ~ "P3",
    grepl("P4", sequence_class) ~ "P4",
    grepl("P5", sequence_class) ~ "P5",
    grepl("P6", sequence_class) ~ "P6",
    grepl("P7", sequence_class) ~ "P7",
    grepl("P8", sequence_class) ~ "P8",
    TRUE ~ NA_character_
  )) %>%
  # Keep only the rows that matched P1-P8
  filter(!is.na(category)) %>%
  # Strictly define the order for the circle
  mutate(category = factor(category, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"))) %>%
  group_by(category) %>%
  summarise(total_count = n(), .groups = "drop")
# --- 2. DEDICATED PALETTE (Purple/Blue Gradient) ---
p_colors <- c(
  P1 = "#cdb4db", P2 = "#ffc8dd", P3 = "#ffafcc", P4 = "#bde0fe", 
  P5 = "#a2d2ff", P6 = "#af4d98", P7 = "#80669d", P8 = "#822faf"
)
# --- 3. CIRCULAR BARPLOT ---
# --- 1. DEFINE GRID RANGES ---
# Adjust the 'to' value based on your maximum G4 count
grid_values <- seq(0, 50000, by = 10000)
# --- 2. UPDATED CIRCULAR BARPLOT ---
ggplot(g4_palindromes_summary, aes(x = category, y = total_count, fill = category)) +
  # Add background concentric circles (Grid lines)
  geom_hline(yintercept = grid_values, color = "gray90", size = 0.3) +
  
  # Draw the actual bars
  geom_col(width = 0.85, color = "white", linewidth = 0.4) + 
  
  # Labels for the rings (optional, shows the 1000, 2000, etc. scale)
  geom_text(data = data.frame(y = grid_values), 
            aes(x = 0.5, y = y, label = y), 
            inherit.aes = FALSE, color = "gray60", size = 2.5, vjust = -0.5) +
  
  # Transform to circular
  coord_polar(start = 0) + 
  
  # Ensure the Y-axis (radius) has enough room for labels
  expand_limits(y = max(g4_palindromes_summary$total_count) * 1.2) +
  
  scale_fill_manual(values = p_colors) +
  
  # Text labels for the specific counts on top of bars
  geom_text(aes(label = total_count), 
            vjust = -0.8, 
            size = 3.5, 
            fontface = "bold") +
  
  labs(
    title = "G4 Counts: Palindromes P1-P8",
    subtitle = "Concentric rings at 1,000 G4 intervals",
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_blank(),
    panel.grid = element_blank(), # We use our manual geom_hline instead
    legend.position = "none"
  )

### ----- || G4 Sequence/length variability across palindromes || -----
## sequence length
# --- 1. PREPARE THE DATA ---
# Filter for P1-P8 and ensure the correct order
p_length_data <- gq_annotated %>%
  mutate(category = case_when(
    grepl("P1", sequence_class) ~ "P1",
    grepl("P2", sequence_class) ~ "P2",
    grepl("P3", sequence_class) ~ "P3",
    grepl("P4", sequence_class) ~ "P4",
    grepl("P5", sequence_class) ~ "P5",
    grepl("P6", sequence_class) ~ "P6",
    grepl("P7", sequence_class) ~ "P7",
    grepl("P8", sequence_class) ~ "P8",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(category)) %>%
  # Use factor levels to keep P1-P8 in order on the X-axis
  mutate(category = factor(category, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")))
# --- 2. CREATE THE VIOLIN PLOT ---
ggplot(p_length_data, aes(x = category, y = Length, fill = category)) +
  # Main distribution shape
  geom_violin(trim = FALSE, alpha = 0.7, color = "white") +
  # Add individual points (jittered) to see each sample
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  # Add a narrow boxplot inside to highlight the median and quartiles
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.4) +
  # Use the same color palette as before for consistency
  scale_fill_manual(values = p_colors) + 
  labs(
    title = "Sequence Length Variability: Palindromes P1-P8",
    subtitle = "Distribution of lengths across all individuals",
    x = "Palindrome Class",
    y = "Sequence Length (bp)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
## sequence
# --- 1. CALCULATE UNIQUE SEQUENCES PER PALINDROME ---
sequence_variability <- gq_annotated %>%
  mutate(category = case_when(
    grepl("P1", sequence_class) ~ "P1",
    grepl("P2", sequence_class) ~ "P2",
    grepl("P3", sequence_class) ~ "P3",
    grepl("P4", sequence_class) ~ "P4",
    grepl("P5", sequence_class) ~ "P5",
    grepl("P6", sequence_class) ~ "P6",
    grepl("P7", sequence_class) ~ "P7",
    grepl("P8", sequence_class) ~ "P8",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  # Assuming the column containing the actual DNA sequence is named 'sequence'
  # We count how many distinct/unique strings exist in that column
  summarise(unique_seq_count = n_distinct(Sequence), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")))
# --- 2. PLOT SEQUENCE VARIABILITY ---
ggplot(sequence_variability, aes(x = category, y = unique_seq_count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = unique_seq_count), vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = p_colors) +
  labs(
    title = "Sequence Polymorphism in Y-Palindromes",
    subtitle = "Number of unique DNA sequences identified across all samples",
    x = "Palindrome Class",
    y = "Count of Unique Sequences"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )
### ----- || G4s density along the AZFc region || -----
# 1. Filter G4s
gq_filtered = g4hunter_df[g4hunter_df$score > 1, ]

g4_distribution_data <- lapply(names(combined_ann_ldf), function(indiv) {
  
  g4_df <- gq_filtered[gq_filtered$sample == indiv, ]
  ann_df <- combined_ann_ldf[[indiv]] %>% 
    filter(!is.na(start), !is.na(end))
  
  if(nrow(g4_df) == 0 || nrow(ann_df) == 0) return(NULL)
  
  # CRITICAL FIX: Ensure naming consistency (remove "chr" from BOTH)
  g4_seq <- gsub("chr", "", as.character(g4_df$chr))
  ann_seq <- gsub("chr", "", as.character(ann_df$chr))
  
  # Create G4 Ranges
  g4_gr <- GRanges(seqnames = g4_seq, 
                   ranges = IRanges(g4_df$start, g4_df$end))
  
  NUM_BINS <- 5 
  
  # 2. Optimized Binning Logic
  all_bins <- ann_df %>%
    rowwise() %>%
    do({
      region <- .
      span <- max(0, region$end - region$start)
      
      # Generate exact break points
      b_breaks <- seq(region$start, region$end, length.out = NUM_BINS + 1)
      
      data.frame(
        bin_id = 1:NUM_BINS,
        start = head(b_breaks, -1),
        end = tail(b_breaks, -1),
        blocks_order = region$blocks_order,
        palindrome = region$parent_palindrome,
        arm = region$arm,
        chr = gsub("chr", "", as.character(region$chr)),
        stringsAsFactors = FALSE
      )
    }) %>% ungroup()
  
  # 3. Create GRanges for all bins at once
  bins_gr <- GRanges(seqnames = all_bins$chr, 
                     ranges = IRanges(floor(all_bins$start), ceiling(all_bins$end)))
  
  # Synchronize seqlevels to avoid "no common levels" warning
  seqlevels(bins_gr) <- unique(c(seqlevels(bins_gr), seqlevels(g4_gr)))
  
  # 4. Perform Overlap
  all_bins$g4_count <- countOverlaps(bins_gr, g4_gr)
  
  return(all_bins %>% mutate(sample = indiv))
})

# Combine everything
final_g4_dist <- bind_rows(g4_distribution_data)
# --- 2. MIRRORED SYMMETRY PLOT ---
mirrored_dist <- final_g4_dist %>%
  # Flip bin IDs for the Left arm so it mirrors the Right
  mutate(adjusted_bin = ifelse(arm == "left", -bin_id, bin_id))
ggplot(mirrored_dist, aes(x = adjusted_bin, y = g4_count, color = palindrome)) +
  geom_smooth(method = "loess", span = 0.5, se = TRUE) +
  # Add a vertical line at 0 to represent the palindrome center
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~palindrome, ncol = 2) +
  scale_color_manual(values = p_colors) +
  labs(
    title = "Palindrome Symmetry: G4 Distribution",
    subtitle = "Left Arm (negative axis) vs Right Arm (positive axis)",
    x = "← Left Arm | Right Arm →",
    y = "Smoothed G4 Count"
  ) +
  theme_minimal()
### ----- || Analysis: methylation on nonB vs. B-DNA || -----
# 1. Prepare G4 data as a list for fast access
g4_dt <- as.data.table(g4hunter_df)
colnames(g4_dt) <- c("chr","start",'end','length','score','sample')
g4_list <- split(g4_dt, by = "sample")

# 2. Identify samples present in both Methylation (ldf) and G4 data
sample_list <- intersect(names(ldf), names(g4_list))

# 3. Process sample-by-sample
g4_meth_summary_list <- lapply(sample_list, function(s) {
  
  # Convert current sample's methylation dataframe to data.table
  m_sub <- as.data.table(ldf[[s]])
  
  # Ensure standard column names (Assuming: chr, start, end, meth)
  setnames(m_sub, old = colnames(m_sub)[1:4], new = c("chr", "start", "end", "meth"))
  
  # Get G4s for this sample
  g_sub <- g4_list[[s]]
  
  # Set keys for fast overlap
  setkey(m_sub, start, end)
  setkey(g_sub, start, end)
  
  # Overlap methylation sites with G4 regions
  ov <- foverlaps(m_sub, g_sub, type = "any")
  
  # Categorize: If 'score' is NA, it didn't hit a G4 region
  ov[, category := ifelse(is.na(score), "B-DNA", "G4-DNA")]
  
  # 4. Summarize: Fix the 's' length issue by assigning it after aggregation
  # We aggregate by category only first
  summary <- ov[, .(
    avg_meth = mean(meth, na.rm = TRUE),
    median_meth = median(meth, na.rm = TRUE),
    avg_g4_stability = mean(score, na.rm = TRUE),
    site_count = .N
  ), by = .(category)]
  
  # Now add the sample ID to the result (broadcasts 's' to the 2 resulting rows)
  summary[, sample := s]
  
  return(summary)
})

# 5. Combine and Plot
final_g4_meth_summary <- rbindlist(g4_meth_summary_list)

# Visualization
ggplot(final_g4_meth_summary, aes(x = category, y = avg_meth, fill = category)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_fill_manual(values = c("B-DNA" = "#999999", "G4-DNA" = "#e63946")) +
  theme_minimal() +
  labs(
    title = "Methylation Levels: B-DNA vs G4-DNA",
    subtitle = "Aggregated sample-by-sample using individual G4 profiles",
    y = "Average Methylation",
    x = "Genomic DNA Category"
  )
### ----- || Analysis: methylation vs. stability in G4s || -----
# 1. Filter to keep only the G4-DNA category (where stability exists)
g4_stability_plot_data <- final_g4_meth_summary[category == "G4-DNA"]

# 2. Build the 2D Scatterplot
ggplot(g4_stability_plot_data, aes(x = avg_g4_stability, y = avg_meth)) +
  # Add points colored by sample (optional)
  geom_point(size = 3, alpha = 0.7, color = "#e63946") +
  # Add a linear regression trend line
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Correlation: G4 Stability vs. Methylation",
    subtitle = "Each point represents a sample average",
    x = "Average G4 Stability (G4Hunter Score)",
    y = "Average Methylation Level"
  ) +
  # Add Pearson correlation coefficient to the plot
  annotate("text", 
           x = min(g4_stability_plot_data$avg_g4_stability), 
           y = max(g4_stability_plot_data$avg_meth), 
           label = paste("r =", round(cor(g4_stability_plot_data$avg_g4_stability, 
                                          g4_stability_plot_data$avg_meth), 2)),
           hjust = 0, vjust = 1, fontface = "italic")
### ----- || Analysis: G4 methylation in different sequence classes|| -----
meth_master <- rbindlist(ldf, idcol = "sample")
sample_list <- unique(meth_master$sample)
g4_dt <- as.data.table(g4hunter_df)
colnames(g4_dt) <- c("chr","start",'end','length','score','sample')
g4_list <- split(g4_dt, by = "sample")
# 1. Process sample-by-sample (Lightweight version)
comparison_list <- lapply(sample_list, function(s) {
  
  # A. Prepare DataTables - Select only needed columns immediately
  m_sub <- as.data.table(ldf[[s]])[, 1:4, with=FALSE]
  setnames(m_sub, c("chr", "start", "end", "meth"))
  
  g_sub <- g4_list[[s]][, .(start, end, score)]
  ann_sub <- as.data.table(ldf_ann[[s]])[, .(start, end, sequence_class)]
  
  # B. First Overlap: Distinguish B-DNA from G4-DNA
  setkey(m_sub, start, end)
  setkey(g_sub, start, end)
  meth_overlap <- foverlaps(m_sub, g_sub, type = "any")
  
  # Define DNA Type and clean up columns to save RAM
  meth_overlap[, dna_type := ifelse(is.na(score), "B-DNA", "G4-DNA")]
  meth_overlap <- meth_overlap[!is.na(start), .(chr, start, end, meth, dna_type)]
  
  # C. Second Overlap: Assign to Sequence Classes
  setkey(meth_overlap, start, end)
  setkey(ann_sub, start, end)
  final_ov <- foverlaps(meth_overlap, ann_sub, type = "any")
  
  # D. Group sequence classes efficiently
  final_ov[, class_group := case_when(
    grepl("PAR", sequence_class) ~ "PAR",
    grepl("XDR", sequence_class) ~ "XDR",
    grepl("HET", sequence_class) ~ "HET",
    grepl("^P",    sequence_class) ~ "P",
    grepl("CEN", sequence_class) ~ "CEN",
    TRUE ~ "Other"
  )]
  
  # E. Summarize IMMEDIATELY - This is the "lightweight" part
  # Each sample now returns at most 12 rows (6 classes * 2 DNA types)
  summary <- final_ov[, .(
    avg_meth = mean(meth, na.rm = TRUE),
    site_count = .N
  ), by = .(class_group, dna_type)]
  
  summary[, sample := s]
  return(summary)
})

# 2. Combine results - This table will be very small and fast to plot
final_comparison_summary <- rbindlist(comparison_list)

# 3. Generate the Side-by-Side Boxplot
plot_data <- final_comparison_summary[!is.na(avg_meth)]
plot_data$class_group <- factor(plot_data$class_group, levels = c("PAR", "XDR", "HET", "P", "CEN", "Other"))

ggplot(plot_data, aes(x = class_group, y = avg_meth, fill = dna_type)) +
  # Disegna i violin plot affiancati (position_dodge)
  geom_violin(alpha = 0.8, position = position_dodge(width = 0.8), draw_quantiles = 0.5) +
  # Aggiunge i punti individuali sopra i violini
  geom_jitter(aes(group = dna_type), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1), 
              alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = c("B-DNA" = "#999999", "G4-DNA" = "#e63946")) +
  theme_minimal() +
  labs(
    title = "Methylation: G4-DNA vs B-DNA per Sequence Class",
    subtitle = "Violin plot showing sample distribution",
    x = "Sequence Class Group",
    y = "Average Methylation Level",
    fill = "DNA Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### ----- || Analysis: methlyation by nonB-DNA per sequence class || -----
sample_list <- unique(meth_master$sample)
meth_class_summary_list <- lapply(sample_list, function(s) {
  # Subset data
  m_sub <- meth_master[sample == s]
  n_sub <- nonb_master[sample == s]
  c_sub <- classes_master[sample == s]
  
  # Set keys for the "right" tables
  setkey(n_sub, sample, start, end)
  setkey(c_sub, sample, start, end)
  
  # --- FIRST OVERLAP: Methylation onto nonB-DNA ---
  ov_nonb <- foverlaps(m_sub, n_sub, 
                       by.x = c("sample", "start", "end"), 
                       by.y = c("sample", "start", "end"), 
                       type = "any")
  
  # FIX: foverlaps created i.start and i.end. 
  # We must remove them before the NEXT overlap to avoid name conflicts.
  ov_nonb_only <- ov_nonb[!is.na(dna_type)]
  ov_nonb_only[, `:=`(i.start = NULL, i.end = NULL)] 
  
  # --- SECOND OVERLAP: Map results to Sequence Classes ---
  ov_final <- foverlaps(ov_nonb_only, c_sub, 
                        by.x = c("sample", "start", "end"), 
                        by.y = c("sample", "start", "end"), 
                        type = "any")
  
  # Cleanup class names
  ov_final[is.na(sequence_class), sequence_class := "Other"]
  
  # Summarize immediately
  summary <- ov_final[, .(
    avg_meth = mean(meth, na.rm = TRUE),
    site_count = .N
  ), by = .(sample, dna_type, sequence_class)]
  
  return(summary)
})
# Combine and Filter
meth_class_agg <- rbindlist(meth_class_summary_list)[site_count >= 10]
# --- PLOTTING ---
# Assuming meth_class_agg is already created from your lapply
# 1. Broad Categorization
meth_class_agg[, plot_class := fcase(
  grepl("XTR", sequence_class), "X-transp",
  grepl("XDR", sequence_class), "X-deg",
  grepl("^P[1-8]$", sequence_class), "palindromic",
  grepl("PAR", sequence_class), "PAR",
  default = "other"
)]
# Set factor levels for consistent plotting order
meth_class_agg$plot_class <- factor(meth_class_agg$plot_class, 
                                    levels = c("X-deg", "inverted", "palindromic", "PAR", "other"))
# 2. Subset for the Palindrome Focus (Density Plot)
# Plot B: Palindrome Density Profiles (Stacked)
p2_density <- ggplot(pal_data, aes(x = avg_meth, fill = dna_type)) +
  geom_density(alpha = 0.5) +
  # Use facet_grid to stack them vertically by subclass
  facet_grid(sequence_class ~ ., scales = "free_y") + 
  theme_minimal() +
  labs(
    title = "Methylation Density by Palindrome Sub-class",
    x = "Average Methylation (%)",
    y = "Density",
    fill = "nonB-DNA Type"
  ) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0), # Keep P1-P8 labels horizontal
    panel.spacing = unit(0.2, "lines")      # Tighten space between plots
  )

print(p2_density)

                    
