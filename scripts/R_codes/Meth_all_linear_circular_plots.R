## ------ || METHYLATION ON ALL INDIVIDUALS - THE LINEAR PLOT || ------
# 0. Prepare Genome Skeleton (Sequence Classes)
seq_classes <- c("PAR1", "XDR1", "XTR1", "IR3", "AMPL", "XTR2", "XDR2", "IR1", "AMPL", "IR4", 
                 "IR3", "SAT", "CEN", "SAT", "DYZ17", "XDR3", "AMPL", "P8", "XDR4", "AMPL", 
                 "XDR5", "AMPL", "P6", "XDR6", "AMPL", "P5", "P4", "XDR7", "DYZ19", "XDR8", 
                 "AMPL", "IR4", "P3", "IR1", "P2", "P2", "P1", "IR1-g2", "IR1-b3", "P1", 
                 "IR1-g3", "IR1-b4", "HET", "PAR2", "XTR", "PAR2")
metadata = read.csv("Work/annotation/sample_metadata.csv", sep = ",", stringsAsFactors = F, header = T)
ordered_metadata <- metadata %>% 
  dplyr::filter(ID %in% names(ldf)) %>% 
  dplyr::arrange(haplogroup_medium, ID)
# 1. Load Phylogeny Order
phylog_annotation <- read.delim("Work/annotation/HPRC_HGSVC3_sample_annotations_20March2025.txt", stringsAsFactors = F, header = T, sep = "\t")
phylo_order <- phylog_annotation[,2] # Your IDs in phylogenetic order
# Load the file
nexus_data <- read.nexus("Work/annotation/142males_HGSVC_HPRC_CEPH_241125_150M-FOR_SHARING.nex")
# Check if it's a single tree (phylo) or multiple (multiPhylo)
if (inherits(nexus_data, "multiPhylo")) {
  # If multiple, pick the first one or a consensus tree
  my_tree <- nexus_data[[1]] 
} else {
  my_tree <- nexus_data
}
# 1. Define Resolution (5000 relative bins across the Y)
n_bins <- 5000 
common_samples = intersect(names(ldf), names(ldf_ann))
common_samples <- intersect(common_samples, phylog_annotation$sample)
# 2. Process each sample - Data Frame Version
meth_list <- map(common_samples, function(s) {
  meth_df <- ldf[[s]]     # This is a Data Frame
  ann_df  <- ldf_ann[[s]]  # This is a Data Frame
  
  # 1. Safety check for empty or missing data
  if (is.null(meth_df) || nrow(meth_df) == 0 || is.null(ann_df) || nrow(ann_df) == 0) {
    return(rep(NA, n_bins))
  }
  
  # 2. Calculate relative positions
  # Accessing columns directly by name or index since it's a data frame
  max_len <- max(ann_df$end, na.rm = TRUE)
  rel_pos <- meth_df$start / max_len
  
  # 3. Assign each site to a bin (1 to 5000)
  bin_idx <- cut(rel_pos, breaks = seq(0, 1, length.out = n_bins + 1), 
                 labels = FALSE, include.lowest = TRUE)
  
  # 4. Extract Methylation (4th column)
  # Force to a plain numeric vector to avoid any hidden attributes
  meth_values <- as.numeric(as.vector(meth_df[[4]]))
  
  # 5. Aggregate per bin
  bin_values <- tapply(meth_values, bin_idx, mean, na.rm = TRUE)
  
  # 6. Fill the 5000-element vector
  full_vec <- rep(NA, n_bins)
  if (length(bin_values) > 0) {
    indices <- as.numeric(names(bin_values))
    # Safety: ensure indices are within 1:5000
    valid <- indices >= 1 & indices <= n_bins
    full_vec[indices[valid]] <- bin_values[valid]
  }
  
  return(full_vec)
})

meth_matrix <- do.call(rbind, meth_list)
rownames(meth_matrix) <- common_samples

# 3. SYNC TREE AND MATRIX 
# This is the most important step: keep only samples present in both
common_ids <- intersect(rownames(meth_matrix), my_tree$tip.label)
meth_matrix <- meth_matrix[common_ids, ]
my_tree <- keep.tip(my_tree, common_ids)

# Convert tree to a dendrogram format for Heatmap
dend <- as.dendrogram(my_tree)
# 4. Create Top Annotation (Using first sample as structural reference)
ref_id <- common_samples[1] 
ref_ann <- ldf_ann[[ref_id]]
ref_ann = ref_ann[which(ref_ann$chr=="chrY"),]
# 2. Map the 5000 bins to the Sequence Classes
# This creates the vector needed for the inner ring
bin_coords <- seq(min(ref_ann$start), max(ref_ann$end), length.out = n_bins)
idx_class  <- findInterval(bin_coords, ref_ann$start)
idx_class[idx_class == 0] <- 1
bin_classes <- ref_ann$sequence_class[idx_class]# Map the 5000 bins to the Sequence Classes of the reference

# 5. Final "Pretty" Plot
# Example: Creating a simple annotation based on column names or categories
top_ann <- HeatmapAnnotation(
    Group = anno_block(gp = gpar(fill = "grey")), # Or map to a vector
    show_annotation_name = TRUE
)
#pdf("Normalized_Phylo_Heatmap_Final.pdf", width = 18, height = 12)
h = Heatmap(meth_matrix, 
        name = "Methylation", 
        
        # --- Tree Integration ---
        cluster_rows = dend,           # Use your .nex tree for the Y-axis
        row_dend_side = "left",        # Place the tree on the left
        row_dend_width = unit(4, "cm"), 
        
        # --- Aesthetics ---
        cluster_columns = FALSE,       # Keep genomic order on X-axis
        top_annotation = top_ann,
        col = colorRamp2(c(0, 50, 100), c("#f7fbff", "#6baed6", "#08306b")), # Prettier Blue gradient
        
        # --- Clean up ---
        show_row_names = FALSE,        # Labels usually overlap at 142 samples
        show_column_names = FALSE,
        use_raster = TRUE,             # Keeps file size small/fast
        raster_quality = 5,
        
        # Optional: Split by continent if you want the gaps shown in your image
        # row_split = ordered_metadata$continent[match(common_ids, ordered_metadata$ID)]
)
draw(h)
#dev.off()
## ------ || METHYLATION ON ALL INDIVIDUALS - CIRCULAR VERSION || ------
# 1. Sync tree and methylation matrix
common_ids <- intersect(rownames(meth_matrix), my_tree$tip.label)
my_tree <- keep.tip(my_tree, common_ids)
meth_matrix <- meth_matrix[common_ids, , drop = FALSE]

# Add column names if missing
if (is.null(colnames(meth_matrix))) {
  colnames(meth_matrix) <- paste0("bin_", seq_len(ncol(meth_matrix)))
}

n_bins <- ncol(meth_matrix)

# 2. Build annotation bins from one reference sample
ref_id <- common_ids[1]
ref_ann <- ldf_ann[[ref_id]] %>%
  as.data.frame() %>%
  dplyr::filter(chr == "chrY") %>%
  dplyr::arrange(start)

ref_max <- max(ref_ann$end, na.rm = TRUE)
bin_coords <- seq(0, ref_max, length.out = n_bins)

idx_class <- findInterval(bin_coords, ref_ann$start)
idx_class[idx_class == 0] <- 1
idx_class[idx_class > nrow(ref_ann)] <- nrow(ref_ann)

bin_classes <- ref_ann$sequence_class[idx_class]

# 3. Long-format methylation dataframe
meth_df <- as.data.frame(meth_matrix) %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    cols = -sample,
    names_to = "bin",
    values_to = "methylation"
  ) %>%
  mutate(
    sample = factor(sample, levels = my_tree$tip.label),
    bin_num = as.numeric(gsub("bin_", "", bin))
  )

# 4. Single annotation spoke dataframe
#    This will place the annotation at the angle of one chosen sample.
#    Change ann_tip to move the spoke around the circle.
ann_tip <- my_tree$tip.label[1]

annotation_spoke_df <- data.frame(
  sample = ann_tip,
  bin_num = seq_len(n_bins),
  sequence_class = bin_classes,
  stringsAsFactors = FALSE
) %>%
  mutate(
    sample = factor(sample, levels = my_tree$tip.label)
  )

# 5. Colors
class_colors <- c(
  PAR1     = "#f4a261",
  XDR1     = "#e76f51",
  XTR1     = "#ffbe0b",
  IR3      = "#90be6d",
  AMPL     = "#c77dff",
  XTR2     = "#ffd166",
  XDR2     = "#d62828",
  IR1      = "#43aa8b",
  IR4      = "#4d908e",
  SAT      = "#adb5bd",
  CEN      = "#6c757d",
  DYZ17    = "#457b9d",
  P8       = "#3a86ff",
  XDR4     = "#f72585",
  XDR5     = "#b5179e",
  P6       = "#7209b7",
  XDR6     = "#560bad",
  P5       = "#480ca8",
  P4       = "#4361ee",
  XDR7     = "#4895ef",
  DYZ19    = "#4cc9f0",
  XDR8     = "#2a9d8f",
  P3       = "#8ac926",
  P2       = "#ff595e",
  P1       = "#1982c4",
  `IR1-g2` = "#a8dadc",
  `IR1-b3` = "#bde0fe",
  `IR1-g3` = "#caffbf",
  `IR1-b4` = "#ffc6ff",
  HET      = "#9d4edd",
  PAR2     = "#f28482",
  XTR      = "#f6bd60"
)

present_classes <- unique(annotation_spoke_df$sequence_class)
class_colors <- class_colors[names(class_colors) %in% present_classes]

# 6. Circular tree - Use xlim to shrink the tree's relative size
# Increasing the max xlim (e.g., to 5 or 10) makes the tree occupy less radial space
p <- ggtree(my_tree, layout = "circular", size = 0.2) + 
  xlim(0, 5) 

# 7. Methylation annulus - Increased pwidth
p1 <- p +
  geom_fruit(
    data = meth_df,
    geom = geom_tile,
    mapping = aes(
      y = sample,
      x = bin_num,
      fill = methylation
    ),
    offset = 0.02,
    pwidth = 0.8, # Increased from 0.42 to 0.8 to make it "bigger"
    tilewidth = 1,
    tileheight = 0.95
  ) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#6baed6", "#08306b"),
    limits = c(0, 100),
    na.value = "white",
    name = "Meth"
  )

# 8. Reset fill scale
p1 <- p1 + ggnewscale::new_scale_fill()

# 9. Final Plot with expansion settings
p2 <- p1 +
  geom_fruit(
    data = annotation_spoke_df,
    geom = geom_tile,
    mapping = aes(
      y = sample,
      x = bin_num,
      fill = sequence_class
    ),
    offset = 0.1,
    pwidth = 0.3, 
    tilewidth = 1,
    tileheight = 0.95
  ) +
  scale_fill_manual(
    values = class_colors,
    na.value = "white",
    name = "ChrY annotation"
  ) +
  # --- THE "MAKE IT BIGGER" SETTINGS ---
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.x = unit(0.1, 'cm'),
    # This removes the empty space around the circle
    plot.margin = margin(-1, -1, -1, -1, "cm") 
  ) +
  # Forces the coordinate system to expand to the edges
  coord_polar(theta = "y", start = 0, clip = "off")

# To see it "Big", it's best to save it directly to a file:
ggsave("ChrY_Methylation_Plot.pdf", plot = p2, width = 12, height = 10, device = "pdf")
