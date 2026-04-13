# ------------- || LIBRARIES || -----------------
library(tibble)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(phylogram)
library(ape)
library(zoo)
library(RColorBrewer)
library(scales)
library(purrr)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(corrplot)
library(data.table)
library(Gviz)
library(rtracklayer)
library(plyranges)
library(ggridges) 
library(patchwork) 
library(valr)
library(circlize)
library(magrittr)
library(scales)
library(ComplexHeatmap)
library(grid)
library(magick)
library(Biostrings)
library(stringr)

# ------------ || DATA UPLOAD || -------------
## methylation score data
#filenames = list.files(path = "Work/Data/BedGraph/BedGraph_Pille_Hallast/m_C0G", full.names = T)
###filenames = filenames[file.size(filenames) > 0]
ldf = lapply(filenames,read.table)
## annotation
#filenames_ann = list.files(path = "Work/annotation/seq_classes_Pille_Hallast/2026-01_final/t2tv2/", full.names = T)
ldf_ann = lapply(filenames_ann,read.table)

# ------------- || FILE ORGANIZATION || ---------------
# 1. Define Paths
input_path <- "Work/annotation/seq_classes_Pille_Hallast/2026-01_final/t2tv2/"
output_path <- "Work/annotation/seq_classes_Pille_Hallast/2026-01_final/t2tv2/new/" # New folder name
# 2. Create the new folder if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
# 3. Load files and annotation
ann = read.table("Work/annotation/seq_classes_Pille_Hallast/ann_colors_Y_chrom_pangenome.txt", 
                 header = TRUE, stringsAsFactors = FALSE)
# 4. Process the data (Vectorized match is faster than a for-loop with grep)
ldf_ann_2 = lapply(ldf_ann, function(df) {
  # Map colors based on matching V4 to ann$region
  df$color <- ann$color[match(df$V4, ann$region)]
  
  # Add BedGraph specific columns
  df$thickStart = df$V2
  df$thickEnd = df$V3
  
  # Order columns: V1, V2, V3, V4, V5, V6, thickStart, thickEnd, color
  return(df[, c("V1","V2","V3","V4","V5","V6","thickStart","thickEnd","color")])
})
# 5. Define new file paths in the new folder
# basename() extracts just "file.bed" from "Downloads/BedGraph/file.bed"
new_filenames <- file.path(output_path, basename(filenames_ann))
# 6. Save files to the new folder
# Use the length-safe check to avoid the "longer argument" warning
if (length(ldf_ann_3) == length(new_filenames)) {
  mapply(function(data, file_path) {
    write.table(
      x = data,
      file = file_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )
  }, ldf_ann_2, new_filenames)
} else {
  stop("Mismatch between number of processed dataframes and filenames.")
}
## -------- || take out ERR regions from ldf_ann || ----------
### Function to remove rows where sequence_class is 'ERR'
remove_err_rows <- function(df) {
  df[df$V4 != "ERRBASE", ]
}
### Apply the function to all dataframes in the list
ldf_ann_2 <- lapply(ldf_ann, remove_err_rows)
ldf_ann = ldf_ann_2
rm(ldf_ann_2)
### Alternatively, using an anonymous function directly within lapply:
### my_list_filtered <- lapply(my_list, function(df) df[df$sequence_class != "ERR", ])
## -------- || change seqnames in BedGraph files || ----------
## get sample names
samples = gsub(".*[//]|[_].*$", "", filenames)
samples_ann = gsub(".*[//]|[.].*$", "", filenames_ann)
names(ldf) = samples
names(ldf_ann) = samples_ann
## change seqnames in BedGraph files
ldf = lapply(ldf, function(df) {
  df[[1]] = gsub(".*[_]", "", df[[1]])
  return(df)
})
## ---------- || rename colnames of ldf and ldf_ann || -------------------
### filter for chrom, start, end, score (as %) from the list of bedgraph files obtained after calling
# Create the filtered list from the raw list
ldf_2 <- lapply(ldf, function(df) {
  # 1. Select specific columns (chrom, start, end, score)
  # Score is column 11 in your bedgraph structure
  new_df <- df[, c(1, 2, 3, 5, 11, 12)]
  # 2. Rename the columns
  colnames(new_df) <- c("chrom", "start", "end", "total_bases","score", "mod_bases")
  return(new_df)
})
ldf_3 <- lapply(ldf, function(df) {
  # 1. Select specific columns (chrom, start, end, score)
  # Score is column 11 in your bedgraph structure
  new_df <- df[, c(1, 2, 3, 11)]
  # 2. Rename the columns
  colnames(new_df) <- c("chrom", "start", "end", "score")
  return(new_df)
})
ldf = ldf_3
rm(ldf_3)

### --- 2. Rename ALL ldf_ann dataframes using index position ---
rename_ldf_ann_by_index <- function(df) {
  if (ncol(df) < 9) {
    # If fewer than 9 cols (e.g. standard bed format), just name the first 4
    df <- df %>% rename(chrom = names(df)[1], start = names(df)[2], end = names(df)[3], class = names(df)[4])
  } else {
    # If 9 columns are present:
    df <- setNames(df, c("chrom", "start", "end", "class", "score", "strand", "thickStart", "thickEnd", "color"))
  }
  return(df)
}
ldf_ann <- lapply(ldf_ann, rename_ldf_ann_by_index)
### Verification step:
print(names(ldf[[1]]))
print(names(ldf_ann[[1]]))
## ----- || Clean Chromosome Names Explicitly || -----
# We target the first column directly to ensure we don't touch other data
for (i in seq_along(ldf)) {
  ldf[[i]][[1]] <- gsub("^.*_chrY$", "chrY", as.character(ldf[[i]][[1]]))
  ldf_2[[i]][[1]] <- gsub("^.*_chrY$", "chrY", as.character(ldf_2[[i]][[1]]))
  ldf_ann[[i]][[1]] <- gsub("^.*_chrY$", "chrY", as.character(ldf_ann[[i]][[1]]))
}

# ------------- || CHECK VARIATIONS OF COVERAGE IN THE BEDGRAPH FILES || -------------
## ----- || check coverage on P7 for each individual || -----
common_names <- intersect(names(ldf_2), names(ldf_ann))
# 2. Subset lists so they match perfectly
ldf_2_subset <- ldf_2[common_names]
ldf_ann_subset <- ldf_ann[common_names]
# 3. Proceed with map2_dfr now that sizes match
plot_df <- map2_dfr(ldf_2_subset, ldf_ann_subset, function(df_data, df_ann) {
  
  # Filter coverage (total_bases) for rows where sequence_class is "P7"
  # Note: Ensure the row count matches within each dataframe pair
  p7_coverage <- df_data$total_bases[df_ann$sequence_class == "P7"]
  
  return(data.frame(coverage = p7_coverage))
  
}, .id = "Individual")
# 4. Create the density plot
# Use the same 'plot_df' created in your previous steps
ggplot(plot_df, aes(x = coverage, fill = Individual)) +
  # Add alpha (transparency) so you can see overlapping distributions
  geom_density(alpha = 0.4) + 
  theme_minimal() +
  labs(
    title = "P7 Coverage Distribution by Individual",
    x = "Total Bases (Coverage)",
    y = "Density",
    fill = "Individual"
  ) +
  theme(legend.position = "none") # Re-enable legend to identify individuals
## ----- || check coverage on specific samples (HG02055) with Gviz || -----
# 1. Extract data for HG02055
target_id <- "HG02055"
df_data <- ldf_2[[target_id]]
df_ann  <- ldf_ann[[target_id]]
# 2. Filter for sequence class "P7"
p7_mask <- df_ann$sequence_class == "P7"
p7_ranges <- df_data[p7_mask, ]
# 3. Convert to GRanges (Assumes your dataframes have chr, start, end)
# Adjust column names if yours are different
gr_p7 <- GRanges(
  seqnames = p7_ranges$chr,
  ranges = IRanges(start = p7_ranges$start, end = p7_ranges$end),
  coverage = p7_ranges$total_bases
)
# 4. Initialize Gviz Tracks
# Genome Axis shows the genomic coordinates
axis_track <- GenomeAxisTrack()
# DataTrack with polygon type for the coverage
data_track <- DataTrack(
  range = gr_p7,
  data = "coverage",
  name = paste("P7 Coverage:", target_id),
  type = "polygon",        # Fills the area under the curve
  fill.mountain = "orange", # Color of the polygon fill
  col.mountain = "red",     # Color of the polygon border
  baseline = 0
)
# 5. Plot the tracks
# You can use 'from' and 'to' to zoom into a specific high-variability region
plotTracks(
  list(axis_track, data_track),
  main = paste("Coverage Polygon for", target_id),
  cex.main = 1.2
)
## ----- || check MAPQ along P7 in HG02055 || -----
data_mapq = read.csv("Work/Data/HG02055_chrY.windows.csv", sep = ",", stringsAsFactors = F, header = T)
data_mapq_P7 = data_mapq[which(data_mapq$Start>=17856764 & data_mapq$End<=17886841),]
boxplot(data_mapq_P7$Mean_MAPQ)

# ------------- || METHYLATION AVERAGE ACROSS SEQUENCE CLASSES || -------------
## create violin plots showing average methylation in palindromic sequences vs. everything else
## Assuming your list of data frames is stored in a variable called 'my_bed_list'
# 1. Clean Chromosome Names Explicitly
process_individual_data_frame = function(sample_id) {
  # 1. Isolate and Clean Methylation Data
  meth_df <- ldf[[sample_id]][, 1:4]
  colnames(meth_df) <- c("chr", "temp_start", "temp_end", "meth")
  
  # 2. Isolate and Clean Annotation Data 
  annot_df <- ldf_ann[[sample_id]][, 1:4]
  colnames(annot_df) <- c("chr", "temp_start", "temp_end", "class")
  
  # 3. Create GRanges
  meth_gr = makeGRangesFromDataFrame(meth_df, seqnames.field="chr", 
                                     start.field="temp_start", end.field="temp_end",
                                     keep.extra.columns=TRUE)
  
  annot_gr = makeGRangesFromDataFrame(annot_df, seqnames.field="chr", 
                                      start.field="temp_start", end.field="temp_end",
                                      keep.extra.columns=TRUE)
  
  # 4. Standard check for chromosome levels
  common_chrs <- intersect(seqlevels(meth_gr), seqlevels(annot_gr))
  if(length(common_chrs) == 0) return(NULL) 
  
  meth_gr <- keepSeqlevels(meth_gr, common_chrs, pruning.mode="tidy")
  annot_gr <- keepSeqlevels(annot_gr, common_chrs, pruning.mode="tidy")
  
  # 5. Overlap and Extract
  overlaps = findOverlaps(meth_gr, annot_gr, type = "any")
  if (length(overlaps) == 0) return(NULL)
  
  # 6. Build Tibble including coordinates
  return(tibble(
    chr   = as.character(seqnames(meth_gr))[queryHits(overlaps)],
    start = start(meth_gr)[queryHits(overlaps)],  # Extracts the start coordinate
    end   = end(meth_gr)[queryHits(overlaps)],    # Extracts the end coordinate
    meth  = mcols(meth_gr)$meth[queryHits(overlaps)],
    class = mcols(annot_gr)$class[subjectHits(overlaps)],
    individual_id = sample_id
  ))
}
# 3. Combine and Plot
all_individuals_data = map_dfr(samples, process_individual_data_frame)
# Check if data exists before plotting
if (nrow(all_individuals_data) == 0) {
  stop("No overlaps found. Check your chromosome names with: 
       head(ldf[[1]][,1]) and head(ldf_ann[[1]][,1])")
}
plot_data = all_individuals_data %>%
  mutate(
    group = case_when(
      grepl("P[1-8]", class) ~ "palindromic",
      grepl("IR[1-4]", class) ~ "inverted",
      grepl("XDR", class) ~ "X-deg",
      grepl("PAR", class) ~ "PAR", 
      TRUE ~ "other"
    ),
    group = factor(group, levels = c("X-deg", "inverted", "palindromic", "PAR", "other"))
  ) %>%
  group_by(individual_id, group) %>%
  # FIX: Use 'meth' here because that is what the function returns
  summarise(avg_methylation_per_group = mean(meth, na.rm = TRUE), .groups = "drop")
# 4. Plotting
ggplot(plot_data, aes(x = group, y = avg_methylation_per_group, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, fill = "white") +
  theme_minimal()
### ----- || Methylation in all palindromes || ------
# 2. Filter for palindromes and calculate mean PER REGION
# Corrected Pattern Logic
palindrome_region_data <- all_individuals_data %>%
  filter(grepl(paste(c("P[1-8]", "red"), collapse = "|"), class)) %>%
  group_by(individual_id, class) %>%
  summarise(avg_meth = mean(meth, na.rm = TRUE), .groups = "drop")
palindrome_region_data$class[grep("P1",palindrome_region_data$class)] = "P1"
palindrome_region_data$class[grep("P2",palindrome_region_data$class)] = "P2"
palindrome_region_data$class[grep("red",palindrome_region_data$class)] = "P2"
palindrome_region_data$class[grep("P3",palindrome_region_data$class)] = "P3"
palindrome_region_data$class[grep("P5",palindrome_region_data$class)] = "P5"
# 3. Create the Density Plot
ggplot(palindrome_region_data, aes(x = avg_meth, fill = class, color = class)) +
  # alpha sets transparency (0.3 - 0.5 is usually best for overlaps)
  geom_density(alpha = 0.4, size = 0.8) +
  # Standardize X-axis for methylation (0 to 1)
  # Optional: Define your own colors if you have many classes
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Methylation Density Overlay: Palindromic Regions",
    subtitle = "Comparing distribution of P1 through P8 on one axis",
    x = "Average Methylation Level",
    y = "Density",
    fill = "Palindrome Class",
    color = "Palindrome Class"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

