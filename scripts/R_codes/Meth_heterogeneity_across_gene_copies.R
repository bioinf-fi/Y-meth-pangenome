## ------ || METHYLATION HETEROGENEITY ACROSS GENE COPIES || ------
# ASSUME YOUR INPUTS ARE LOADED AS:
# 1. methylation_heterogeneity_list: A named list of data frames
#    (e.g., list(HG01433 = df1, ...), with 'genes' as rownames/a column, and 
#     a column for heterogeneity value, e.g., 'heterogeneity_score').
# 2. gene_copies_wide: A data frame with individuals as rows and genes as columns.
#    (e.g., columns are 'GeneA', 'GeneB', etc., and rownames/a column called 'individual' 
#     stores the individual names).
## ----- || Upload data || -----
## meth. heterogeneity - genes
filenames_het_ampl = list.files(path = "Work/Pille_Hallast_methylation/modbamtools/AmplGenes/", full.names = T)
ldf_het_ampl = lapply(filenames_het_ampl,read.table)
## ------ || organize ldf_het_ampl || ----------
samples_het_ampl = gsub(".*[//]|[_].*$", "", filenames_het_ampl)
names(ldf_het_ampl) = samples_het_ampl
### Function to remove rows where sequence_class is 'ERR'
remove_err_rows <- function(df) {
  df[df$V4 != "ERR", ]
}
### Apply the function to all dataframes in the list
ldf_het_2 <- lapply(ldf_het_ampl, remove_err_rows)
ldf_het_ampl = ldf_het_2
rm(ldf_het_2)
# ----- || Step 1: Process Methylation Data || -----
# Combine the list of data frames into one, automatically adding the list name 
# as a new column named 'individual'.
### Define the function to process each data frame and add the sample ID
process_het_file <- function(df, sample_id) {
  df$individual_id <- sample_id
  df <- df %>%
    rename(
      gene = V4, 
      score = V6,
      start = V2,
      end = V3
    ) %>%
    select(individual_id, gene, score, start, end)
  return(df)
}
het_long_ampl <- map2_dfr(ldf_het_ampl, samples_het_ampl, process_het_file) 
# Keep only necessary columns
# ----- || Step 2: Process Gene Copies Data || -----
# Reshape the wide format gene copies data into a long format.
gene_copies_long <- gene_copies %>%
  # Assuming individuals are in a column named 'individual'
  pivot_longer(
    cols = -Sample.ID , # Reshape all columns except 'individual'
    names_to = "gene",  # Put gene names into a new 'gene' column
    values_to = "n.copies" # Put copy numbers into a new 'n.copies' column
  )
colnames(gene_copies_long)[1] = "individual_id"
# ----- || Step 3: Merge the two data frames || -----
# Perform a full join on both 'individual' and 'gene' columns.
# This creates your final desired data frame.
final_merged_data <- full_join(
  het_long_ampl, 
  gene_copies_long, 
  by = c("individual_id", "gene")
)
# View the result
View(final_merged_data)
### ----- || Part 4: Plotting || -----
##### plot methylation heterogeneity in all genes
ggplot(final_merged_data, aes(x = gene, y = score, fill = gene)) +
geom_boxplot() +
  theme_minimal()
##### plot methylation heterogeneity across gene copies
#final_merged_data = na.omit(final_merged_data)
ggplot(data = final_merged_data, aes(x = n.copies, y = score)) +
  geom_point(alpha = 0.5) +
  # Option A: Add a linear trend line
  geom_smooth(method = "gam", color = "blue", se = TRUE, na.rm = TRUE) +
  # Option B: Add a curved (LOESS/GAM) trend line (default)
  # geom_smooth(method = "auto", color = "red", na.rm = TRUE) +
  labs(title = "Methylation Heterogeneity vs Gene Copy Number", 
       x = "Number of Gene Copies", 
       y = "Methylation Heterogeneity") +
  theme_minimal()




