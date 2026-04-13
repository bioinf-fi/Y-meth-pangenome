## -------- | METHYLATION HETEROGENEITY ACROSS SEQUENCE CLASSES  | -----------------
## ----- || Upload data || -----
## meth. heterogeneity - seq. classes
filenames_het = list.files(path = "Work/Pille_Hallast_methylation/modbamtools/SeqClasses/", full.names = T)
ldf_het = lapply(filenames_het,read.table)
## ------ || organize ldf_het || ----------
samples_het = gsub(".*[//]|[.].*$", "", filenames_het)
names(ldf_het) = samples_het
### Function to remove rows where sequence_class is 'ERR'
remove_err_rows <- function(df) {
  df[df$V4 != "ERR", ]
}
### Apply the function to all dataframes in the list
ldf_het_2 <- lapply(ldf_het, remove_err_rows)
ldf_het = ldf_het_2
rm(ldf_het_2)
### Alternatively, using an anonymous function directly within lapply:
### my_list_filtered <- lapply(my_list, function(df) df[df$sequence_class != "ERR", ])

## ----- || Methylation heterogeneity across sequence classes || -----
### 1. Define the vectors for the different region categories
palindromic = c("P1","P2","P3","P4","P5","P6","P7","P8")
inverted = c("IR1","IR2","IR3","IR4")
### Assuming these vectors are defined from your existing data
#x_deg = unique(het_long[grep("XDR", het_long$region_group),"region_group"]) 
#PAR = unique(het_long[grep("PAR", het_long$region_group),"region_group"])
### --- Data Preparation (Assuming these steps are working for you) ---
### Define the function to process each data frame and add the sample ID
process_het_file <- function(df, sample_id) {
  df$individual_id <- sample_id
  df <- df %>%
    rename(
      class = V4, 
      score = V10,
      start = V2,
      end = V3
    ) %>%
    select(individual_id, class, score, start, end)
  return(df)
}
### Use map2_dfr to apply the function to all data frames and samples, combining them by rows
het_long <- map2_dfr(ldf_het, samples_het, process_het_file) 
### Assuming 'het_long' is available
### 2. Process the data and add the 'group' factor
plot_data = het_long %>%
  # Use case_when() to assign the correct category based on the 'class' column
  mutate(
    group = case_when(
      # Check for specific memberships using grepl(..., class)
      grepl(paste(palindromic, collapse = "|"), class) ~ "palindromic",
      grepl(paste(inverted, collapse = "|"), class) ~ "inverted",
      # Check for substrings "XDR" or "PAR" using grepl(substring, class)
      grepl("XDR", class) ~ "X-deg", # Renamed to match desired plot label "X-deg"
      grepl("PAR", class) ~ "PAR", 
      TRUE ~ "other" # Default condition for everything else
    )
  ) %>%
  # *** Add this step to order the 'group' column as a factor ***
  mutate(group = factor(group, levels = c("X-deg", "inverted", "palindromic", "PAR", "other"))) %>%
  # Group by individual and by the new 'group' column
  group_by(individual_id, group) %>%
  # Calculate the mean score for each group
  summarise(avg_heterogeneity_per_group = mean(score, na.rm = T),
            .groups = "drop" # Automatically removes grouping after summarizing
  )
### --- Plotting with the correct order ---
ggplot(plot_data, aes(x = group, y = avg_heterogeneity_per_group, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, fill = "white", outlier.shape = NA) +
  scale_y_continuous(limits = c(0,500)) + # Changed ylim() to scale_y_continuous() for better ggplot practice
  theme_minimal() +
  labs(
    x = "Region Category",
    y = "Average Methylation Heterogeneity",
    title = "Heterogeneity in Genomic Regions (Ordered)"
  )
## ------- || METHYLATION HETEROGENEITY ACROSS ALL PALINDROMES || ------------

### The script expects columns named "chrom", "start", "end", "score" in ldf 
### and "chrom", "start", "end", "class" (etc) in ldf_ann. 
### Make sure you use the *penultimate* full code block that handles this mapping correctly.


### 1. Define the exact order of the sequence segments
### Note: This vector should contain every single unique region ID present in your het_long$region_group column,
### in the precise order they appear physically on the Y chromosome map.
### The list provided in the prompt is a good start, but needs to be a complete, exhaustive list of all unique region IDs.
### I will use a representative subset here; replace with your final complete 'order' vector.

sequences <- c("P8","P7", "P6", "P5", "P4", "blue", "teal1",  "teal2", "P3-spacer","green","red", "P2-spacer","gray","yellow")
het_long_new = het_long[which(het_long$class %in% sequences),]

### The complete sequence order is defined here for reference
sequence_class_order <- c("P8", "P8", "P7", "P7", "P6", "P6", "P5", "P5", "P4", "P4", 
                          "blue", "teal1", "P3-spacer", "teal2", "blue", "green", "red", "P2-spacer", "red", 
                          "gray", "blue", "yellow", "green", "red", "red", "green", "yellow", 
                          "blue", "gray")
                          
### 1. Define the mapping for the major palindromic arms and simple numeric sequence index
arm_mapping <- data.frame(
  class = sequence_class_order,
  major_arm = c(
    rep("P8", 2), rep("P7", 2), rep("P6", 2), rep("P5", 2), rep("P4", 2),
    rep("P3", 5),                
    "Green_Region",              
    rep("P2", 3),                
    rep("P1", 10)               
  ),
  sequence_index_simple = seq_along(sequence_class_order)
) %>% distinct(class, major_arm, sequence_index_simple)


### 2. Prepare the raw data by joining with the mapping table
### We assume het_long_new has columns: individual_id, class, score, start, end
het_raw_indexed <- het_long_new %>%
  filter(class %in% sequence_class_order) %>%
  inner_join(arm_mapping, by = "class") %>%
  arrange(major_arm, sequence_index_simple) %>%
  mutate(major_arm = factor(major_arm, levels = c("P8", "P7", "P6", "P5", "P4", "P3", "Green_Region", "P2", "P1")))


### 3. Plot the raw data points using geom_jitter, separated by major_arm using facet_wrap
ggplot(het_raw_indexed, aes(x = sequence_index_simple, y = score)) +
  # Plot all points, colored by individual ID
  geom_jitter(aes(color = individual_id), width = 0.2, alpha = 0.7, size = 1.5) +
  # Add a smoothed line overlay for the trend within each facet
  # We needed to explicitly map x and y aesthetics here:
  geom_smooth(aes(x = sequence_index_simple, y = score), 
              method = "loess", color = "black", size = 1, se = TRUE, inherit.aes = FALSE, 
              data = het_raw_indexed %>% group_by(major_arm, sequence_index_simple) %>% summarise(score = mean(score), .groups = "drop")) +
  # Separate the plots into different panels (facets) for each major arm
  facet_wrap(~ major_arm, scales = "free_x", nrow = 2) + 
  # Set the x-axis labels to the actual region names using a manual break/label lookup
  scale_x_continuous(
    breaks = arm_mapping$sequence_index_simple,
    labels = arm_mapping$class
  ) +
  labs(
    title = "Methylation Heterogeneity Points per Major Palindromic Arm",
    x = "Genomic Region (Sequence Index and Class)",
    y = "Methylation Heterogeneity Score",
    color = "Individual ID"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
