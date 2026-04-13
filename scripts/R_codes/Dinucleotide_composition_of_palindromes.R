## --------- || CHECK DINUCLEOTIDE COMPOSITION OF PALINDROMES || -----------
# 1. Load and process (Normalization Fix included)
data_list <- map_df(fasta_files, function(f) {
  p_id <- str_extract(basename(f), "P[1-8]")
  seqs <- readDNAStringSet(f)
  
  # Calculate raw counts first to handle "N" bases correctly
  counts <- dinucleotideFrequency(seqs)
  
  # Convert to DF and normalize rows so the 16 pairs sum exactly to 1.0
  df <- as.data.frame(counts)
  row_totals <- rowSums(df)
  df_norm <- df / row_totals
  
  df_norm$p_class <- p_id
  return(df_norm)
})

# 2. Summarize for plotting
p_summary <- data_list %>%
  group_by(p_class) %>%
  summarise(across(AA:TT, mean, na.rm = TRUE)) %>%
  pivot_longer(-p_class, names_to = "dinucleotide", values_to = "mean_frequency")

# 3. Define a Pastel Palette (16 distinct pastel colors)
# We use a combination of qualitative palettes to get 16 unique soft colors
muted_palette <- c(brewer.pal(12, "Set3"), brewer.pal(4, "Paired"))

# 4. Final Plot
ggplot(p_summary, aes(x = p_class, y = mean_frequency, fill = dinucleotide)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.1) +
  scale_fill_manual(values = muted_palette) +
  theme_minimal() +
  labs(
    title = "Normalized Dinucleotide Composition (P1-P8)",
    subtitle = "Frequencies rescaled to 100% excluding N-bases",
    x = "Palindrome Class",
    y = "Relative Frequency",
    fill = "Dinucleotide"
  ) +
  theme(
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8),
    panel.grid.major.x = element_blank()
  )
