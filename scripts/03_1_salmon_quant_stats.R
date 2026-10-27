# Fig S1a: Salmon mapping statistics per tissue (values from the MultiQC report)
# Input:  salmon_quant/PRJNA72713/quant_stats_15_tissues.tsv
# Output: p_salmon_quant_stats (ggplot), used by 03_2

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

stats <- read_tsv(file.path(quant_dir, "quant_stats_15_tissues.tsv"), show_col_types = FALSE)
check_tissues(str_replace_all(stats$Tissue, " ", "_"))

# Keep the file order (sorted by M Aligned), top to bottom
plot_data <- stats %>%
  mutate(Tissue = factor(Tissue, levels = rev(Tissue))) %>%
  pivot_longer(cols = -Tissue, names_to = "Metric", values_to = "Value") %>%
  mutate(LabelText = ifelse(Metric == "% Aligned", paste0(sprintf("%.1f", Value), "%"), as.character(Value)),
         Metric = factor(Metric, levels = c("% Aligned", "M Aligned")))

# MultiQC-style table: light green for % Aligned, light purple for M Aligned
cols <- c("% Aligned" = "#dbe9dd", "M Aligned" = "#e7d4e0")

p_salmon_quant_stats <- ggplot(plot_data, aes(x = Value, y = Tissue)) +
  geom_col(aes(fill = Metric), width = 0.8) +
  geom_text(aes(label = LabelText), x = 0.5, hjust = 0, size = 3, color = "black") +
  facet_wrap(~Metric, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL, x = NULL) +
  ggtitle("Quantification stats") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 9, color = "black", hjust = 0),
        panel.grid.major.y = element_line(color = "gray90"),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 10, face = "bold", color = "#337ab7"),
        plot.title = element_text(size = 11),
        panel.spacing = unit(0.5, "lines"),
        plot.margin = margin(10, 10, 10, 10))
