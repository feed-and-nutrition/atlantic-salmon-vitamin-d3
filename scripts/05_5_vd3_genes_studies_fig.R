# Fig 4 and Fig S5: VD3-related gene expression in independent RNA-seq studies
# Uses the replicate studies (replicate_studies in 00_config.R). Excluded samples (00_config.R) are dropped.
# All samples of a tissue are used, whatever their treatment, sex, or stage.
# Fig 4 pools the studies per tissue, with each study weighted equally so a large study does not dominate.
# Fig S5 keeps the studies apart (per-study medians and cyp27b1 in individual samples).
# Input:  vd3_genes/tpm_vd3_genes_studies.tsv (from 05_2), vd3_genes/tpm_vd3_genes.tsv (Fig 3 gene groups)
# Output: vd3_genes/studies/vd3_genes_tissues.{svg,png} (Fig 4)
#         vd3_genes/studies/vd3_genes_per_study.{svg,png} (Fig S5)
#         vd3_genes/studies/tissue_means.tsv (Fig 4 values per gene and tissue)
#         vd3_genes/studies/study_medians.tsv (n, median and maximum TPM per gene, study, and tissue)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_dir <- file.path(vd3_dir, "studies")
dir.create(out_dir, showWarnings = FALSE)

fig_studies <- replicate_studies

vd3 <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"), col_types = cols(gene_id = col_character())) %>%
  filter(study %in% fig_studies, !excluded)
stopifnot(setequal(unique(vd3$tissue), intersect(tissue_order, vd3$tissue)))
gene_levels <- pathway_order(unique(vd3$name))  # pathway order (No in salmon_vd3_orthologs_v1.tsv)

medians <- vd3 %>%
  group_by(name, study, tissue) %>%
  summarise(n = n(), median_tpm = median(tpm), max_tpm = max(tpm), .groups = "drop")
write_tsv(mutate(medians, across(c(median_tpm, max_tpm), ~ round(.x, 3))), file.path(out_dir, "study_medians.tsv"))

# a: median log2(TPM + 1) per study and tissue, white x where the median TPM <= 1 (as in Fig 3)
row_labels <- medians %>%
  distinct(study, tissue, n) %>%
  mutate(row = paste0(tissue_label(tissue), " (", n, ")"))
tile_df <- medians %>%
  left_join(row_labels, by = c("study", "tissue", "n")) %>%
  mutate(study = factor(study, levels = fig_studies),
         row = fct_rev(fct_reorder(row, match(tissue, tissue_order))),
         name = factor(name, levels = gene_levels),
         log_tpm = log2(median_tpm + 1))

p_tile <- ggplot(tile_df, aes(name, row, fill = log_tpm)) +
  geom_tile() +
  geom_point(data = filter(tile_df, median_tpm <= 1), colour = "white", shape = 4, size = 1.2) +
  facet_grid(study ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_viridis_c(name = "log2(median\nTPM + 1)") +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "Median expression of VD3-related genes per study and tissue") +
  theme_grey(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"),
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        strip.placement = "outside",
        panel.spacing.y = unit(0.6, "mm"),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.5, "cm"))

# b: every sample for the two cyp27b1 paralogs, so absence in kidney and head kidney is visible per sample
study_colors <- set_names(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#882255"), fig_studies)
study_shapes <- set_names(c(21, 22, 23, 24, 25, 21, 22), fig_studies)
point_df <- vd3 %>%
  filter(name %in% c("cyp27b1.a", "cyp27b1.b")) %>%
  mutate(tissue = factor(tissue_label(tissue), levels = tissue_label(intersect(tissue_order, vd3$tissue))),
         study = factor(study, levels = fig_studies),
         log_tpm = log2(tpm + 1))

p_points <- ggplot(point_df, aes(tissue, log_tpm, fill = study, shape = study)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  geom_point(position = position_jitter(width = 0.25, height = 0, seed = 1),
             size = 1.3, stroke = 0.2, colour = "grey15", alpha = 0.85) +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = study_colors, name = NULL) +
  scale_shape_manual(values = study_shapes, name = NULL) +
  labs(x = NULL, y = "log2(TPM + 1)", title = "cyp27b1 paralogs in individual samples") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "italic"),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.title.position = "plot",
        legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2)))

p_s6 <- plot_grid(p_tile, p_points, ncol = 1, rel_heights = c(1.55, 1), labels = c("a", "b"), label_size = 11)
save_fig(p_s6, file.path(out_dir, "vd3_genes_per_study"), 180, 230)

# Fig 4: tiles as in Fig 3, split into the Fig 3 expression groups.
# Colour = mean log2(TPM + 1), averaged over studies. White x = mean of 1 or lower (TPM <= 1), as in Fig 3.
# tissue_means.tsv also keeps the percentage of samples above 1 TPM.
groups <- expr_level_groups(read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))$log_tpm)
tissue_means <- vd3 %>%
  group_by(name, tissue, study) %>%
  summarise(n = n(), mean_log = mean(log2(tpm + 1)), pct_tpm1 = 100 * mean(tpm > 1), .groups = "drop") %>%
  group_by(name, tissue) %>%
  summarise(n_studies = n(), n = sum(n), mean_log = mean(mean_log), pct_tpm1 = mean(pct_tpm1), .groups = "drop")
write_tsv(mutate(tissue_means, across(c(mean_log, pct_tpm1), ~ round(.x, 2))), file.path(out_dir, "tissue_means.tsv"))

# Tissue labels give the number of samples and of studies, e.g. "Liver (34, 3)"
tissue_labels <- tissue_means %>%
  distinct(tissue, n, n_studies) %>%
  arrange(match(tissue, tissue_order)) %>%
  mutate(label = paste0(tissue_label(tissue), " (", n, ", ", n_studies, ")"))
stopifnot(!anyDuplicated(tissue_labels$tissue))

tile_panel <- function(level, bottom) {
  d <- tissue_means %>%
    filter(name %in% groups[[level]]) %>%
    mutate(tissue = factor(tissue, levels = tissue_labels$tissue, labels = tissue_labels$label),
           name = fct_rev(factor(name, levels = gene_levels)))
  ggplot(d, aes(tissue, name, fill = mean_log)) +
    geom_tile() +
    geom_point(data = filter(d, mean_log <= 1), colour = "white", shape = 4, size = 1.2) +
    scale_fill_viridis_c(name = "Mean\nlog2(TPM + 1)", limits = c(0, NA)) +
    scale_y_discrete(position = "right") +
    labs(x = NULL, y = NULL, title = paste(level, "expression")) +
    theme_grey(base_size = 9) +
    theme(axis.text.x = if (bottom) element_text(angle = 45, hjust = 1, vjust = 1) else element_blank(),
          axis.ticks.x = if (bottom) element_line() else element_blank(),
          axis.text.y = element_text(face = "italic"),
          plot.title = element_text(size = 9, hjust = 0),
          plot.margin = margin(5.5, 5.5, 5.5, 30),  # room for the a/b/c labels and the first rotated tissue label
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key.height = unit(if (level == "High") 0.22 else 0.35, "cm"))
}

p_fig4 <- plot_grid(tile_panel("High", FALSE), tile_panel("Medium", FALSE), tile_panel("Low", TRUE),
                    ncol = 1, rel_heights = c(length(groups$High) + 3.5, length(groups$Medium) + 1.5,
                                              length(groups$Low) + 6.5),
                    labels = c("a", "b", "c"), label_size = 11, align = "v", axis = "lr")
save_fig(p_fig4, file.path(out_dir, "vd3_genes_tissues"), 180, 125)
