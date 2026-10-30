# Fig 4 and Fig S4: VD3-related gene expression in independent RNA-seq studies
# Uses the replicate studies (replicate_studies in 00_config.R). Excluded samples (00_config.R) are dropped.
# All samples of a tissue are used, whatever their treatment, sex, or stage.
# Fig 4 pools the studies per tissue, with each study weighted equally so a large study does not dominate.
# Fig S4 keeps the studies apart (per-study medians and cyp27b1 in individual samples).
# Input:  vd3_genes/tpm_vd3_genes_studies.tsv (from 05_2), vd3_genes/tpm_vd3_genes.tsv (Fig 3 gene groups)
# Output: vd3_genes/studies/vd3_genes_tissues.{svg,png} (Fig 4)
#         vd3_genes/studies/vd3_genes_per_study.{svg,png} (Fig S4)
#         vd3_genes/studies/tissue_means.tsv (Fig 4 values per gene and tissue)
#         vd3_genes/studies/atlas_comparison.tsv, atlas_peaks.tsv (Fig 4 tiles that differ from Fig 3, atlas peaks)
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
  facet_grid(study ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = as_labeller(\(x) paste0(study_names[x], "\n", x))) +  # name over project number
  scale_fill_viridis_c(name = "log2(median\nTPM + 1)") +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "Median expression of VD3-related genes per study and tissue") +
  theme_grey(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"),
        strip.text.y.left = element_text(angle = 0, hjust = 0, size = 7.5, lineheight = 0.85, margin = margin(1, 4, 1, 4)),
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
  scale_fill_manual(values = study_colors, labels = study_names[fig_studies], name = NULL) +
  scale_shape_manual(values = study_shapes, labels = study_names[fig_studies], name = NULL) +
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

# Comparison with the atlas (Fig 3), in the 14 tissues it shares with the studies. A tile differs when the studies and
# the atlas lie on opposite sides of 1 TPM (log2(TPM + 1) = 1) and differ at least twofold (1 log2 unit), so that
# values near the cutoff do not count. The peak column shows whether the atlas peak tissue of each gene is also the
# highest of these tissues in the studies. Tissues without an atlas sample (heart excluded) are marked below the heatmaps.
atlas <- read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))$log_tpm %>%
  as_tibble(rownames = "name") %>%
  pivot_longer(-name, names_to = "label", values_to = "atlas_log")
compare <- tissue_means %>%
  mutate(label = tissue_label(tissue)) %>%
  left_join(atlas, by = c("name", "label")) %>%
  mutate(differs = !is.na(atlas_log) & (mean_log > 1) != (atlas_log > 1) & abs(mean_log - atlas_log) >= 1)
peaks <- compare %>%
  filter(!is.na(atlas_log)) %>%
  group_by(name) %>%
  summarise(atlas_peak = label[which.max(atlas_log)], studies_peak = label[which.max(mean_log)], .groups = "drop") %>%
  mutate(same_peak = atlas_peak == studies_peak)
write_tsv(mutate(select(compare, name, tissue, mean_log, atlas_log, differs), across(where(is.double), ~ round(.x, 2))),
          file.path(out_dir, "atlas_comparison.tsv"))
write_tsv(peaks, file.path(out_dir, "atlas_peaks.tsv"))
message(sum(compare$differs), " of ", sum(!is.na(compare$atlas_log)), " shared tiles differ from the atlas, ",
        sum(peaks$same_peak), " of ", nrow(peaks), " atlas peaks reproduced")
new_tissues <- tissue_labels$tissue[!tissue_label(tissue_labels$tissue) %in% atlas$label]
x_levels <- c(tissue_labels$label, "Atlas peak")
diff_colour <- "#FF4500"

tile_panel <- function(level, bottom) {
  d <- compare %>%
    filter(name %in% groups[[level]]) %>%
    mutate(tissue = factor(tissue, levels = tissue_labels$tissue, labels = tissue_labels$label),
           name = fct_rev(factor(name, levels = gene_levels)))
  y_levels <- c(if (bottom) "new", levels(droplevels(d$name)))  # the bottom panel ends with a row marking the new tissues
  d <- mutate(d, name = factor(name, levels = y_levels))
  pk <- peaks %>%
    filter(name %in% groups[[level]]) %>%
    mutate(x = factor("Atlas peak", levels = x_levels), name = factor(name, levels = y_levels),
           peak = if_else(same_peak, "also highest in the studies", "not highest in the studies"))
  new_row <- tibble(x = factor(tissue_labels$label[match(new_tissues, tissue_labels$tissue)], levels = x_levels),
                    name = factor("new", levels = y_levels))
  p <- ggplot(d, aes(factor(tissue, levels = x_levels), name)) +
    geom_tile(aes(fill = mean_log)) +
    geom_point(data = filter(d, mean_log <= 1), colour = "white", shape = 4, size = 1.2) +
    geom_tile(data = filter(d, differs), aes(colour = "differs from atlas"), fill = NA,
              linewidth = 0.7, width = 0.9, height = 0.9) +
    geom_point(data = pk, aes(x, name, shape = peak), size = 2, stroke = 0.6, colour = "grey15", fill = "grey15",
               inherit.aes = FALSE) +
    scale_fill_viridis_c(name = "Mean\nlog2(TPM + 1)", limits = c(0, NA)) +
    scale_colour_manual(values = diff_colour, name = NULL,
                        guide = if (bottom) guide_legend(position = "bottom") else "none") +
    scale_shape_manual(values = c(`also highest in the studies` = 16, `not highest in the studies` = 1),
                       limits = c("also highest in the studies", "not highest in the studies"),
                       guide = "none") +  # filled and open circles are explained in the figure legend
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(position = "right", drop = FALSE, breaks = setdiff(y_levels, "new")) +
    labs(x = NULL, y = NULL, title = paste(level, "expression")) +
    theme_grey(base_size = 9) +
    theme(axis.text.x = if (bottom) element_text(angle = 45, hjust = 1, vjust = 1) else element_blank(),
          axis.ticks.x = if (bottom) element_line() else element_blank(),
          axis.text.y = element_text(face = "italic"),
          panel.background = element_blank(),  # white behind the Atlas peak column
          panel.grid = element_blank(),
          plot.title = element_text(size = 9, hjust = 0),
          plot.margin = margin(8, 5.5, 2, 30),  # room for the a/b/c labels and the first rotated tissue label
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key = element_blank(),
          legend.key.height = unit(if (level == "High") 0.22 else 0.3, "cm"),
          legend.key.width = unit(0.35, "cm"),
          legend.spacing.y = unit(0.5, "mm"),
          legend.margin = margin(0, 0, 0, 0))
  # tissues without an atlas sample, labelled "new" below the bottom heatmap (explained in the figure legend)
  if (bottom) p <- p +
    geom_tile(data = new_row, aes(x, name), fill = "grey25", height = 0.7, inherit.aes = FALSE) +
    geom_text(data = new_row, aes(x, name), label = "new", colour = "white", size = 2.2, inherit.aes = FALSE)
  p
}

# Panel heights proportional to the number of rows (panel c has the extra "new" row), so all tiles are the same height
p_fig4 <- patchwork::wrap_plots(tile_panel("High", FALSE), tile_panel("Medium", FALSE), tile_panel("Low", TRUE),
                                ncol = 1, heights = c(length(groups$High), length(groups$Medium), length(groups$Low) + 1)) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1), plot.tag.location = "plot",
        plot.tag.position = c(-0.03, 1))  # the tag shares the line of the panel title, in the left margin (as in Fig 3)
save_fig(p_fig4, file.path(out_dir, "vd3_genes_tissues"), 180, 125)
