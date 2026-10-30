# Fig 3 (absolute expression) and Fig S5 (Z-scores) of VD3-related genes, split by level
# Fig 3 uses fixed orders, as Fig 4: genes in pathway order (rows of tpm_vd3_genes.tsv) and tissues by organ
# system (tissue_order in 00_config.R). Fig S5 is clustered.
# Input:  vd3_genes/tpm_vd3_genes.tsv (from 00_1)
# Output: vd3_genes/vd3_genes_all_expr.{svg,png} (Fig 3), vd3_genes/vd3_genes_all_expr_z.{svg,png} (Fig S5)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

vd3 <- read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))
long <- to_long(vd3$log_tpm)
tau <- apply(vd3$log_tpm, 1, tau_index)
peak <- apply(vd3$tpm, 1, max)

groups <- expr_level_groups(vd3$log_tpm)
legend_breaks <- list(High = NA, Medium = c(-3, -1.5, 0, 1.5, 3), Low = NA)

panels <- imap(groups, function(genes, level) {
  z <- vd3$z[genes, ]
  annotation <- pattern_annotation(vd3$sd[genes], tau[genes], peak[genes])
  colors <- modifyList(pattern_colors, list(Pattern = pattern_colors$Pattern[levels(annotation$Pattern)]))
  extra <- if (anyNA(legend_breaks[[level]])) list() else list(legend_breaks = legend_breaks[[level]])
  ph <- do.call(z_heatmap, c(list(z, vd3$sd,
                                  paste0("Relative tissue-specific expression of VD3-related genes (", level, " expression)"),
                                  annotation = annotation, annotation_colors = colors),
                             extra))
  tissues <- tissue_label(intersect(tissue_order, str_to_lower(str_replace_all(colnames(z), " ", "_"))))
  stopifnot(setequal(tissues, colnames(z)))
  expr <- expr_tile_plot(long, intersect(pathway_order(rownames(vd3$log_tpm)), genes), tissues,
                         paste(level, "expression"))
  expr$layers[[2]]$aes_params$size <- 1.2  # crosses as in Fig 4
  # Font sizes and styles as in Fig 4 (05_5): 9 pt text, 7 pt legends, italic gene names, white background
  expr <- expr + labs(x = NULL, y = NULL) +  # no empty axis-title lines
    theme(text = element_text(size = 9),
          axis.text.y = element_text(face = "italic"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 9, hjust = 0),
          plot.title.position = "panel",
          plot.margin = margin(8, 5.5, 2, 30),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key.height = unit(if (level == "High") 0.22 else 0.3, "cm"),
          legend.key.width = unit(0.35, "cm")) +
    # tissue names only under the bottom panel, as in Fig 4 (all panels share the tissues)
    if (level == "Low") theme() else theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  list(z = ph$gtable, expr = expr)
})

# Panel heights proportional to the number of genes, so that the tiles are the same height in all three panels
p_vd3_all <- patchwork::wrap_plots(panels$High$expr, panels$Medium$expr, panels$Low$expr, ncol = 1,
                                   heights = lengths(groups[c("High", "Medium", "Low")])) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1), plot.tag.location = "plot",
        plot.tag.position = c(-0.03, 1))  # the tag shares the line of the panel title, in the left margin (as in Fig 4)
save_fig(p_vd3_all, file.path(vd3_dir, "vd3_genes_all_expr"), 180, 112)

p_vd3_z_all <- plot_grid(panels$High$z, panels$Medium$z, panels$Low$z,
                         nrow = 3, rel_heights = c(1, 1.4, 1.4), labels = c("a", "b", "c"), label_size = 11)
save_fig(p_vd3_z_all, file.path(vd3_dir, "vd3_genes_all_expr_z"), 180, 160)
