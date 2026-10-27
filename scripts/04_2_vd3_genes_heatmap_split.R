# Fig 3 (absolute expression) and Fig S6 (Z-scores) of VD3-related genes, split by level
# Fig 3 uses fixed orders, as Fig 4: genes in pathway order (rows of tpm_vd3_genes.tsv) and tissues by organ
# system (tissue_order in 00_config.R). Fig S6 is clustered.
# Input:  vd3_genes/tpm_vd3_genes.tsv (from 00_1)
# Output: vd3_genes/vd3_genes_all_expr.{svg,png} (Fig 3), vd3_genes/vd3_genes_all_expr_z.{svg,png} (Fig S6)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

vd3 <- read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))
long <- to_long(vd3$log_tpm)
tau <- apply(vd3$log_tpm, 1, tau_index)

groups <- expr_level_groups(vd3$log_tpm)
legend_breaks <- list(High = NA, Medium = c(-3, -1.5, 0, 1.5, 3), Low = NA)

panels <- imap(groups, function(genes, level) {
  z <- vd3$z[genes, ]
  extra <- if (anyNA(legend_breaks[[level]])) list() else list(legend_breaks = legend_breaks[[level]])
  ph <- do.call(z_heatmap, c(list(z, vd3$sd,
                                  paste0("Relative tissue-specific expression of VD3-related genes (", level, " expression)"),
                                  annotation = pattern_annotation(vd3$sd[genes], tau[genes]),
                                  annotation_colors = pattern_colors),
                             extra))
  tissues <- tissue_label(intersect(tissue_order, str_to_lower(str_replace_all(colnames(z), " ", "_"))))
  stopifnot(setequal(tissues, colnames(z)))
  list(z = ph$gtable,
       expr = expr_tile_plot(long, intersect(pathway_order(rownames(vd3$log_tpm)), genes), tissues,
                             paste0("Expression levels of VD3-related genes (", level, " expression)"),
                             key_height = 0.25))
})

p_vd3_all <- plot_grid(panels$High$expr, panels$Medium$expr, panels$Low$expr,
                       nrow = 3, rel_heights = c(1, 1.4, 1.4), labels = c("a", "b", "c"), label_size = 11)
save_fig(p_vd3_all, file.path(vd3_dir, "vd3_genes_all_expr"), 180, 160)

p_vd3_z_all <- plot_grid(panels$High$z, panels$Medium$z, panels$Low$z,
                         nrow = 3, rel_heights = c(1, 1.4, 1.4), labels = c("a", "b", "c"), label_size = 11)
save_fig(p_vd3_z_all, file.path(vd3_dir, "vd3_genes_all_expr_z"), 180, 160)
