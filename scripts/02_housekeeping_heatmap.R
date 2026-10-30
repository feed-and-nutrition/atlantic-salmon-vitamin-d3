# Housekeeping gene expression and stability (exploratory, Figure S3 until v7 round 8, not in run_all)
# Input:  hk_genes/tpm_hk_genes.tsv (from 00_1)
# Output: hk_genes/hk_genes.{svg,png}, hk_genes/summary_hk_genes.tsv

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

hk <- read_gene_tpm(file.path(hk_dir, "tpm_hk_genes.tsv"))

# ef1a1 and eef1a1l3 are barely expressed, so they are left out of the Z-score panel
low_genes <- c("ef1a1", "eef1a1l3")
z_hk <- hk$z[!(rownames(hk$z) %in% low_genes), ]

ph <- z_heatmap(z_hk, hk$sd, "Relative tissue-specific expression of housekeeping genes")
ord <- ph_order(ph, z_hk)

p_expr <- expr_tile_plot(to_long(hk$log_tpm), c(ord$genes, low_genes), ord$tissues,
                         "Expression levels of housekeeping genes")

p_hk_genes <- plot_grid(p_expr, ph$gtable,
                        nrow = 2, rel_heights = c(1.1, 1), labels = c("a", "b"), label_size = 11)
save_fig(p_hk_genes, file.path(hk_dir, "hk_genes"), 180, 150)

expr_summary_hk <- to_long(hk$log_tpm) %>%
  group_by(gene) %>%
  summarise(zeros = sum(logTPM == 0),
            min_logtpm = round(min(logTPM), 2),
            mean_logtpm = round(mean(logTPM), 2),
            median_logtpm = round(median(logTPM), 2),
            max_logtpm = round(max(logTPM), 2),
            sd = round(sd(logTPM), 2)) %>%
  mutate(stability = ifelse(sd < 1, "Stable", "Variable"))
write_tsv(expr_summary_hk, file.path(hk_dir, "summary_hk_genes.tsv"))
