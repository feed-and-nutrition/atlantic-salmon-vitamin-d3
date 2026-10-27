# Fig S1: mapping statistics (a, from 03_1) and PCA of all genes (b)
# Input:  salmon_quant/PRJNA72713/tpm_atlantic_salmon_15_tissues.tsv.gz (from 00_1)
# Output: salmon_quant/PRJNA72713/salmon_stats.{svg,png} (Fig S1)

source(if (file.exists("03_1_salmon_quant_stats.R")) "03_1_salmon_quant_stats.R" else "scripts/03_1_salmon_quant_stats.R")
library(ggrepel)

if (!file.exists(tpm_all_file)) {
  stop("Missing ", tpm_all_file, ". Run 00_1_export_tpm.R with the 15-tissue txi_gene.rds first.", call. = FALSE)
}
tpm <- read_tsv(tpm_all_file, col_types = cols(ncbi_gene_id = col_character()))
check_tissues(colnames(tpm)[-1])

log_tpm <- log2(as.matrix(tpm[, -1]) + 1)
colnames(log_tpm) <- tissue_label(colnames(log_tpm))

# Samples as rows, genes as columns. Drop genes with zero variance across tissues.
pca_input <- t(log_tpm)
pca_input <- pca_input[, apply(pca_input, 2, var) > 0]
pca_res <- prcomp(pca_input, center = TRUE, scale. = FALSE)

pca_df <- as.data.frame(pca_res$x) %>% mutate(Sample = rownames(.))
var_explained <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = Inf, seed = 1) +
  xlab(paste0("PC1: ", var_explained[1], "% variance")) +
  ylab(paste0("PC2: ", var_explained[2], "% variance")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11)) +
  ggtitle(paste0("PCA of ", nrow(pca_df), " Tissues (log2 TPM)"))

p_salmon_stats <- plot_grid(p_salmon_quant_stats, p_pca, ncol = 2, rel_widths = c(1, 1),
                            labels = c("a", "b"), label_size = 11)
save_fig(p_salmon_stats, file.path(quant_dir, "salmon_stats"), 180, 90)
