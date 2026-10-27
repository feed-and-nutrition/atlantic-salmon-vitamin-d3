# Sample identity checks behind excluded_samples (00_config.R) and Fig S2
# Fig S2: the PRJNA72713 heart sample is mostly liver tissue, so it is excluded from the expression analyses.
# Liver markers: TPM > 100 in PRJEB47409 liver and < 1 in the hearts of PRJNA543940 and PRJNA647285 (medians).
# Heart markers: TPM > 100 in both reference hearts and < 1 in PRJEB47409 liver.
# Index of a sample = median over the markers of (sample TPM / reference TPM).
# Input:  salmon_quant/samples_all.tsv, salmon_quant/<study>/tpm_gene.tsv.gz (from 05_2 and the server)
#         Reference medians use only samples that are not excluded.
# Output: salmon_quant/PRJNA72713/atlas_heart_liver.{svg,png} (Fig S2)
#         salmon_quant/PRJNA72713/marker_index.tsv (liver and heart index of every sample)
#         salmon_quant/tissue_identity.tsv (each sample's closest study-tissue centroid from the other studies)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")
library(ggrepel)
set_colnames <- function(m, nm) { colnames(m) <- nm; m }

atlas <- "PRJNA72713"
samples <- read_tsv(file.path(proj_dir, "salmon_quant", "samples_all.tsv"), show_col_types = FALSE)
read_tpm <- function(study) {
  d <- read_tsv(file.path(proj_dir, "salmon_quant", study, "tpm_gene.tsv.gz"),
                col_types = cols(gene_id = col_character(), .default = col_double()))
  m <- as.matrix(d[, -1])
  rownames(m) <- d$gene_id
  m
}
tpm <- set_names(map(unique(samples$study), read_tpm), unique(samples$study))
genes <- rownames(tpm[[atlas]])
tissue_median <- function(study, tissue) {
  cols <- samples$sample[samples$study == study & samples$tissue == tissue & !samples$excluded]
  apply(tpm[[study]][genes, cols, drop = FALSE], 1, median)
}

liver_ref <- tissue_median("PRJEB47409", "liver")
heart_ref <- tissue_median("PRJNA543940", "heart")
heart_ref2 <- tissue_median("PRJNA647285", "heart")
liver_markers <- genes[liver_ref > 100 & heart_ref < 1 & heart_ref2 < 1]
heart_markers <- genes[heart_ref > 100 & heart_ref2 > 100 & liver_ref < 1]
message(length(liver_markers), " liver markers, ", length(heart_markers), " heart markers")

marker_index <- function(m, markers, ref) apply(m[markers, , drop = FALSE] / ref[markers], 2, median)
index <- imap_dfr(tpm, function(m, study) {
  tibble(study = study, sample = colnames(m),
         liver_index = marker_index(m, liver_markers, liver_ref),
         heart_index = marker_index(m, heart_markers, heart_ref))
}) %>%
  left_join(select(samples, study, sample, tissue), by = c("study", "sample"))
write_tsv(index, file.path(quant_dir, "marker_index.tsv"))

summary_atlas <- index %>% filter(study == atlas) %>% select(tissue, liver_index, heart_index)
other_max <- index %>% filter(tissue != "liver", !(study == atlas & tissue == "heart")) %>% pull(liver_index) %>% max()
message("Atlas heart: liver index ", round(summary_atlas$liver_index[summary_atlas$tissue == "heart"], 3),
        ", heart index ", round(summary_atlas$heart_index[summary_atlas$tissue == "heart"], 3),
        ". Highest liver index of any other non-liver sample: ", signif(other_max, 2))

# Spearman correlation of each sample with the mean log2(TPM + 1) of every tissue in the other studies
# (genes with mean log2(TPM + 1) > 1). Labels differ between studies (kidney vs head kidney, gut vs
# distal intestine), so read "closest" together with r_same, the best match among the same tissue label.
log_all <- do.call(cbind, imap(tpm, function(m, study) set_colnames(log2(m + 1), paste(study, colnames(m), sep = ":"))))
log_all <- log_all[rowMeans(log_all) > 1, samples$id]
centroid_key <- paste(samples$study, samples$tissue, sep = ":")
ok <- !samples$excluded
centroids <- sapply(split(samples$id[ok], centroid_key[ok]), function(ids) rowMeans(log_all[, ids, drop = FALSE]))
identity <- map_dfr(seq_len(nrow(samples)), function(i) {
  cc <- centroids[, !startsWith(colnames(centroids), paste0(samples$study[i], ":")), drop = FALSE]
  r <- cor(log_all[, samples$id[i]], cc, method = "spearman")[1, ]
  same <- r[str_remove(names(r), "^[^:]+:") == samples$tissue[i]]
  tibble(id = samples$id[i], tissue = samples$tissue[i], excluded = samples$excluded[i],
         closest = names(which.max(r)), r_closest = round(max(r), 3),
         r_same = if (length(same)) round(max(same), 3) else NA)
})
write_tsv(identity, file.path(proj_dir, "salmon_quant", "tissue_identity.tsv"))

# Log axes need a floor for samples without marker expression
floor_value <- 1e-4
plot_df <- index %>%
  mutate(group = case_when(tissue == "liver" ~ "Liver", tissue == "heart" ~ "Heart", TRUE ~ "Other tissues"),
         group = factor(group, levels = c("Liver", "Heart", "Other tissues")),
         source = if_else(study == atlas, "PRJNA72713", "Other studies"),
         across(c(liver_index, heart_index), ~ pmax(.x, floor_value)),
         label = if_else(study == atlas & tissue %in% c("heart", "liver"), paste(atlas, tissue_label(tissue)), NA))
group_colors <- c(Liver = "#D55E00", Heart = "#0072B2", `Other tissues` = "grey65")
log_breaks <- 10^(-4:1)
log_labels <- c("<0.0001", "0.001", "0.01", "0.1", "1", "10")

p_scatter <- ggplot(plot_df, aes(liver_index, heart_index, colour = group, shape = source)) +
  geom_point(data = filter(plot_df, source != atlas), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(plot_df, source == atlas), size = 2.2, stroke = 0.8) +
  geom_text_repel(aes(label = label), size = 2.8, colour = "black", min.segment.length = 0,
                  box.padding = 0.5, na.rm = TRUE, show.legend = FALSE, seed = 1) +
  scale_x_log10(breaks = log_breaks, labels = log_labels) +
  scale_y_log10(breaks = log_breaks, labels = log_labels) +
  scale_colour_manual(values = group_colors, name = NULL) +
  scale_shape_manual(values = c(`Other studies` = 16, PRJNA72713 = 5), name = NULL) +
  labs(x = "Liver marker index", y = "Heart marker index",
       title = paste0("Liver and heart markers in all samples (", nrow(plot_df), ")")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11), legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0, "mm"))

atlas_df <- plot_df %>%
  filter(study == atlas) %>%
  mutate(tissue = fct_reorder(tissue_label(tissue), liver_index))
p_atlas <- ggplot(atlas_df, aes(liver_index, tissue, colour = group)) +
  geom_vline(xintercept = other_max, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2) +
  scale_x_log10(breaks = log_breaks, labels = log_labels) +
  scale_colour_manual(values = group_colors, guide = "none") +
  labs(x = "Liver marker index", y = NULL, title = "PRJNA72713 tissues") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11))

p_fig <- plot_grid(p_scatter, p_atlas, ncol = 2, rel_widths = c(1.25, 1), labels = c("a", "b"), label_size = 11)
save_fig(p_fig, file.path(quant_dir, "atlas_heart_liver"), 180, 100)
