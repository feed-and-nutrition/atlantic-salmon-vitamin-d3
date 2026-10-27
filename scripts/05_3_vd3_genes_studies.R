# Expression of VD3-related genes across all RNA-seq studies (replication of the PRJNA72713 atlas)
# Input:  vd3_genes/tpm_vd3_genes_studies.tsv (from 05_2)
# Output: vd3_genes/studies/vd3_genes_studies.{svg,png}   log2(TPM + 1) per sample, by tissue and study
#         vd3_genes/studies/tissue_summary.tsv            median TPM and fraction of samples with TPM > 1
#         vd3_genes/studies/concordance_atlas.tsv         Spearman rho of tissue medians vs the atlas

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_dir <- file.path(vd3_dir, "studies")
dir.create(out_dir, showWarnings = FALSE)
atlas <- "PRJNA72713"

vd3 <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"), col_types = cols(gene_id = col_character())) %>%
  filter(!excluded) %>%
  mutate(log_tpm = log2(tpm + 1))
gene_levels <- unique(vd3$name)
study_levels <- c(atlas, setdiff(unique(vd3$study), atlas))
tissue_levels <- c(expected_tissues, setdiff(sort(unique(vd3$tissue)), expected_tissues))

tissue_summary <- vd3 %>%
  group_by(name, study, tissue) %>%
  summarise(n = n(), median_tpm = round(median(tpm), 3), frac_tpm1 = round(mean(tpm > 1), 2), .groups = "drop")
write_tsv(tissue_summary, file.path(out_dir, "tissue_summary.tsv"))

# Does each study rank the shared tissues like the atlas? Only studies with >= 4 shared tissues.
atlas_profile <- tissue_summary %>% filter(study == atlas) %>% select(name, tissue, atlas_tpm = median_tpm)
concordance <- tissue_summary %>%
  filter(study != atlas) %>%
  inner_join(atlas_profile, by = c("name", "tissue")) %>%
  group_by(name, study) %>%
  filter(n() >= 4) %>%
  summarise(n_tissues = n(),
            rho = round(cor(median_tpm, atlas_tpm, method = "spearman"), 2),
            top_atlas = tissue[which.max(atlas_tpm)],
            top_study = tissue[which.max(median_tpm)],
            .groups = "drop") %>%
  arrange(factor(name, levels = gene_levels), study)
write_tsv(concordance, file.path(out_dir, "concordance_atlas.tsv"))

study_colors <- set_names(c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                            "#CC79A7", "grey50")[seq_along(study_levels)], study_levels)

plot_df <- vd3 %>%
  mutate(name = factor(name, levels = gene_levels),
         study = factor(study, levels = study_levels),
         tissue = factor(tissue_label(tissue), levels = tissue_label(tissue_levels)))

p <- ggplot(plot_df, aes(tissue, log_tpm, colour = study)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  geom_point(data = filter(plot_df, study != atlas), position = position_jitter(width = 0.25, height = 0, seed = 1),
             size = 0.5, alpha = 0.7) +
  geom_point(data = filter(plot_df, study == atlas), shape = 23, size = 1.4, fill = "white", stroke = 0.6) +
  facet_wrap(~name, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = study_colors, name = NULL) +
  labs(x = NULL, y = "log2(TPM + 1)") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
        legend.position = "top",
        strip.text = element_text(face = "italic")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2, alpha = 1)))
save_fig(p, file.path(out_dir, "vd3_genes_studies"), 180, 260)
