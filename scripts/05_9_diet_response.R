# Dietary vitamin D3 and the VD3-related genes in PRJNA1160017 skin (tracker item 38): juvenile salmon of two strains
# (Native Irish, Norwegian commercial) fed 866 or 1779 ng/g vitamin D3, three fish per diet and strain. Not in the
# figures, and not one of the replicate studies of Fig 4. The stomach samples are left out, because no study using them
# has been published and they are of two kinds (gastric mucosa in fish 1, muscular wall without pepsinogen in fish 2 and 3).
# Skin samples whose median Spearman correlation with the other skin samples is below 0.8 are flagged as atypical (the
# three Norwegian fish fed 866 ng/g, 0.71 to 0.76, the others 0.88 or more), and skin is fitted with and without them.
# Model per gene: log2(TPM + 1) ~ strain + log2(dose), genes with median TPM >= 1, BH within each sample set.
# Sheet Dietary vitamin D3 of Table S9 (05_7).
# Input:  salmon_quant/samples_all.tsv, vd3_genes/tpm_vd3_genes_studies.tsv (05_2),
#         salmon_quant/PRJNA1160017/{samples,tpm_gene}.tsv(.gz) (05_1 on the server)
# Output: vd3_genes/studies/diet_PRJNA1160017/{samples,medians,fits}.tsv

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

study <- "PRJNA1160017"
out_dir <- file.path(vd3_dir, "studies", "diet_PRJNA1160017")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tpm <- read_tsv(file.path(proj_dir, "salmon_quant", study, "tpm_gene.tsv.gz"),
                col_types = cols(gene_id = col_character(), .default = col_double()))
m <- as.matrix(tpm[, -1])
rownames(m) <- tpm$gene_id
sheet <- read_tsv(file.path(proj_dir, "salmon_quant", study, "samples.tsv"), col_types = cols(.default = col_character()))
samples <- read_tsv(file.path(proj_dir, "salmon_quant", "samples_all.tsv"), show_col_types = FALSE) %>%
  filter(study == !!study, tissue == "skin") %>%
  left_join(select(sheet, sample = txi_column, strain), by = "sample") %>%
  mutate(dose = as.numeric(str_remove(treatment, "^vd3_")))

# Median Spearman correlation of each sample with the other skin samples (genes with TPM > 1 in any skin sample)
m <- m[, samples$sample]
log_m <- log2(m[apply(m, 1, max) > 1, ] + 1)
r <- cor(log_m, method = "spearman")
diag(r) <- NA
samples <- samples %>%
  mutate(r_within = apply(r, 1, median, na.rm = TRUE)[sample],
         group = if_else(r_within < 0.8, "skin atypical", "skin"))
write_tsv(mutate(select(samples, sample, strain, fish_id, dose, m_counts, r_within, group), r_within = round(r_within, 3)),
          file.path(out_dir, "samples.tsv"))
print(count(samples, group, dose, strain), n = Inf)

d <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"), show_col_types = FALSE) %>%
  filter(study == !!study) %>%
  mutate(sample = str_remove(id, "^[^:]+:")) %>%
  inner_join(select(samples, sample, strain, dose, group), by = "sample")

# Median TPM per diet and strain (three fish each), with the atypical samples
medians <- d %>%
  group_by(name, dose, strain) %>%
  summarise(tpm = median(tpm), atypical = sum(group == "skin atypical"), .groups = "drop")
write_tsv(mutate(medians, across(where(is.numeric), ~ round(.x, 2))), file.path(out_dir, "medians.tsv"))

fit_group <- function(x) {
  x <- filter(x, median(tpm) >= 1)
  if (!nrow(x)) return(tibble())
  f <- lm(log2(tpm + 1) ~ strain + log2(dose), data = x)
  co <- summary(f)$coefficients
  tibble(n = nrow(x), median_tpm = median(x$tpm), fold_per_doubling = 2^co["log2(dose)", 1],
         p = co["log2(dose)", 4], p_strain = co[2, 4])
}
sets <- list(skin = "skin", `skin incl. atypical` = c("skin", "skin atypical"))
fits <- imap_dfr(sets, function(groups, label) {
  d %>% filter(group %in% groups) %>% group_by(name) %>% group_modify(~ fit_group(.x)) %>% ungroup() %>%
    mutate(samples = label, .after = name)
}) %>%
  group_by(samples) %>%
  mutate(fdr = p.adjust(p, "BH")) %>%
  ungroup()
write_tsv(mutate(fits, across(where(is.double), ~ signif(.x, 3))), file.path(out_dir, "fits.tsv"))
print(arrange(fits, p), n = Inf, width = Inf)
