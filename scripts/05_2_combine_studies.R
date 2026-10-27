# Combine the per-study exports from the server (05_1) into one sample table and
# long tables of VD3-related and housekeeping gene TPM across all studies.
# Input:  ncbi_samples/sample_info.txt (study order and sheet names)
#         salmon_quant/<Project_no>/{samples.tsv,tpm_gene.tsv.gz,counts_gene.tsv.gz}
# Output: salmon_quant/samples_all.tsv (one row per sample, harmonised tissue, QC, exclusion flag)
#         vd3_genes/tpm_vd3_genes_studies.tsv, hk_genes/tpm_hk_genes_studies.tsv (long)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

studies <- read_csv(file.path(proj_dir, "ncbi_samples", "sample_info.txt"), show_col_types = FALSE) %>%
  rename_with(str_trim) %>%
  mutate(across(everything(), str_trim),
         dataset = str_remove(sample_info, "\\.txt$"))

# Sheet tissue labels -> lowercase names used in expected_tissues. Gonads are split by sex.
# PRJNA647285 "Kidney" is head kidney: these samples match PRJEB47409 head kidney best (salmon_quant/tissue_identity.tsv,
# from 05_4). The study (Gervais et al. 2021, doi:10.1186/s12864-021-07671-6) used uninfected 0 dpi fish.
harmonise_tissue <- function(tissue, sex, study) {
  t <- str_replace_all(str_to_lower(tissue), " ", "_")
  case_when(t == "gonad" & sex == "M" ~ "testis",
            t == "gonad" & sex == "F" ~ "ovary",
            t == "kidney" & study == "PRJNA647285" ~ "head_kidney",
            TRUE ~ t)
}

read_study <- function(project) {
  dir <- file.path(proj_dir, "salmon_quant", project)
  samples <- read_tsv(file.path(dir, "samples.tsv"), col_types = cols(.default = col_character()))
  tpm <- read_tsv(file.path(dir, "tpm_gene.tsv.gz"), col_types = cols(gene_id = col_character(), .default = col_double()))
  counts <- read_tsv(file.path(dir, "counts_gene.tsv.gz"), col_types = cols(gene_id = col_character(), .default = col_double()))
  qc <- tibble(txi_column = colnames(tpm)[-1],
               m_counts = round(colSums(counts[, -1]) / 1e6, 2),
               n_genes_tpm1 = colSums(tpm[, -1] > 1))
  list(samples = left_join(samples, qc, by = "txi_column"), tpm = tpm)
}

data <- set_names(map(studies$Project_no, read_study), studies$Project_no)

samples_all <- map_dfr(data, "samples") %>%
  left_join(select(studies, study = Project_no, dataset), by = "study") %>%
  mutate(tissue_raw = tissue,
         tissue = harmonise_tissue(tissue, sex, study),
         sex = na_if(na_if(sex, ""), "NA"),
         stage = na_if(stage, ""),
         id = paste(study, txi_column, sep = ":"),
         excluded = id %in% names(excluded_samples)) %>%
  select(id, study, dataset, sample = txi_column, tissue, tissue_raw, sex, stage, treatment, fish_id,
         m_counts, n_genes_tpm1, excluded)
write_tsv(samples_all, file.path(proj_dir, "salmon_quant", "samples_all.tsv"))

# Long TPM table for a gene list (name, gene_id), one row per gene and sample
export_genes_long <- function(genes, out_file) {
  out <- imap_dfr(data, function(d, project) {
    d$tpm %>%
      filter(gene_id %in% genes$gene_id) %>%
      pivot_longer(-gene_id, names_to = "sample", values_to = "tpm") %>%
      mutate(id = paste(project, sample, sep = ":"))
  }) %>%
    inner_join(select(genes, name, gene_id), by = "gene_id") %>%
    inner_join(select(samples_all, -sample), by = "id") %>%
    select(name, gene_id, study, dataset, id, tissue, sex, stage, treatment, excluded, tpm) %>%
    arrange(factor(name, levels = genes$name), study, tissue)
  write_tsv(out, out_file)
}

vd3_genes <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(gene_id = col_character())) %>%
  filter(!(name %in% names(excluded_genes))) %>%
  arrange(No)
export_genes_long(vd3_genes, file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"))

hk_genes <- read_tsv(file.path(hk_dir, "hk_genes.tsv"), col_types = cols(gene_id = col_character())) %>%
  filter(!(gene_id %in% c("100196887", "106576970"))) %>%
  arrange(no)
export_genes_long(hk_genes, file.path(hk_dir, "tpm_hk_genes_studies.tsv"))

samples_all %>% count(study, tissue) %>% pivot_wider(names_from = study, values_from = n, values_fill = 0) %>%
  print(n = Inf, width = Inf)
