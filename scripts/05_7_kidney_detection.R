# Table S11: detection power for cyp27b1 in the kidney and head kidney samples (atlas and replicate studies).
# For each sample: library size (estimated reads), estimated reads and TPM of both cyp27b1 paralogs, and the TPM
# that one read of each paralog corresponds to. Where a sample has no read of the gene, that value is estimated
# from the samples of the same study that have reads, scaled by library size (TPM per read is about inversely
# proportional to the number of reads).
# Needs python3 for the Excel file.
# Input:  salmon_quant/samples_all.tsv (05_2), salmon_quant/<study>/{tpm_gene,counts_gene}.tsv.gz,
#         vd3_genes/tpm_vd3_genes.tsv (gene IDs)
# Output: vd3_genes/studies/kidney_cyp27b1.tsv, manuscript/v6/Table S11 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_file <- file.path(vd3_dir, "studies", "kidney_cyp27b1.tsv")
table_file <- file.path(fig_dir, "Table S11 Detection of cyp27b1 in kidney and head kidney samples.xlsx")

genes <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes.tsv"), col_types = cols(Gene_id = col_character())) %>%
  filter(Name %in% c("cyp27b1.a", "cyp27b1.b")) %>%
  select(name = Name, gene_id = Gene_id)
samples <- read_tsv(file.path(proj_dir, "salmon_quant", "samples_all.tsv"), show_col_types = FALSE) %>%
  filter(!excluded)

read_long <- function(study, file, value) {
  read_tsv(file.path(proj_dir, "salmon_quant", study, file),
           col_types = cols(gene_id = col_character(), .default = col_double())) %>%
    filter(gene_id %in% genes$gene_id) %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = value) %>%
    mutate(study = study)
}
studies <- unique(filter(samples, tissue %in% c("kidney", "head_kidney"))$study)
long <- map_dfr(studies, function(s) {
  inner_join(read_long(s, "tpm_gene.tsv.gz", "tpm"), read_long(s, "counts_gene.tsv.gz", "reads"),
             by = c("gene_id", "sample", "study"))
}) %>%
  inner_join(genes, by = "gene_id") %>%
  inner_join(select(samples, study, sample, tissue, tissue_raw, m_counts), by = c("study", "sample"))

# TPM of one read: observed where the gene has reads, else the study median of (TPM / reads) x (library size ratio)
long <- long %>%
  group_by(study, name) %>%
  mutate(per_read_obs = ifelse(reads >= 0.5, tpm / reads, NA),
         per_read = ifelse(!is.na(per_read_obs), per_read_obs,
                           median(per_read_obs * m_counts, na.rm = TRUE) / m_counts)) %>%
  ungroup()

kidney <- long %>%
  filter(tissue %in% c("kidney", "head_kidney")) %>%
  select(study, sample, tissue, tissue_raw, m_counts, name, reads, tpm, per_read) %>%
  pivot_wider(names_from = name, values_from = c(reads, tpm, per_read), names_glue = "{name}_{.value}") %>%
  mutate(study = factor(study, levels = unique(samples$study)),
         tissue = factor(tissue, levels = c("kidney", "head_kidney"))) %>%
  arrange(tissue, study, sample)
write_tsv(kidney, out_file)

fmt <- function(x, d) formatC(x, format = "f", digits = d)
table <- kidney %>%
  transmute(Study = as.character(study), Sample = sample,
            Tissue = ifelse(tissue_raw == tissue_label(as.character(tissue)) |
                              str_to_lower(tissue_raw) == str_replace(as.character(tissue), "_", " "),
                            tissue_label(as.character(tissue)),
                            paste0(tissue_label(as.character(tissue)), " (labelled ", tissue_raw, ")")),
            `Library size (M reads)` = m_counts,
            `cyp27b1.a reads` = fmt(cyp27b1.a_reads, 1), `cyp27b1.a TPM` = fmt(cyp27b1.a_tpm, 2),
            `cyp27b1.a TPM per read` = fmt(cyp27b1.a_per_read, 3),
            `cyp27b1.b reads` = fmt(cyp27b1.b_reads, 1), `cyp27b1.b TPM` = fmt(cyp27b1.b_tpm, 2),
            `cyp27b1.b TPM per read` = fmt(cyp27b1.b_per_read, 3))
tmp <- tempfile(fileext = ".tsv")
write_tsv(table, tmp, na = "")
note <- paste("Reads are Salmon estimated counts. TPM per read is the TPM that one read of the gene gives in the sample.",
              "Where the sample has no read, it is estimated from the samples of the same study with reads, scaled by",
              "library size. The PRJNA647285 samples labelled kidney were analysed as head kidney (Methods).")
cat("\n", note, "\n", file = tmp, append = TRUE, sep = "")
system2("python3", c(file.path(proj_dir, "scripts", "tsv_to_xlsx.py"), shQuote(tmp), shQuote(table_file)))
message("Wrote ", out_file, " and ", table_file)
