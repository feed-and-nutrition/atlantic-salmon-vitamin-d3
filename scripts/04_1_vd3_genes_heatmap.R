# Tables for VD3-related gene expression (basis of Table 2)
# Input:  vd3_genes/tpm_vd3_genes.tsv (from 00_1)
#         vd3_genes/tpm_vd3_genes_studies.tsv (from 05_2, for the replicated PRJNA359635 tissues)
# Output: vd3_genes/summary_vd3_genes.tsv (level, peak TPM, SD, tau, Table 2 pattern per gene)
#         vd3_genes/top3_tissues_vd3_genes.tsv (highest tissue, and the next two if their TPM is above 1)
#         vd3_genes/pattern_sensitivity.tsv (Table 2 pattern under other cutoffs and in PRJNA359635)
#         manuscript/v7/Table S8 ... .xlsx (sheet TPM: TPM per gene and atlas tissue with SD and tau, via
#         vd3_genes/table_s8.tsv; sheet Sensitivity: pattern_sensitivity.tsv)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

vd3 <- read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))
long <- to_long(vd3$log_tpm)

expr_summary_vd3 <- long %>%
  group_by(gene) %>%
  summarise(zeros = sum(logTPM == 0),
            min_logtpm = round(min(logTPM), 2),
            mean_logtpm = round(mean(logTPM), 2),
            median_logtpm = round(median(logTPM), 2),
            max_logtpm = round(max(logTPM), 2),
            sd = round(sd(logTPM), 2),
            tau = round(tau_index(logTPM), 2),
            peak_tpm = max(2^logTPM - 1)) %>%
  mutate(stability = ifelse(sd < 1, "Stable", "Variable"),
         pattern = expr_pattern(sd, tau, peak_tpm),
         peak_tpm = round(peak_tpm, 2))
write_tsv(expr_summary_vd3, file.path(vd3_dir, "summary_vd3_genes.tsv"))

# Table 2 tissues: the highest, then the second and third only if their TPM is above 1
df_top3 <- long %>%
  group_by(gene) %>%
  slice_max(logTPM, n = 3) %>%
  filter(row_number() == 1 | expressed) %>%
  summarise(top3 = paste0(tissue, collapse = ","),
            top3v = paste0(round(logTPM, 2), collapse = ","))
write_tsv(df_top3, file.path(vd3_dir, "top3_tissues_vd3_genes.tsv"))

# Sensitivity of the Table 2 pattern: other SD and tau cutoffs, tau of untransformed TPM, tau with values below
# 1 TPM set to 0, and the mean log2(TPM + 1) of the four pooled samples per tissue of PRJNA359635 (10 tissues)
patterns <- function(log_tpm, tpm = 2^log_tpm - 1, floor_tpm = 0, ...) {
  sd <- apply(log_tpm, 1, sd)
  tau <- apply(ifelse(tpm < floor_tpm, 0, log_tpm), 1, tau_index)
  expr_pattern(sd, tau, apply(tpm, 1, max), ...)
}
replicated <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"), show_col_types = FALSE) %>%
  filter(study == "PRJNA359635", !excluded) %>%
  group_by(name, tissue) %>% summarise(logTPM = mean(log2(tpm + 1)), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = logTPM)
rep_log <- as.matrix(select(replicated, -name))
rownames(rep_log) <- replicated$name
genes <- rownames(vd3$log_tpm)
tau_tpm <- apply(vd3$tpm, 1, tau_index)
sensitivity <- expr_summary_vd3 %>%
  transmute(Gene = gene, `Peak TPM` = peak_tpm, SD = sd, tau, `Table 2 (SD >= 1 or tau >= 0.8)` = pattern) %>%
  mutate(`SD >= 0.75` = patterns(vd3$log_tpm, sd_cut = 0.75)[Gene],
         `SD >= 1.25` = patterns(vd3$log_tpm, sd_cut = 1.25)[Gene],
         `tau >= 0.7` = patterns(vd3$log_tpm, tau_cut = 0.7)[Gene],
         `tau >= 0.9` = patterns(vd3$log_tpm, tau_cut = 0.9)[Gene],
         `tau of TPM` = expr_pattern(apply(vd3$log_tpm, 1, sd), tau_tpm, apply(vd3$tpm, 1, max))[Gene],
         `tau with TPM below 1 set to 0` = patterns(vd3$log_tpm, floor_tpm = 1)[Gene],
         `PRJNA359635 SD` = round(apply(rep_log, 1, sd), 3)[Gene],
         `PRJNA359635 tau` = round(apply(rep_log, 1, tau_index), 3)[Gene],
         `PRJNA359635 pattern` = patterns(rep_log)[Gene])
changed <- select(sensitivity, `SD >= 0.75`:`tau with TPM below 1 set to 0`, `PRJNA359635 pattern`) %>%
  map(~ sensitivity$Gene[.x != sensitivity$`Table 2 (SD >= 1 or tau >= 0.8)`])
print(changed)
note <- paste("Pattern of each gene as in Table 2 (atlas, 14 tissues, SD and tau of log2(TPM + 1)) and with one",
              "change at a time: other SD or tau cutoffs, tau calculated from TPM, tau after setting values below",
              "1 TPM to 0, or the mean log2(TPM + 1) of the four pooled samples per tissue of PRJNA359635 (10 tissues).",
              "Genes below 1 TPM in all tissues are not classified.")
sensitivity_tsv <- file.path(vd3_dir, "pattern_sensitivity.tsv")
sensitivity %>%
  mutate(across(everything(), as.character)) %>%
  add_row(Gene = c("", note)) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  write_tsv(sensitivity_tsv)

# Table S8: TPM in the analysed atlas tissues. "# Silent" = tissues with a TPM of 0. SD and tau are the
# Table 2 pattern values
table_tsv <- file.path(vd3_dir, "table_s8.tsv")
vd3$info %>%
  transmute(Name, `NCBI Gene ID` = Gene_id, `NCBI Type` = Type, `NCBI Gene Symbol` = Gene_symbol,
            `NCBI Gene Name` = Gene_name, `# Silent` = rowSums(vd3$tpm == 0)) %>%
  left_join(select(expr_summary_vd3, Name = gene, `SD of log2(TPM + 1)` = sd, `tau of log2(TPM + 1)` = tau),
            by = "Name") %>%
  bind_cols(as_tibble(round(vd3$tpm, 6))) %>%
  write_tsv(table_tsv)
table_file <- file.path(fig_dir, "Table S8 TPM values and expression patterns of VD3-related genes across 14 tissues.xlsx")
sheets <- file.path(tempdir(), c("TPM.tsv", "Sensitivity.tsv"))
invisible(file.copy(c(table_tsv, sensitivity_tsv), sheets, overwrite = TRUE))
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(sheets),
                              shQuote(table_file))) == 0)
