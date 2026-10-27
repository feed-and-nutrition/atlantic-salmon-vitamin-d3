# Tables for VD3-related gene expression (basis of Table 2)
# Input:  vd3_genes/tpm_vd3_genes.tsv (from 00_1)
# Output: vd3_genes/summary_vd3_genes.tsv (level, SD, tau, Table 2 pattern per gene)
#         vd3_genes/top3_tissues_vd3_genes.tsv (highest tissue, and the next two if their TPM is above 1)
#         manuscript/v6/Table S8 ... .xlsx (TPM per gene and atlas tissue with SD and tau, via vd3_genes/table_s8.tsv)

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
            tau = round(tau_index(logTPM), 2)) %>%
  mutate(stability = ifelse(sd < 1, "Stable", "Variable"),
         pattern = expr_pattern(sd, tau))
write_tsv(expr_summary_vd3, file.path(vd3_dir, "summary_vd3_genes.tsv"))

# Table 2 tissues: the highest, then the second and third only if their TPM is above 1
df_top3 <- long %>%
  group_by(gene) %>%
  slice_max(logTPM, n = 3) %>%
  filter(row_number() == 1 | expressed) %>%
  summarise(top3 = paste0(tissue, collapse = ","),
            top3v = paste0(round(logTPM, 2), collapse = ","))
write_tsv(df_top3, file.path(vd3_dir, "top3_tissues_vd3_genes.tsv"))

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
table_file <- file.path(fig_dir, "Table S8 TPM values of VD3-related genes across 14 tissues.xlsx")
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(table_tsv),
                              shQuote(table_file))) == 0)
