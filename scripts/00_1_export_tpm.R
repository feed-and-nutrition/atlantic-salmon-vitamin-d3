# Export gene-level TPM from tximport output into the project.
# Input:  salmon_quant/PRJNA72713/tpm_gene.tsv.gz (exported on the server by 05_1), or else
#         <VD3_TXI_DIR>/txi_gene.rds (15 tissues of PRJNA72713)
# Output: salmon_quant/PRJNA72713/tpm_atlantic_salmon_15_tissues.tsv.gz (all genes)
#         hk_genes/tpm_hk_genes.tsv, vd3_genes/tpm_vd3_genes.tsv (genes of interest)
# All downstream scripts read these tables, so only this script needs the txi object.

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

server_tpm_file <- file.path(quant_dir, "tpm_gene.tsv.gz")
if (file.exists(server_tpm_file)) {
  tpm_server <- read_tsv(server_tpm_file, col_types = cols(gene_id = col_character()))
  expr_gene <- as.matrix(tpm_server[, -1])
  rownames(expr_gene) <- tpm_server$gene_id
} else {
  if (is.na(txi_dir)) stop("Neither ", server_tpm_file, " (05_1) nor VD3_TXI_DIR is available")
  expr_gene <- readRDS(file.path(txi_dir, "txi_gene.rds"))$abundance
}
check_tissues(colnames(expr_gene))
expr_gene <- expr_gene[, expected_tissues]

as_tibble(expr_gene, rownames = "ncbi_gene_id") %>% write_tsv(tpm_all_file)

# Per-gene tables with annotation columns (tissue columns in title case, e.g. Head_kidney)
export_genes <- function(genes, out_file) {
  tpm <- as_tibble(expr_gene, rownames = "gene_id")
  out <- genes %>%
    dplyr::select(name, gene_id, type, gene_symbol, gene_name) %>%
    inner_join(tpm, by = "gene_id")
  colnames(out) <- str_to_title(colnames(out))
  write_tsv(out, out_file)
}

hk_genes <- read_tsv(file.path(hk_dir, "hk_genes.tsv"), col_types = cols(gene_id = col_character())) %>%
  filter(!(gene_id %in% c("100196887", "106576970"))) %>%
  arrange(no)
export_genes(hk_genes, file.path(hk_dir, "tpm_hk_genes.tsv"))

vd3_genes <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(gene_id = col_character())) %>%
  filter(!(name %in% names(excluded_genes))) %>%
  arrange(No)
export_genes(vd3_genes, file.path(vd3_dir, "tpm_vd3_genes.tsv"))
