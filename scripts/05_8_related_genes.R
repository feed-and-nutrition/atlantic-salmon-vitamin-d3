# Table S12: expression of genes related to the VD3 pathway but outside the 16 genes (tracker item 45).
# Groups:
#   25-hydroxylase candidates: CYP3A4 and CYP2J2, which also 25-hydroxylate vitamin D3 in mammals, the four
#     cyp27a2/a3-like genes named "sterol 26-hydroxylase" (Fig 2a) and LOC106586122 (cyp2x, excluded_genes)
#   Cholesterol synthesis: HMGCR and the post-squalene enzymes, to place dhcr7 against the pathway that makes 7-DHC
#   Renal uptake: LRP2 (megalin), CUBN (cubilin) and AMN (amnionless), the receptor complex that takes up GC-bound
#     25(OH)D3 in the mammalian proximal tubule
# Salmon genes: the human protein (longest isoform) was searched with BLASTp (E <= 1e-10) against the RS_2025_05
# proteome, and a hit gene was kept when its longest protein had its best human hit in the query gene (any CYP3A for
# CYP3A4). The search is cached in vd3_genes/related_genes/genes.tsv. Without it, it needs blastp and makeblastdb and
# the RS_2025_05 and GRCh38.p14 protein and feature tables in annotation/ (07_1, 06_2).
# Expression: atlas TPM (14 tissues) and the replicate studies pooled per tissue as in Fig 4 (mean log2(TPM + 1) per
# study, then the mean of the studies), given back as TPM.
# Needs python3 for the Excel file.
# Input:  salmon_quant/PRJNA72713/tpm_atlantic_salmon_15_tissues.tsv.gz (00_1), salmon_quant/samples_all.tsv (05_2),
#         salmon_quant/<study>/tpm_gene.tsv.gz, msa/cyp27a1/cyp27a1_accessions.tsv, vd3_genes/salmon_vd3_orthologs_v1.tsv
# Output: vd3_genes/related_genes/{genes,atlas_tpm,studies_tpm}.tsv, manuscript/v6/Table S12 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_dir <- file.path(vd3_dir, "related_genes")
dir.create(out_dir, showWarnings = FALSE)
genes_file <- file.path(out_dir, "genes.tsv")
table_file <- file.path(fig_dir, paste("Table S12 Expression of 25-hydroxylase candidates, cholesterol synthesis genes,",
                                       "and renal uptake receptors.xlsx"))

queries <- tribble(
  ~group,                       ~human,    ~accept,
  "25-hydroxylase candidates",  "CYP3A4",  "CYP3A\\d+",
  "25-hydroxylase candidates",  "CYP2J2",  "CYP2J2",
  "25-hydroxylase candidates",  "CYP27A1", "CYP27A1",
  "Cholesterol synthesis",      "HMGCR",   "HMGCR",
  "Cholesterol synthesis",      "SQLE",    "SQLE",
  "Cholesterol synthesis",      "LSS",     "LSS",
  "Cholesterol synthesis",      "CYP51A1", "CYP51A1",
  "Cholesterol synthesis",      "TM7SF2",  "TM7SF2",
  "Cholesterol synthesis",      "LBR",     "LBR",
  "Cholesterol synthesis",      "MSMO1",   "MSMO1",
  "Cholesterol synthesis",      "NSDHL",   "NSDHL",
  "Cholesterol synthesis",      "HSD17B7", "HSD17B7",
  "Cholesterol synthesis",      "EBP",     "EBP",
  "Cholesterol synthesis",      "SC5D",    "SC5D",
  "Cholesterol synthesis",      "DHCR7",   "DHCR7",
  "Cholesterol synthesis",      "DHCR24",  "DHCR24",
  "Renal uptake",               "LRP2",    "LRP2",
  "Renal uptake",               "CUBN",    "CUBN",
  "Renal uptake",               "AMN",     "AMN")

vd3 <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
  select(name, gene_id)

# --- Salmon genes (reciprocal BLASTp against human)
if (!file.exists(genes_file)) {
  ann_dir <- file.path(proj_dir, "annotation")
  files <- c(salmon = file.path(ann_dir, "RS_2025_05", "GCF_905237065.1_Ssal_v3.1"),
             human = file.path(ann_dir, "human_GRCh38.p14", "GCF_000001405.40_GRCh38.p14"))
  stopifnot(file.exists(paste0(files, "_protein.faa.gz")), file.exists(paste0(files, "_feature_table.txt.gz")))
  read_features <- function(prefix) {
    file <- paste0(prefix, "_feature_table.txt.gz")
    cols <- str_split(str_remove(readLines(file, n = 1), "^# "), "\t")[[1]]
    read_tsv(file, skip = 1, col_names = cols, col_types = cols(.default = "c"),
             col_select = c(feature, class, chromosome, product_accession, name, symbol, GeneID))
  }
  read_fasta <- function(file) {
    lines <- readLines(file)
    head <- grepl("^>", lines)
    tibble(protein = str_extract(lines[head], "^>(\\S+)", group = 1),
           seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
  }
  ft <- map(files, read_features)
  gene_info <- map(ft, ~ .x %>% filter(feature == "gene", class == "protein_coding") %>%
                     transmute(gene_id = GeneID, symbol, product = name, chr = chromosome) %>%
                     distinct(gene_id, .keep_all = TRUE))  # human genes repeat on the alternate loci
  prot <- map(names(files), function(sp) {
    read_fasta(paste0(files[[sp]], "_protein.faa.gz")) %>%
      inner_join(distinct(filter(ft[[sp]], feature == "CDS", !is.na(product_accession)),
                          protein = product_accession, gene_id = GeneID), by = "protein")
  }) %>% set_names(names(files))
  longest <- map(prot, ~ .x %>% group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup())
  db_dir <- file.path(ann_dir, "blastdb")
  make_db <- function(sp) {
    db <- file.path(db_dir, paste0("related_", sp))
    if (!file.exists(paste0(db, ".pin"))) {
      fa <- tempfile(fileext = ".faa")
      writeLines(paste0(">", prot[[sp]]$protein, "\n", prot[[sp]]$seq), fa)
      stopifnot(system2("makeblastdb", c("-in", fa, "-dbtype", "prot", "-out", shQuote(db)), stdout = FALSE) == 0)
    }
    db
  }
  blastp <- function(query, sp, max_target) {
    fa <- tempfile(fileext = ".faa")
    writeLines(paste0(">", query$protein, "\n", query$seq), fa)
    out <- system2("blastp", c("-query", fa, "-db", shQuote(make_db(sp)), "-evalue", "1e-10",
                               "-max_target_seqs", max_target, "-num_threads", 4,
                               "-outfmt", shQuote("6 qseqid sseqid pident bitscore")), stdout = TRUE)
    read_tsv(I(paste(c("qseqid\tsseqid\tpident\tbitscore", out), collapse = "\n")), col_types = "ccdd") %>%
      inner_join(select(prot[[sp]], sseqid = protein, gene_id), by = "sseqid")
  }
  human_query <- longest$human %>%
    inner_join(select(gene_info$human, gene_id, human = symbol), by = "gene_id") %>%
    filter(human %in% queries$human)
  stopifnot(setequal(human_query$human, queries$human))
  fwd <- blastp(human_query, "salmon", 500) %>%
    inner_join(select(human_query, qseqid = protein, human), by = "qseqid") %>%
    distinct(human, gene_id)
  cyp2x <- excluded_genes[["cyp2r1.c"]]
  rev <- blastp(filter(longest$salmon, gene_id %in% c(fwd$gene_id, cyp2x)), "human", 5) %>%
    group_by(qseqid) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    inner_join(select(gene_info$human, gene_id, best_human = symbol), by = "gene_id") %>%
    inner_join(select(longest$salmon, qseqid = protein, gene_id), by = "qseqid", suffix = c("_human", "")) %>%
    transmute(gene_id, best_human, identity = round(pident, 1))
  cyp27 <- read_tsv(file.path(msa_dir, "cyp27a1", "cyp27a1_accessions.tsv"), col_types = cols(.default = "c")) %>%
    filter(str_detect(note, "closest zebrafish protein")) %>%
    transmute(gene_id, note = paste0(str_extract(note, "cyp27a[23]$"), "-like (Figure 2a)"))
  # LOC106586122 (cyp2x) is added by hand, because its best human hit is not a query gene
  genes <- fwd %>%
    inner_join(queries, by = "human") %>%
    bind_rows(tibble(group = "25-hydroxylase candidates", human = "CYP2R1", gene_id = cyp2x, accept = ".*")) %>%
    inner_join(rev, by = "gene_id") %>%
    filter(str_detect(best_human, paste0("^(", accept, ")$"))) %>%
    left_join(gene_info$salmon, by = "gene_id") %>%
    left_join(vd3, by = "gene_id") %>%
    left_join(cyp27, by = "gene_id") %>%
    # The two cyp27a1 orthologs are among the 16 genes
    filter(human != "CYP27A1" | is.na(name)) %>%
    mutate(note = if_else(gene_id == cyp2x, "cyp2x (Figure 2c), excluded from the 16 genes", note),
           name = if_else(is.na(name) | gene_id == cyp2x, symbol, name)) %>%
    select(group, query = human, name, gene_id, chr, product, best_human, identity, note) %>%
    arrange(match(group, unique(queries$group)), match(query, c(queries$human, "CYP2R1")), chr, name)
  stopifnot(!anyDuplicated(genes$gene_id))
  write_tsv(genes, genes_file, na = "")
}
genes <- read_tsv(genes_file, col_types = cols(.default = "c"))

# --- Expression
atlas <- read_tsv(file.path(quant_dir, "tpm_atlantic_salmon_15_tissues.tsv.gz"),
                  col_types = cols(ncbi_gene_id = col_character())) %>%
  rename(gene_id = ncbi_gene_id) %>%
  filter(gene_id %in% genes$gene_id) %>%
  select(gene_id, all_of(analysis_tissues))
stopifnot(nrow(atlas) == nrow(genes))
samples <- read_tsv(file.path(proj_dir, "salmon_quant", "samples_all.tsv"), show_col_types = FALSE) %>%
  filter(study %in% replicate_studies, !excluded)
studies <- map_dfr(replicate_studies, function(s) {
  read_tsv(file.path(proj_dir, "salmon_quant", s, "tpm_gene.tsv.gz"),
           col_types = cols(gene_id = col_character(), .default = col_double())) %>%
    filter(gene_id %in% genes$gene_id) %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "tpm") %>%
    mutate(id = paste(s, sample, sep = ":"))
}) %>%
  inner_join(select(samples, id, study, tissue), by = "id") %>%
  group_by(gene_id, tissue, study) %>%
  summarise(mean_log = mean(log2(tpm + 1)), .groups = "drop") %>%
  group_by(gene_id, tissue) %>%
  summarise(tpm = 2^mean(mean_log) - 1, .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = tpm) %>%
  select(gene_id, any_of(tissue_order))
write_tsv(mutate(atlas, across(-gene_id, ~ round(.x, 3))), file.path(out_dir, "atlas_tpm.tsv"))
write_tsv(mutate(studies, across(-gene_id, ~ round(.x, 3))), file.path(out_dir, "studies_tpm.tsv"))

# --- Table S12
fmt <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 1))
sheets <- file.path(tempdir(), c("Genes.tsv", "Atlas TPM.tsv", "Studies TPM.tsv"))
key <- genes %>% transmute(Group = group, Gene = name, `Gene ID` = gene_id)
genes %>%
  transmute(Group = group, `Human query` = query, Gene = name, `Gene ID` = gene_id, Chromosome = chr,
            `RefSeq description` = product, `Best human BLASTp hit` = best_human, `Identity (%)` = identity,
            Note = note) %>%
  write_tsv(sheets[1], na = "")
note <- paste("Human proteins (longest isoform) were searched with BLASTp (E <= 1e-10) against the RS_2025_05 proteome.",
              "A salmon gene was kept when its longest protein had its best human hit in the query gene (any CYP3A for",
              "CYP3A4). The cyp27a1 and dhcr7 orthologs among the 16 genes are listed for comparison where they are hits.")
cat("\n", note, "\n", file = sheets[1], append = TRUE, sep = "")
key %>%
  left_join(atlas, by = c(`Gene ID` = "gene_id")) %>%
  mutate(across(any_of(analysis_tissues), fmt)) %>%
  rename_with(tissue_label, any_of(analysis_tissues)) %>%
  write_tsv(sheets[2], na = "")
cat("\nTPM in the 14 tissues of the PRJNA72713 atlas (heart excluded).\n", file = sheets[2], append = TRUE)
n_studies <- samples %>% distinct(study, tissue) %>% count(tissue)
key %>%
  left_join(studies, by = c(`Gene ID` = "gene_id")) %>%
  mutate(across(any_of(tissue_order), fmt)) %>%
  rename_with(~ paste0(tissue_label(.x), " (", n_studies$n[match(.x, n_studies$tissue)], ")"), any_of(tissue_order)) %>%
  write_tsv(sheets[3], na = "")
note <- paste("TPM in the seven replicate studies, pooled per tissue as in Figure 4: the mean of log2(TPM + 1) over the",
              "samples of each study, averaged over the studies (number in brackets) and given back as TPM.")
cat("\n", note, "\n", file = sheets[3], append = TRUE, sep = "")
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(sheets),
                              shQuote(table_file))) == 0)
print(as.data.frame(select(genes, group, query, name, chr, best_human, identity, note)))
message("Wrote ", genes_file, " and ", table_file)
