# Table S11, sheet Kidney: detection power for cyp27b1 in the kidney and head kidney samples (atlas and replicate
# studies). For each sample: library size (estimated reads), estimated reads and TPM of both cyp27b1 paralogs, and the TPM
# that one read of each paralog corresponds to. Where a sample has no read of the gene, that value is estimated
# from the samples of the same study that have reads, scaled by library size (TPM per read is about inversely
# proportional to the number of reads).
# Sheet Freshwater and seawater: TPM of the VD3-related genes in the liver and kidney of PRJNA359635 (tracker item
# 38), one pool of three fish per sex before (freshwater) and after (seawater) smoltification, with the difference of
# the mean log2(TPM + 1).
# Sheet Genome search: the mitochondrial P450 loci found by the CYP27B1 genome search of 07_5 (not in run_all).
# Needs python3 for the Excel file.
# Input:  salmon_quant/samples_all.tsv (05_2), salmon_quant/<study>/{tpm_gene,counts_gene}.tsv.gz,
#         vd3_genes/tpm_vd3_genes.tsv (gene IDs), vd3_genes/tpm_vd3_genes_studies.tsv (05_2),
#         vd3_genes/cyp27b1_search/loci.tsv (07_5), msa/cyp27a1/cyp27a1_accessions.tsv (names of the cyp27a2/a3-like genes)
# Output: vd3_genes/studies/{kidney_cyp27b1,freshwater_seawater}.tsv, manuscript/v6/Table S11 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_file <- file.path(vd3_dir, "studies", "kidney_cyp27b1.tsv")
fw_sw_file <- file.path(vd3_dir, "studies", "freshwater_seawater.tsv")
table_file <- file.path(fig_dir, paste("Table S11 Detection of cyp27b1 in kidney samples, CYP27B1 genome search, and",
                                       "VD3-related genes in freshwater and seawater fish.xlsx"))

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
sheet_dir <- tempfile()
dir.create(sheet_dir)
sheets <- file.path(sheet_dir, c("Kidney.tsv", "Genome search.tsv", "Freshwater and seawater.tsv"))
write_tsv(table, sheets[1], na = "")
note <- paste("Reads are Salmon estimated counts. TPM per read is the TPM that one read of the gene gives in the sample.",
              "Where the sample has no read, it is estimated from the samples of the same study with reads, scaled by",
              "library size. The PRJNA647285 samples labelled kidney were analysed as head kidney (Methods).")
cat("\n", note, "\n", file = sheets[1], append = TRUE, sep = "")

# --- Genome search (07_5): loci of the mitochondrial P450 families, with the names used in this study
names_used <- bind_rows(
  read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
    transmute(gene_id, name),
  read_tsv(file.path(proj_dir, "msa", "cyp27a1", "cyp27a1_accessions.tsv"), col_types = cols(.default = "c")) %>%
    filter(str_detect(label, "^Salmon_LOC")) %>%
    transmute(gene_id, name = paste0(str_remove(label, "^Salmon_"), " (", str_extract(note, "cyp27a[23]"), "-like)")))
read_tsv(file.path(vd3_dir, "cyp27b1_search", "loci.tsv"), col_types = cols(chr = "c", .default = col_guess())) %>%
  filter(class != "Other P450") %>%
  mutate(assembly = factor(assembly, levels = c("Ssal_v3.1", "fSalSal1", "fEsoLuc1")),
         ids = str_extract_all(coalesce(genes, ""), "\\((\\d+),", simplify = FALSE) %>% map(~ str_extract(.x, "\\d+")),
         used = map_chr(ids, ~ paste(na.omit(names_used$name[match(.x, names_used$gene_id)]), collapse = "; "))) %>%
  arrange(assembly, desc(class == "CYP27B1"), class, chr, start) %>%
  transmute(Assembly = as.character(assembly), Species = species,
            Chromosome = if_else(is.na(chr), "unplaced", if_else(species == "Atlantic salmon", paste0("ssa", str_remove(chr, "^ssa")), chr)),
            Start = start, End = end, `Best tBLASTn query` = best_query, `Query bit score` = best_query_bits,
            Family = class, `Family bit score` = round(class_bits, 1), `Second family` = second_class,
            `Second family bit score` = round(second_bits, 1), `Annotated genes` = coalesce(genes, "none"),
            `Gene in this study` = used) %>%
  write_tsv(sheets[2], na = "")
note <- paste("Human, zebrafish, spotted gar and northern pike CYP27B1 were searched with tBLASTn (E <= 1e-5) against",
              "the salmon assemblies of RS_2025_05 (Ssal_v3.1) and RS_2026_09 (fSalSal1) and the northern pike genome",
              "(one cyp27b1 gene, no Ss4R duplication). Hits within 50 kb form a locus, which was assigned by BLASTx to the",
              "mitochondrial P450 family with the highest summed bit score. Loci assigned to microsomal P450s are not shown.")
cat("\n", note, "\n", file = sheets[2], append = TRUE, sep = "")

# --- Freshwater and seawater (PRJNA359635)
fw_sw <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"), show_col_types = FALSE) %>%
  filter(study == "PRJNA359635", !excluded, tissue %in% c("liver", "kidney"), !name %in% excluded_genes) %>%
  mutate(group = paste(str_to_sentence(str_replace(stage, "saltwater", "seawater")), if_else(sex == "F", "female", "male")))
stopifnot(count(fw_sw, name, tissue)$n == 4)
fw_sw <- fw_sw %>%
  group_by(name, tissue) %>%
  mutate(log2fc = mean(log2(tpm[stage == "saltwater"] + 1)) - mean(log2(tpm[stage == "freshwater"] + 1))) %>%
  ungroup() %>%
  select(name, tissue, group, tpm, log2fc) %>%
  pivot_wider(names_from = group, values_from = tpm) %>%
  mutate(tissue = factor(tissue, levels = c("liver", "kidney")),
         name = factor(name, levels = pathway_order(unique(name)))) %>%
  arrange(tissue, name)
write_tsv(fw_sw, fw_sw_file)
fw_sw %>%
  transmute(Tissue = tissue_label(as.character(tissue)), Gene = as.character(name),
            across(c(`Freshwater female`, `Freshwater male`, `Seawater female`, `Seawater male`), ~ fmt(.x, 2),
                   .names = "{.col} (TPM)"),
            `Difference in mean log2(TPM + 1), seawater minus freshwater` = fmt(log2fc, 2)) %>%
  write_tsv(sheets[3])
note <- paste("PRJNA359635 (Robertson et al. 2017): each sample is a pool of three fish of one sex, sampled in",
              "freshwater before smoltification (mean length 18.6 cm) or in seawater after smoltification (25.8 cm).")
cat("\n", note, "\n", file = sheets[3], append = TRUE, sep = "")
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(sheets),
                              shQuote(table_file))) == 0)
print(as.data.frame(filter(fw_sw, tissue == "liver")))
message("Wrote ", out_file, ", ", fw_sw_file, " and ", table_file)
