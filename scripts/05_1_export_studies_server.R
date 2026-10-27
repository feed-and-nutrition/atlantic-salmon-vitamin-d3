# Runs on the machine that holds the Salmon/tximport results, which may not be the analysis machine.
# Exports gene-level TPM and counts for every study in ncbi_samples/sample_info.txt, so that
# all filtering, tissue harmonisation, and plotting can be done locally from the copied results.
# Base R only, so it runs with whatever R is installed there.
#
# Copy to that machine (from the project root):
#   tar czf vd3_server.tgz scripts/05_1_export_studies_server.R ncbi_samples/sample_info.txt salmon_quant/samp_info
# There, in an empty directory, with VD3_MAPPED_DIR set to the directory that holds one directory per study:
#   tar xzf vd3_server.tgz && VD3_MAPPED_DIR=/path/to/quant Rscript scripts/05_1_export_studies_server.R
# Copy back and unpack in the project root:
#   tar czf vd3_results.tgz salmon_quant      (there)
#   tar xzf vd3_results.tgz                   (here)
#
# Input:  $VD3_MAPPED_DIR/<location>/txi_gene.rds for each study (location column of sample_info.txt)
#         salmon_quant/samp_info/<sample_info> sample sheets
#         optional: Salmon aux_info/meta_info.json files under <location>
# Output: salmon_quant/<Project_no>/samples.tsv          sample sheet plus the matching txi column
#         salmon_quant/<Project_no>/tpm_gene.tsv.gz      gene_id x txi columns (TPM)
#         salmon_quant/<Project_no>/counts_gene.tsv.gz   gene_id x txi columns (estimated counts)
#         salmon_quant/<Project_no>/salmon_meta.tsv      reads processed and mapping rate (if meta_info.json found)
#         salmon_quant/export_log.tsv                    one row per study
# VD3_MAPPED_DIR defaults to the working directory. An absolute location is used as it is.

root <- if (file.exists("ncbi_samples/sample_info.txt")) "." else ".."
studies <- read.csv(file.path(root, "ncbi_samples/sample_info.txt"), strip.white = TRUE,
                    stringsAsFactors = FALSE)
mapped_dir <- Sys.getenv("VD3_MAPPED_DIR", unset = ".")
studies$location <- ifelse(startsWith(studies$location, "/"), studies$location, file.path(mapped_dir, studies$location))

write_matrix <- function(m, file) {
  con <- gzfile(file, "w")
  on.exit(close(con))
  write.table(data.frame(gene_id = rownames(m), m, check.names = FALSE),
              con, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Salmon meta_info.json: sample = directory above aux_info. Parsed with regex to avoid jsonlite.
read_salmon_meta <- function(location) {
  files <- list.files(location, pattern = "^meta_info\\.json$", recursive = TRUE, full.names = TRUE)
  files <- files[basename(dirname(files)) == "aux_info"]
  if (!length(files)) return(NULL)
  field <- function(txt, key) {
    hit <- regmatches(txt, regexpr(paste0('"', key, '":\\s*[0-9.eE+-]+'), txt))
    if (length(hit)) as.numeric(sub(".*:\\s*", "", hit)) else NA
  }
  do.call(rbind, lapply(files, function(f) {
    txt <- paste(readLines(f, warn = FALSE), collapse = "")
    data.frame(salmon_dir = basename(dirname(dirname(f))),
               num_processed = field(txt, "num_processed"),
               num_mapped = field(txt, "num_mapped"),
               percent_mapped = field(txt, "percent_mapped"))
  }))
}

export_study <- function(project, sheet_file, location) {
  out_dir <- file.path(root, "salmon_quant", project)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sheet <- read.delim(file.path(root, "salmon_quant/samp_info", sheet_file), colClasses = "character",
                      na.strings = character(0))
  sheet <- sheet[rowSums(sheet != "") > 0, , drop = FALSE]  # drop blank lines
  txi_file <- file.path(location, "txi_gene.rds")
  if (!file.exists(txi_file)) stop("missing ", txi_file)
  txi <- readRDS(txi_file)
  txi_cols <- colnames(txi$abundance)

  # Use the sheet column that names the most txi columns (sample, read_group, or read_file)
  id_cols <- intersect(c("sample", "read_group", "read_file", "lims", "fish_id"), colnames(sheet))
  hits <- sapply(id_cols, function(col) sum(sheet[[col]] %in% txi_cols))
  match_col <- id_cols[which.max(hits)]
  sheet$txi_column <- ifelse(sheet[[match_col]] %in% txi_cols, sheet[[match_col]], NA)
  sheet <- cbind(study = project, sheet)
  write.table(sheet, file.path(out_dir, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  write_matrix(txi$abundance, file.path(out_dir, "tpm_gene.tsv.gz"))
  write_matrix(txi$counts, file.path(out_dir, "counts_gene.tsv.gz"))

  meta <- read_salmon_meta(location)
  if (!is.null(meta)) {
    write.table(meta, file.path(out_dir, "salmon_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  data.frame(study = project, location = location, n_sheet = nrow(sheet), n_txi = length(txi_cols),
             matched_by = match_col, n_matched = sum(!is.na(sheet$txi_column)),
             unmatched_sheet = paste(sheet[[match_col]][is.na(sheet$txi_column)], collapse = ","),
             unmatched_txi = paste(setdiff(txi_cols, sheet$txi_column), collapse = ","),
             n_genes = nrow(txi$abundance), counts_from_abundance = txi$countsFromAbundance,
             n_salmon_meta = if (is.null(meta)) 0 else nrow(meta), status = "ok")
}

log <- do.call(rbind, lapply(seq_len(nrow(studies)), function(i) {
  s <- studies[i, ]
  message("== ", s$Project_no, " (", s$location, ")")
  tryCatch(export_study(s$Project_no, s$sample_info, s$location), error = function(e) {
    message("   FAILED: ", conditionMessage(e))
    data.frame(study = s$Project_no, location = s$location, n_sheet = NA, n_txi = NA, matched_by = NA,
               n_matched = NA, unmatched_sheet = NA, unmatched_txi = NA, n_genes = NA,
               counts_from_abundance = NA, n_salmon_meta = NA, status = conditionMessage(e))
  })
}))
write.table(log, file.path(root, "salmon_quant/export_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(log[, c("study", "n_sheet", "n_txi", "n_matched", "matched_by", "n_salmon_meta", "status")], row.names = FALSE)
