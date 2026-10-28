# Northern pike as the reference without the salmonid genome duplication (Ss4R): tissue expression of the pike
# orthologs of the VD3-related genes, from the tissue RNA-seq of the pike genome project (BioProject PRJNA221548,
# Rondeau et al. 2014: 13 tissues, one sample each, 2 x 100 nt, HiSeq 2000). Its tissues are those of the salmon atlas,
# including kidney, which the PhyloFish pike samples (PRJNA256965) lack.
# The first 4 million read pairs of each run are quantified with Salmon (-l A --seqBias --gcBias) against an index of
# all fEsoLuc1 RefSeq transcripts, and the TPM of the transcripts are summed per gene. 4 million pairs give about 10
# reads per gene at 1 TPM, enough for the profiles compared here.
# Pike orthologs: the best BLASTp hit of each salmon protein (the protein of the longest CDS, cached by 06_1) in the
# fEsoLuc1 proteome (BLAST database of 07_3). Both genes of each salmon pair must hit the same pike gene.
# Downloads about 9 GB (reads, kept in annotation/pike_rnaseq/, git-ignored) and the pike transcripts, and builds the
# Salmon index there. Not in run_all. Needs fastq-dump (sra-tools), salmon (1.10) and blastp on the PATH.
# Input:  vd3_genes/paralogs/seq/ (06_1), annotation/gc_search/ pike feature table and annotation/blastdb/ (07_3)
# Output: vd3_genes/paralogs/pike_orthologs.tsv, vd3_genes/paralogs/pike_tpm.tsv (pike orthologs, TPM per tissue),
#         annotation/pike_rnaseq/pike_gene_tpm.tsv.gz (all pike genes)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

n_reads <- 4000000L  # integer, so it is not passed as 4e+06
threads <- 8
downloads <- 4  # runs downloaded at the same time
bioproject <- "PRJNA221548"
work_dir <- file.path(proj_dir, "annotation", "pike_rnaseq")
for (d in c("reads", "quant")) dir.create(file.path(work_dir, d), recursive = TRUE, showWarnings = FALSE)
prefix <- "GCF_011004845.1_fEsoLuc1.pri"
ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/004/845/GCF_011004845.1_fEsoLuc1.pri"
feature_table <- file.path(proj_dir, "annotation", "gc_search", paste0(prefix, "_feature_table.txt.gz"))
proteome_db <- file.path(proj_dir, "annotation", "blastdb", "fEsoLuc1_proteome")
rna <- file.path(work_dir, paste0(prefix, "_rna.fna.gz"))
index <- file.path(work_dir, "salmon_index")
out_dir <- file.path(vd3_dir, "paralogs")
stopifnot(file.exists(feature_table), file.exists(paste0(proteome_db, ".pin")))
# Runs a command, with its output in log when given (shown if the command fails)
run <- function(cmd, args, log = "") {
  if (system2(cmd, args, stdout = log, stderr = log) != 0) {
    if (nzchar(log)) message(paste(tail(readLines(log), 20), collapse = "\n"))
    stop(cmd, " failed")
  }
}

# --- Runs: SRA run info, cached as for the salmon studies (05_6)
runinfo_file <- file.path(proj_dir, "ncbi_samples", "runinfo", paste0(bioproject, ".csv"))
if (!file.exists(runinfo_file)) {
  download.file(paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=", bioproject), runinfo_file,
                quiet = TRUE)
}
runs <- read_csv(runinfo_file, show_col_types = FALSE) %>%
  filter(LibraryStrategy == "RNA-Seq", ScientificName == "Esox lucius") %>%
  transmute(run = Run, tissue = str_to_sentence(str_replace(LibraryName, " ", "_")))
stopifnot(nrow(runs) == 13, !anyDuplicated(runs$tissue))

# --- Salmon index of the pike transcripts
if (!file.exists(rna)) download.file(file.path(ftp, basename(rna)), rna, mode = "wb")
if (!file.exists(file.path(index, "versionInfo.json"))) {
  run("salmon", c("index", "-t", shQuote(rna), "-i", shQuote(index), "-p", threads),
      file.path(work_dir, "salmon_index.log"))
}

# --- Reads and quantification
reads <- function(r) file.path(work_dir, "reads", paste0(r, "_", 1:2, ".fastq.gz"))
missing <- runs$run[!map_lgl(runs$run, ~ all(file.exists(reads(.x))))]
invisible(parallel::mclapply(missing, function(r) {
  tmp <- file.path(work_dir, "reads", paste0(r, ".part"))
  dir.create(tmp, showWarnings = FALSE)
  run("fastq-dump", c("-X", n_reads, "--split-files", "--gzip", "-O", shQuote(tmp), r), file.path(tmp, "log"))
  file.rename(file.path(tmp, paste0(r, "_", 1:2, ".fastq.gz")), reads(r))
}, mc.cores = downloads))
for (r in runs$run) {
  quant <- file.path(work_dir, "quant", r)
  if (!file.exists(file.path(quant, "quant.sf"))) {
    run("salmon", c("quant", "-i", shQuote(index), "-l", "A", "--seqBias", "--gcBias", "-p", threads,
                    "-1", shQuote(reads(r)[1]), "-2", shQuote(reads(r)[2]), "-o", shQuote(quant)),
        paste0(quant, ".log"))
  }
}

# Transcript -> gene, from the feature table (all RNA features with a product accession)
ft <- read_tsv(feature_table, show_col_types = FALSE,
               col_select = c(`# feature`, product_accession, symbol, GeneID, name))
tx2gene <- ft %>%
  filter(!`# feature` %in% c("gene", "CDS"), !is.na(product_accession)) %>%
  distinct(transcript = product_accession, gene_id = as.character(GeneID))
gene_tpm <- map_dfr(set_names(runs$run, runs$tissue), function(r) {
  read_tsv(file.path(work_dir, "quant", r, "quant.sf"), show_col_types = FALSE) %>%
    inner_join(tx2gene, by = c(Name = "transcript")) %>%
    group_by(gene_id) %>%
    summarise(tpm = sum(TPM), reads = sum(NumReads), .groups = "drop")
}, .id = "tissue")
gene_tpm %>%
  select(-reads) %>%
  pivot_wider(names_from = tissue, values_from = tpm) %>%
  write_tsv(file.path(work_dir, "pike_gene_tpm.tsv.gz"))
mapping <- map_dfr(runs$run, function(r) {
  info <- jsonlite::read_json(file.path(work_dir, "quant", r, "aux_info", "meta_info.json"))
  tibble(run = r, processed = info$num_processed, mapping_rate = round(info$percent_mapped, 1))
})
print(as.data.frame(left_join(runs, mapping, by = "run")))

# --- Pike orthologs: best BLASTp hit of each salmon protein
read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(header = sub("^>", "", lines[head]),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}
seq_dir <- file.path(out_dir, "seq")
genes <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
  transmute(name, transcript = str_split(str_trim(`mRNA refseq`), ",")) %>%
  unnest(transcript)
proteins <- read_fasta(file.path(seq_dir, "protein.faa")) %>%
  transmute(transcript = str_match(header, "^lcl\\|(.+?)_prot_")[, 2], protein = seq) %>%
  inner_join(genes, by = "transcript") %>%
  group_by(name) %>%
  slice_max(nchar(protein), n = 1, with_ties = FALSE) %>%
  ungroup()
query <- tempfile(fileext = ".faa")
writeLines(paste0(">", proteins$name, "\n", proteins$protein), query)
hits <- read.table(text = system2("blastp", c("-query", query, "-db", shQuote(proteome_db), "-evalue", "1e-10",
                                              "-max_target_seqs", "5", "-outfmt",
                                              shQuote("6 qseqid sseqid pident length bitscore")), stdout = TRUE),
                   col.names = c("name", "protein_acc", "identity", "length", "bits")) %>%
  as_tibble() %>%
  mutate(protein_acc = str_remove_all(protein_acc, "^ref\\||\\|$"))
products <- ft %>%
  filter(`# feature` == "CDS", !is.na(product_accession)) %>%
  distinct(protein_acc = product_accession, pike_gene_id = as.character(GeneID), pike_symbol = symbol,
           pike_product = name)
orthologs <- hits %>%
  inner_join(products, by = "protein_acc") %>%
  group_by(name) %>%
  slice_max(bits, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(name, pike_gene_id, pike_symbol, pike_protein = protein_acc, pike_product = str_remove(pike_product, " isoform X\\d+$"),
            identity = round(identity, 1), bits)
write_tsv(orthologs, file.path(out_dir, "pike_orthologs.tsv"))
print(as.data.frame(orthologs))
# Both copies of each Ss4R pair point to one pike gene
pair_of <- function(name) str_remove(name, "\\.[ab]$")
ss4r <- filter(orthologs, str_detect(name, "\\.[ab]$"))
stopifnot(setequal(orthologs$name, proteins$name),
          all(tapply(ss4r$pike_gene_id, pair_of(ss4r$name), n_distinct) == 1))

gene_tpm %>%
  filter(gene_id %in% orthologs$pike_gene_id) %>%
  select(pike_gene_id = gene_id, tissue, tpm) %>%
  mutate(tpm = round(tpm, 3)) %>%
  pivot_wider(names_from = tissue, values_from = tpm) %>%
  right_join(distinct(orthologs, pike_gene_id, pike_symbol), ., by = "pike_gene_id") %>%
  write_tsv(file.path(out_dir, "pike_tpm.tsv"))
