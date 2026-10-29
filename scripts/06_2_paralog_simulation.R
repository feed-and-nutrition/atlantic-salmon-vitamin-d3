# Can Salmon tell the paralogs apart? Mixtures of reads simulated from both genes of a pair, quantified against the
# whole transcriptome. For each of the eight paralog pairs, 100,000 read pairs are simulated with wgsim from the
# transcripts with the longest CDS (those compared in Table S11, 06_1), with the first gene contributing 100%, 90%, 50%,
# 10% or 0% of them: 2 x 100 nt as in the atlas (HiSeq 2000), fragments of 250 +/- 50 nt, 1% sequencing errors and
# 0.5% sequence differences from the reference (wgsim -r, 15% of them indels), as a fish differs from the reference. They are quantified with Salmon as in the Methods (-l A --seqBias --gcBias) against an index
# of all RS_2025_05 transcripts, and the estimated reads are summed per gene.
# Downloads the RS_2025_05 transcripts to annotation/ (git-ignored), where the Salmon index is also built (about
# 30 min on 8 threads, once). Not in run_all. Needs salmon (1.10) and wgsim on the PATH.
# Input:  vd3_genes/paralogs/paralog_similarity.tsv (06_1), annotation/RS_2025_05/ feature table (07_1)
# Output: vd3_genes/paralogs/simulation.tsv (true and estimated share of the first gene among the reads of the pair,
#         and the share of all reads assigned to other genes)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

n_pairs <- 100000L  # integer, so it is not passed as 1e+05
threads <- 8
ann_dir <- file.path(proj_dir, "annotation", "RS_2025_05")
prefix <- "GCF_905237065.1_Ssal_v3.1"
ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/237/065/GCF_905237065.1_Ssal_v3.1"
rna <- file.path(ann_dir, paste0(prefix, "_rna.fna.gz"))
index <- file.path(ann_dir, "salmon_index")
# Runs a command, with its output in log when given (shown if the command fails)
run <- function(cmd, args, log = "") {
  if (system2(cmd, args, stdout = log, stderr = log) != 0) {
    if (nzchar(log)) message(paste(tail(readLines(log), 20), collapse = "\n"))
    stop(cmd, " failed")
  }
}

if (!file.exists(rna)) download.file(file.path(ftp, basename(rna)), rna, mode = "wb")
if (!file.exists(file.path(index, "versionInfo.json"))) {
  run("salmon", c("index", "-t", shQuote(rna), "-i", shQuote(index), "-p", threads))
}

# Transcript -> gene, from the RS_2025_05 feature table (all RNA features with a product accession)
tx2gene <- read_tsv(file.path(ann_dir, paste0(prefix, "_feature_table.txt.gz")), show_col_types = FALSE,
                    col_select = c(`# feature`, product_accession, GeneID)) %>%
  filter(`# feature` != "gene", `# feature` != "CDS", !is.na(product_accession)) %>%
  distinct(transcript = product_accession, gene_id = as.character(GeneID))

pairs <- read_tsv(file.path(vd3_dir, "paralogs", "paralog_similarity.tsv"), show_col_types = FALSE)
ids <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes.tsv"), show_col_types = FALSE) %>%
  transmute(gene = Name, gene_id = as.character(Gene_id))
pairs <- pairs %>%
  select(gene_a, gene_b, transcript_a, transcript_b) %>%
  left_join(rename(ids, gene_a = gene, id_a = gene_id), by = "gene_a") %>%
  left_join(rename(ids, gene_b = gene, id_b = gene_id), by = "gene_b")
stopifnot(!anyNA(pairs))

# Transcript sequences, read once from the reference
fa <- readLines(gzfile(rna))
hdr <- which(startsWith(fa, ">"))
acc <- str_extract(fa[hdr], "^>(\\S+)", group = 1)
seq_of <- function(tx) {
  i <- match(tx, acc)
  stopifnot(!is.na(i))
  end <- if (i < length(hdr)) hdr[i + 1] - 1 else length(fa)
  c(paste0(">", tx), fa[(hdr[i] + 1):end])
}

work <- file.path(tempdir(), "paralog_sim")
dir.create(work, showWarnings = FALSE)
# Simulated read pairs from one transcript, as a list of the two FASTQ files (none when n is 0)
reads_from <- function(d, transcript, n, seed) {
  if (n == 0) return(NULL)
  fasta <- file.path(d, paste0(transcript, ".fa"))
  writeLines(seq_of(transcript), fasta)
  fq <- file.path(d, paste0(transcript, c("_1.fq", "_2.fq")))
  run("wgsim", c("-N", n, "-1", 100, "-2", 100, "-d", 250, "-s", 50, "-e", 0.01, "-r", 0.005, "-R", 0.15,
                 "-S", seed, shQuote(fasta), shQuote(fq[1]), shQuote(fq[2])), log = file.path(d, "wgsim.log"))
  fq
}
simulate <- function(gene_a, gene_b, transcript_a, transcript_b, id_a, id_b, share_a, seed) {
  d <- file.path(work, paste0(gene_a, "_", share_a))
  dir.create(d, showWarnings = FALSE)
  n_a <- as.integer(round(n_pairs * share_a))
  fq <- list(reads_from(d, transcript_a, n_a, seed), reads_from(d, transcript_b, n_pairs - n_a, seed + 1000L))
  fq <- compact(fq)
  r1 <- paste(shQuote(map_chr(fq, 1)), collapse = " ")
  r2 <- paste(shQuote(map_chr(fq, 2)), collapse = " ")
  run("salmon", c("quant", "-i", shQuote(index), "-l", "A", "-1", r1, "-2", r2,
                  "--seqBias", "--gcBias", "-p", threads, "-o", shQuote(file.path(d, "quant"))),
      log = file.path(d, "salmon.log"))
  quant <- read_tsv(file.path(d, "quant", "quant.sf"), show_col_types = FALSE) %>%
    left_join(tx2gene, by = c(Name = "transcript"))
  stopifnot(!anyNA(quant$gene_id))
  reads <- function(id) sum(quant$NumReads[quant$gene_id %in% id])
  total <- sum(quant$NumReads)
  tibble(gene_a, gene_b, true_share_a = 100 * share_a,
         est_share_a = round(100 * reads(id_a) / (reads(id_a) + reads(id_b)), 2),
         other_pct = abs(round(100 * (total - reads(id_a) - reads(id_b)) / total, 2)), assigned = round(total))
}

sim <- pairs %>%
  crossing(share_a = c(1, 0.9, 0.5, 0.1, 0)) %>%
  arrange(match(gene_a, pairs$gene_a), desc(share_a)) %>%
  mutate(seed = row_number()) %>%
  pmap(simulate) %>%
  list_rbind()
write_tsv(sim, file.path(vd3_dir, "paralogs", "simulation.tsv"))
print(as.data.frame(sim))
