# Table S11, sheet Similarity: can Salmon tell the paralogs apart? Sequence identity and shared read-length windows of paralog pairs.
# Sequences: the RefSeq transcripts (RS_2025_05 versions, as quantified by Salmon) listed in
# salmon_vd3_orthologs_v1.tsv, with their CDS and proteins, from NCBI E-utilities. Cached in vd3_genes/paralogs/seq/
# (delete that folder to download again). Each gene is represented by its transcript with the longest CDS.
# Identity: best blastn (CDS) or blastp (protein) hit, with its coverage of the shorter sequence.
# Shared windows: share of the 100-nt windows of each gene's transcript that occur exactly in any transcript of
# its paralog. Reads from these windows cannot be assigned by sequence (100 nt is conservative for these data,
# and read pairs span longer fragments).
# Needs blastn and blastp on the PATH and internet access on the first run.
# Input:  vd3_genes/salmon_vd3_orthologs_v1.tsv, vd3_genes/tpm_vd3_genes.tsv, vd3_genes/paralogs/simulation.tsv (06_2)
# Output: vd3_genes/paralogs/paralog_similarity.tsv, vd3_genes/paralogs/table_s9.tsv (the sheet, written to the
#         Excel file by 06_4)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_dir <- file.path(vd3_dir, "paralogs")
seq_dir <- file.path(out_dir, "seq")
dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- tribble(~gene_a, ~gene_b,
                 "dhcr7.a", "dhcr7.b",
                 "alb1", "alb2",
                 "cyp2r1.a", "cyp2r1.b",
                 "cyp27a1.a", "cyp27a1.b",
                 "cyp27b1.a", "cyp27b1.b",
                 "cyp24a1.a", "cyp24a1.b",
                 "vdra.a", "vdra.b",
                 "vdrb.a", "vdrb.b")

# Transcripts of each gene in the Salmon reference
transcripts <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
  transmute(name, transcript = str_split(str_trim(`mRNA refseq`), ",")) %>%
  unnest(transcript) %>%
  filter(name %in% c(pairs$gene_a, pairs$gene_b))

# --- Sequences
read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(header = sub("^>", "", lines[head]),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}

fasta_files <- set_names(file.path(seq_dir, c("rna.fna", "cds.fna", "protein.faa")), c("fasta", "fasta_cds_na", "fasta_cds_aa"))
for (type in names(fasta_files)) {
  if (!file.exists(fasta_files[[type]])) {
    download.file(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=text&rettype=",
                         type, "&id=", paste(transcripts$transcript, collapse = ",")),
                  fasta_files[[type]], quiet = TRUE)
  }
}

rna <- read_fasta(fasta_files[["fasta"]]) %>%
  transmute(transcript = str_extract(header, "^\\S+"), rna = seq) %>%
  inner_join(transcripts, by = "transcript")
coding <- function(file, col) {
  read_fasta(file) %>% transmute(transcript = str_match(header, "^lcl\\|(.+?)_(?:cds|prot)_")[, 2], !!col := seq)
}
# One transcript per gene: the one with the longest CDS
rep_seq <- coding(fasta_files[["fasta_cds_na"]], "cds") %>%
  inner_join(coding(fasta_files[["fasta_cds_aa"]], "protein"), by = "transcript") %>%
  inner_join(rna, by = "transcript") %>%
  group_by(name) %>%
  slice_max(nchar(cds), n = 1, with_ties = FALSE) %>%
  ungroup()
stopifnot(setequal(rep_seq$name, c(pairs$gene_a, pairs$gene_b)), setequal(rna$transcript, transcripts$transcript))

blast_identity <- function(a, b, program) {
  fa <- tempfile(fileext = ".fa")
  fb <- tempfile(fileext = ".fa")
  writeLines(c(">a", a), fa)
  writeLines(c(">b", b), fb)
  args <- c("-query", fa, "-subject", fb, "-outfmt", shQuote("6 pident qstart qend sstart send bitscore"),
            "-max_hsps", "1")
  if (program == "blastn") args <- c(args, "-task", "blastn")
  hit <- read.table(text = system2(program, args, stdout = TRUE),
                    col.names = c("pident", "qstart", "qend", "sstart", "send", "bits"))
  hit <- hit[which.max(hit$bits), ]
  covered <- if (nchar(a) <= nchar(b)) c(hit$qend - hit$qstart + 1, nchar(a)) else c(abs(hit$send - hit$sstart) + 1, nchar(b))
  c(identity = hit$pident, coverage = 100 * covered[1] / covered[2])
}

seq_of <- function(gene, type) rep_seq[[type]][rep_seq$name == gene]

windows <- function(s, k = 100) substring(s, 1:(nchar(s) - k + 1), k:nchar(s))
shared_windows <- function(gene, paralog) {
  other <- unique(unlist(map(filter(rna, name == paralog)$rna, windows)))
  100 * mean(windows(seq_of(gene, "rna")) %in% other)
}

atlas <- read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))
peak <- function(gene) colnames(atlas$log_tpm)[which.max(atlas$log_tpm[gene, ])]

result <- pairs %>%
  mutate(transcript_a = map_chr(gene_a, seq_of, "transcript"),
         transcript_b = map_chr(gene_b, seq_of, "transcript"),
         protein = map2(gene_a, gene_b, ~ blast_identity(seq_of(.x, "protein"), seq_of(.y, "protein"), "blastp")),
         cds = map2(gene_a, gene_b, ~ blast_identity(seq_of(.x, "cds"), seq_of(.y, "cds"), "blastn")),
         protein_identity = map_dbl(protein, "identity"), protein_coverage = map_dbl(protein, "coverage"),
         cds_identity = map_dbl(cds, "identity"), cds_coverage = map_dbl(cds, "coverage"),
         shared_a = map2_dbl(gene_a, gene_b, shared_windows), shared_b = map2_dbl(gene_b, gene_a, shared_windows),
         atlas_peak_a = map_chr(gene_a, peak), atlas_peak_b = map_chr(gene_b, peak)) %>%
  select(-protein, -cds) %>%
  mutate(across(c(ends_with("identity"), ends_with("coverage"), starts_with("shared")), ~ round(.x, 1)))
write_tsv(result, file.path(out_dir, "paralog_similarity.tsv"))
print(as.data.frame(result))

# Table S11 sheet, with the largest error of the simulated mixtures of 06_2 (simulation.tsv, made once)
simulation <- read_tsv(file.path(out_dir, "simulation.tsv"), show_col_types = FALSE) %>%
  group_by(gene_a, gene_b) %>%
  summarise(sim_error = round(max(abs(est_share_a - true_share_a)), 2), .groups = "drop")
note <- paste("Simulated mixtures: 100,000 read pairs (2 x 100 nt, 1% sequencing errors, 0.5% sequence differences",
              "from the reference) from the two transcripts, with gene a contributing 100%, 90%, 50%, 10% or 0%,",
              "were quantified with Salmon against all RefSeq transcripts. The column gives the largest difference",
              "between the estimated and the true share of gene a.")
table_tsv <- file.path(out_dir, "table_s9.tsv")
result %>%
  left_join(simulation, by = c("gene_a", "gene_b")) %>%
  { stopifnot(!anyNA(.$sim_error)); . } %>%
  transmute(`Gene a` = gene_a, `Gene b` = gene_b, `Transcript a` = transcript_a, `Transcript b` = transcript_b,
            `Protein identity (%)` = protein_identity, `CDS identity (%)` = cds_identity,
            `Shared 100-nt windows, gene a (%)` = shared_a, `Shared 100-nt windows, gene b (%)` = shared_b,
            `Simulated mixtures, largest error of the gene a share (percentage points)` = sim_error,
            `Highest atlas tissue, gene a` = atlas_peak_a, `Highest atlas tissue, gene b` = atlas_peak_b) %>%
  mutate(across(everything(), as.character)) %>%
  add_row(`Gene a` = c("", note)) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  write_tsv(table_tsv)
