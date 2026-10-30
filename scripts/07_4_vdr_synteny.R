# Synteny of the vdr genes (tracker item 40): are salmon vdra.a/vdra.b and vdrb.a/vdrb.b the teleost vdra and vdrb?
# Not part of run_all.R: it reads the feature tables, proteomes and BLAST databases that 07_3 downloads to annotation/
# (it downloads any that are missing) and takes a few minutes.
# 1. Teleost loci against zebrafish: for each vdr gene of salmon (fSalSal1), pike (no Ss4R duplication) and medaka,
#    the longest proteins of the ten protein-coding genes on each side are searched with BLASTp (E <= 1e-10) against
#    the zebrafish proteome. The best match of each is counted when it lies on the chromosome of zebrafish vdra
#    (chr 23) or vdrb (chr 6), and when it lies within 100 genes of them (assignment by the latter).
# 2. Spotted gar (no teleost genome duplication): the ten genes on each side of the single gar vdr are searched
#    against each teleost proteome. Matches within 80% of the best bit score (up to four, to catch the TGD and Ss4R
#    copies) are counted when they lie within 100 genes of a teleost vdr. Gar flanking genes near both teleost copies
#    show that vdra and vdrb come from the teleost genome duplication.
# Needs blastp and makeblastdb.
# Input:  NCBI FTP (feature tables and proteomes, cached in annotation/)
# Output: vd3_genes/vdr_synteny/{zebrafish_matches,gar_matches,summary}.tsv (summary is a sheet of Table S5, 07_2)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

ann_dir <- file.path(proj_dir, "annotation")
blast_dir <- file.path(ann_dir, "blastdb")
out_dir <- file.path(vd3_dir, "vdr_synteny")
dir.create(blast_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
threads <- 8
n_flank <- 10
window <- 100

ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"
asm <- tribble(~label, ~species, ~dir, ~path,
  "fSalSal1", "Atlantic salmon", "RS_2026_09", "965/601/325/GCF_965601325.2_fSalSal1",
  "fEsoLuc1", "Northern pike", "gc_search", "011/004/845/GCF_011004845.1_fEsoLuc1.pri",
  "Olat3.1", "Medaka", "gc_search", "053/564/925/GCF_053564925.1_Olatipes_Hd-rR_3.1",
  "GRCz13ab", "Zebrafish", "gc_search", "052/040/795/GCF_052040795.2_GRCz13ab",
  "fLepOcu1", "Spotted gar", "gc_search", "040/954/835/GCF_040954835.1_fLepOcu1.hap2")
# The vdr genes (Gene IDs as in msa/vdr_extended/vdr_extended_accessions.tsv) and their names in the manuscript
vdr <- tribble(~label, ~gene_id, ~gene,
  "fSalSal1", "106566896", "vdra.a", "fSalSal1", "100136412", "vdra.b",
  "fSalSal1", "106565233", "vdrb.a", "fSalSal1", "100136371", "vdrb.b",
  "fEsoLuc1", "105022772", "LOC105022772", "fEsoLuc1", "105024800", "LOC105024800",
  "Olat3.1", "100169845", "vdr", "Olat3.1", "100169844", "vdralpha",
  "GRCz13ab", "30076", "vdra", "GRCz13ab", "564511", "vdrb",
  "fLepOcu1", "102691130", "LOC102691130")

get_file <- function(label, suffix) {
  a <- asm[asm$label == label, ]
  dest <- file.path(ann_dir, a$dir, paste0(basename(a$path), suffix))
  if (!file.exists(dest)) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    options(timeout = max(3600, getOption("timeout")))
    download.file(file.path(ftp, a$path, paste0(basename(a$path), suffix)), dest, quiet = TRUE, mode = "wb")
  }
  dest
}
read_features <- function(label) {
  file <- get_file(label, "_feature_table.txt.gz")
  cols <- str_split(str_remove(readLines(file, n = 1), "^# "), "\t")[[1]]
  read_tsv(file, skip = 1, col_names = cols, col_types = cols(.default = "c"),
           col_select = c(feature, class, chromosome, genomic_accession, start, end, product_accession, name, symbol,
                          GeneID)) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
}
read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(protein = str_extract(lines[head], "^>(\\S+)", group = 1),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}
write_fasta <- function(x, file) writeLines(paste0(">", x$protein, "\n", x$seq), file)
blastp_fields <- "qseqid saccver pident length qlen evalue bitscore"
read_blast <- function(out) {
  read_tsv(I(paste(c(str_replace_all(blastp_fields, " ", "\t"), out), collapse = "\n")), show_col_types = FALSE)
}
proteome_db <- function(label) {
  db <- file.path(blast_dir, paste0(label, "_proteome"))
  if (!file.exists(paste0(db, ".pin"))) {
    cmd <- sprintf("gzip -dc %s | makeblastdb -in - -dbtype prot -parse_seqids -title %s -out %s",
                   shQuote(get_file(label, "_protein.faa.gz")), label, shQuote(db))
    stopifnot(system(cmd, ignore.stdout = TRUE) == 0)
  }
  db
}

ft <- map(set_names(asm$label), read_features)
coding <- map(ft, ~ .x %>% filter(feature == "gene", class == "protein_coding") %>%
  transmute(gene_id = GeneID, symbol, chr = chromosome, seqid = genomic_accession, start) %>%
  arrange(seqid, start) %>% group_by(seqid) %>% mutate(rank = row_number()) %>% ungroup())
cds <- map(ft, ~ .x %>% filter(feature == "CDS", !is.na(product_accession)) %>%
  distinct(saccver = product_accession, gene_id = GeneID, product = name))
longest <- function(label) {
  read_fasta(get_file(label, "_protein.faa.gz")) %>% inner_join(cds[[label]], by = c(protein = "saccver")) %>%
    group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup()
}
# Each vdr gene with its sequence and rank among the protein-coding genes of its sequence
vdr <- vdr %>%
  mutate(species = asm$species[match(label, asm$label)],
         map2_dfr(label, gene_id, ~ select(filter(coding[[.x]], gene_id == .y), chr, seqid, rank)))
stopifnot(!anyNA(vdr$rank))
flanking <- function(label, gene_id) {
  v <- vdr[vdr$label == label & vdr$gene_id == gene_id, ]
  genes <- coding[[label]] %>% filter(seqid == v$seqid, abs(rank - v$rank) %in% seq_len(n_flank)) %>%
    mutate(offset = rank - v$rank)
  inner_join(genes, select(longest(label), gene_id, protein, seq), by = "gene_id")
}
search <- function(query, label) {
  write_fasta(query, fa <- tempfile(fileext = ".faa"))
  read_blast(system2("blastp", c("-query", fa, "-db", shQuote(proteome_db(label)), "-evalue", "1e-10",
                                 "-max_target_seqs", 50, "-num_threads", threads,
                                 "-outfmt", shQuote(paste("6", blastp_fields))), stdout = TRUE)) %>%
    inner_join(rename(cds[[label]], hit_gene_id = gene_id), by = "saccver") %>%
    group_by(qseqid, hit_gene_id) %>% summarise(bitscore = max(bitscore), identity = max(pident), .groups = "drop") %>%
    inner_join(rename_with(coding[[label]], ~ paste0("hit_", .x)), by = "hit_gene_id")
}
# The vdr gene of `label` within `window` genes of a hit, if any
near_vdr <- function(label, seqid, rank) {
  v <- vdr[vdr$label == label, ]
  map2_chr(seqid, rank, function(s, r) {
    i <- which(v$seqid == s & abs(v$rank - r) <= window)
    if (length(i)) v$gene[i[1]] else NA_character_
  })
}

# --- 1. Teleost vdr loci against zebrafish
zf <- vdr[vdr$label == "GRCz13ab", ]
zebrafish_matches <- vdr %>% filter(label %in% c("fSalSal1", "fEsoLuc1", "Olat3.1")) %>%
  pmap_dfr(function(label, gene_id, gene, species, ...) {
    q <- flanking(label, gene_id)
    search(q, "GRCz13ab") %>% group_by(qseqid) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup() %>%
      right_join(select(q, qseqid = protein, offset, flank_gene_id = gene_id, flank_symbol = symbol), by = "qseqid") %>%
      mutate(species = species, gene = gene, zebrafish_chr_of = zf$gene[match(hit_chr, zf$chr)],
             near_zebrafish = near_vdr("GRCz13ab", hit_seqid, hit_rank))
  }) %>%
  select(species, gene, offset, flank_symbol, flank_gene_id, hit_symbol, hit_gene_id, hit_chr, hit_rank,
         identity, bitscore, zebrafish_chr_of, near_zebrafish) %>%
  arrange(species, gene, offset)
write_tsv(zebrafish_matches, file.path(out_dir, "zebrafish_matches.tsv"))

# --- 2. Gar vdr neighbourhood against the teleosts
gar_flank <- flanking("fLepOcu1", vdr$gene_id[vdr$label == "fLepOcu1"])
gar_matches <- map_dfr(c("GRCz13ab", "Olat3.1", "fEsoLuc1", "fSalSal1"), function(label) {
  search(gar_flank, label) %>% group_by(qseqid) %>% filter(bitscore >= 0.8 * max(bitscore)) %>%
    slice_max(bitscore, n = 4, with_ties = FALSE) %>% ungroup() %>%
    mutate(species = asm$species[asm$label == label], near_vdr = near_vdr(label, hit_seqid, hit_rank))
}) %>%
  right_join(select(gar_flank, qseqid = protein, offset, gar_gene_id = gene_id, gar_symbol = symbol), by = "qseqid") %>%
  select(species, offset, gar_symbol, gar_gene_id, hit_symbol, hit_gene_id, hit_chr, hit_rank, identity, bitscore,
         near_vdr) %>%
  arrange(species, offset, desc(bitscore))
write_tsv(gar_matches, file.path(out_dir, "gar_matches.tsv"))

# --- Summary: one row per teleost vdr gene
gar_counts <- gar_matches %>% filter(!is.na(near_vdr)) %>% distinct(species, near_vdr, gar_gene_id) %>%
  count(species, gene = near_vdr, name = "gar_flanks_near")
summary <- vdr %>% filter(label != "fLepOcu1") %>%
  select(species, gene, gene_id, chr) %>%
  left_join(zebrafish_matches %>% group_by(species, gene) %>%
              summarise(with_match = sum(!is.na(hit_gene_id)),
                        on_vdra_chr = sum(zebrafish_chr_of %in% "vdra"), on_vdrb_chr = sum(zebrafish_chr_of %in% "vdrb"),
                        near_vdra = sum(near_zebrafish %in% "vdra"), near_vdrb = sum(near_zebrafish %in% "vdrb"),
                        .groups = "drop"),
            by = c("species", "gene")) %>%
  left_join(gar_counts, by = c("species", "gene")) %>%
  mutate(gar_flanks_near = replace_na(gar_flanks_near, 0L),
         assigned = case_when(species == "Zebrafish" ~ gene,
                              near_vdra > near_vdrb ~ "vdra",
                              near_vdrb > near_vdra ~ "vdrb"))
write_tsv(summary, file.path(out_dir, "summary.tsv"))
print(as.data.frame(summary))
print(as.data.frame(count(filter(gar_matches, !is.na(near_vdr)), species, near_vdr)))
