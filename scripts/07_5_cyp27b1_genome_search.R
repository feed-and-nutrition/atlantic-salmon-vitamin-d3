# Search the salmon genomes for unannotated 1-alpha-hydroxylase (CYP27B1) genes (tracker item 38). Not part of
# run_all.R: it uses the genome assemblies and BLAST databases that 07_3 downloads to annotation/ (it downloads any that
# are missing) and takes about half an hour.
# 1. Human, zebrafish, spotted gar and northern pike CYP27B1 are searched with tBLASTn (E <= 1e-5) against the genome
#    assemblies of Atlantic salmon (Ssal_v3.1 of RS_2025_05, fSalSal1 of RS_2026_09) and, as a control with one
#    cyp27b1 gene (no Ss4R duplication), northern pike. Pike is not searched with its own protein.
# 2. Hits on the same sequence within 50 kb form a locus. CYP27B1 queries also find the other mitochondrial P450s, so
#    each locus is searched with BLASTx against a panel of the mitochondrial P450 families (CYP27B1, CYP27A1 and its
#    fish paralogs cyp27a2/a3, CYP24A1, CYP27C1, CYP11A1, CYP11B1) of human and non-salmonid fishes and seven human
#    microsomal P450s (Other P450), and assigned to the family with the highest summed bit score of non-overlapping
#    hits. The annotated genes overlapping each locus
#    are listed, so a CYP27B1 locus without an annotated cyp27b1 gene would be an unannotated copy.
# Needs tblastn, blastx, makeblastdb and blastdbcmd.
# Input:  NCBI FTP (genomes, feature tables, cached in annotation/), E-utilities (panel proteins)
# Output: vd3_genes/cyp27b1_search/loci.tsv

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

ann_dir <- file.path(proj_dir, "annotation")
blast_dir <- file.path(ann_dir, "blastdb")
cache_dir <- file.path(ann_dir, "cyp27b1_search")
out_dir <- file.path(vd3_dir, "cyp27b1_search")
for (d in c(blast_dir, cache_dir, out_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
threads <- 8

ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"
asm <- tribble(~label, ~species, ~dir, ~path,
  "Ssal_v3.1", "Atlantic salmon", "RS_2025_05", "905/237/065/GCF_905237065.1_Ssal_v3.1",
  "fSalSal1", "Atlantic salmon", "RS_2026_09", "965/601/325/GCF_965601325.2_fSalSal1",
  "fEsoLuc1", "Northern pike", "gc_search", "011/004/845/GCF_011004845.1_fEsoLuc1.pri")
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
           col_select = c(feature, class, chromosome, genomic_accession, start, end, name, symbol, GeneID)) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
}
read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(protein = str_extract(lines[head], "^>(\\S+)", group = 1),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}
read_blast <- function(out, fields) {
  read_tsv(I(paste(c(str_replace_all(fields, " ", "\t"), out), collapse = "\n")), show_col_types = FALSE)
}
write_fasta <- function(x, file) writeLines(paste0(">", x$protein, "\n", x$seq), file)

# --- Panel: NCBI Gene IDs, the longest RefSeq protein of each. query: used for tBLASTn
panel <- tribble(~label, ~species, ~class, ~gene_id, ~query,
  "Human_CYP27B1", "Human", "CYP27B1", "1594", TRUE,
  "Zebrafish_cyp27b1", "Zebrafish", "CYP27B1", "793309", TRUE,
  "Gar_cyp27b1", "Spotted gar", "CYP27B1", "102685924", TRUE,
  "Pike_cyp27b1", "Northern pike", "CYP27B1", "105017033", TRUE,
  "Medaka_cyp27b1", "Medaka", "CYP27B1", "101172325", FALSE,
  "Human_CYP27A1", "Human", "CYP27A1", "1593", FALSE,
  "Zebrafish_cyp27a1.1", "Zebrafish", "CYP27A1", "402831", FALSE,
  "Zebrafish_cyp27a2", "Zebrafish", "CYP27A1", "558239", FALSE,
  "Zebrafish_cyp27a3", "Zebrafish", "CYP27A1", "565876", FALSE,
  "Gar_cyp27a2", "Spotted gar", "CYP27A1", "102687850", FALSE,
  "Gar_cyp27a3", "Spotted gar", "CYP27A1", "102689682", FALSE,
  "Human_CYP24A1", "Human", "CYP24A1", "1591", FALSE,
  "Zebrafish_cyp24a1", "Zebrafish", "CYP24A1", "100004700", FALSE,
  "Human_CYP27C1", "Human", "CYP27C1", "339761", FALSE,
  "Zebrafish_cyp27c1", "Zebrafish", "CYP27C1", "558396", FALSE,
  "Human_CYP11A1", "Human", "CYP11A1", "1583", FALSE,
  "Human_CYP11B1", "Human", "CYP11B1", "1584", FALSE,
  "Human_CYP2R1", "Human", "Other P450", "120227", FALSE,
  "Human_CYP2J2", "Human", "Other P450", "1573", FALSE,
  "Human_CYP3A4", "Human", "Other P450", "1576", FALSE,
  "Human_CYP4V2", "Human", "Other P450", "285440", FALSE,
  "Human_CYP17A1", "Human", "Other P450", "1586", FALSE,
  "Human_CYP1A1", "Human", "Other P450", "1543", FALSE,
  "Human_CYP26A1", "Human", "Other P450", "1592", FALSE)
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
panel_fa <- file.path(cache_dir, "panel.faa")
panel_map <- file.path(cache_dir, "panel_proteins.tsv")
if (!file.exists(panel_fa) || !file.exists(panel_map)) {
  links <- jsonlite::fromJSON(paste0(eutils, "/elink.fcgi?dbfrom=gene&db=protein&linkname=gene_protein_refseq",
                                     "&retmode=json&", paste0("id=", panel$gene_id, collapse = "&")),
                              simplifyVector = FALSE)
  uids <- map_dfr(links$linksets, ~ tibble(gene_id = as.character(.x$ids[[1]]),
                                           uid = as.character(unlist(.x$linksetdbs[[1]]$links))))
  download.file(paste0(eutils, "/efetch.fcgi?db=protein&rettype=fasta&id=", paste(uids$uid, collapse = ",")),
                tmp <- tempfile(), quiet = TRUE)
  # efetch returns the proteins in the order requested
  prot <- bind_cols(uids, read_fasta(tmp)) %>%
    group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup()
  write_fasta(prot, panel_fa)
  write_tsv(select(prot, gene_id, protein), panel_map)
}
panel <- panel %>% inner_join(read_tsv(panel_map, col_types = "cc"), by = "gene_id") %>%
  inner_join(read_fasta(panel_fa), by = "protein")
stopifnot(nrow(panel) == 24)
panel_db <- file.path(blast_dir, "cyp27b1_panel")
write_fasta(panel, fa <- tempfile(fileext = ".faa"))
stopifnot(system2("makeblastdb", c("-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", shQuote(panel_db)),
                  stdout = FALSE) == 0)

genome_db <- function(label) {
  db <- file.path(blast_dir, paste0(label, "_genome"))
  if (!file.exists(paste0(db, ".nal")) && !file.exists(paste0(db, ".nin"))) {
    cmd <- sprintf("gzip -dc %s | makeblastdb -in - -dbtype nucl -parse_seqids -title %s -out %s",
                   shQuote(get_file(label, "_genomic.fna.gz")), label, shQuote(db))
    stopifnot(system(cmd, ignore.stdout = TRUE) == 0)
  }
  db
}
tblastn_fields <- "qseqid saccver pident length qstart qend sstart send evalue bitscore"
blastx_fields <- "saccver pident length qstart qend evalue bitscore"
# Summed bit score of HSPs that do not overlap on the locus, best first
summed_bits <- function(qstart, qend, bitscore) {
  o <- order(-bitscore)
  used <- matrix(numeric(0), ncol = 2)
  total <- 0
  for (i in o) {
    a <- min(qstart[i], qend[i]); b <- max(qstart[i], qend[i])
    if (!any(used[, 1] <= b & used[, 2] >= a)) {
      used <- rbind(used, c(a, b))
      total <- total + bitscore[i]
    }
  }
  total
}
search_genome <- function(label) {
  species <- asm$species[asm$label == label]
  db <- genome_db(label)
  # Raw tBLASTn hits are cached
  hits_file <- file.path(cache_dir, paste0(label, "_tblastn.tsv"))
  if (!file.exists(hits_file)) {
    write_fasta(filter(panel, query, species != !!species), fa <- tempfile(fileext = ".faa"))
    writeLines(system2("tblastn", c("-query", fa, "-db", shQuote(db), "-evalue", "1e-5", "-num_threads", threads,
                                    "-outfmt", shQuote(paste("6", tblastn_fields))), stdout = TRUE), hits_file)
  }
  hits <- read_blast(readLines(hits_file), tblastn_fields) %>%
    mutate(sacc = saccver, from = pmin(sstart, send), to = pmax(sstart, send))
  # Loci: hits on the same sequence within 50 kb
  loci <- hits %>% arrange(sacc, from) %>% group_by(sacc) %>%
    mutate(locus = cumsum(from > lag(cummax(to), default = -Inf) + 5e4)) %>%
    group_by(sacc, locus) %>%
    summarise(start = min(from), end = max(to), n_hits = n(), queries = n_distinct(qseqid),
              best_query = panel$label[match(qseqid[which.max(bitscore)], panel$protein)],
              best_query_bits = max(bitscore), best_query_identity = pident[which.max(bitscore)],
              .groups = "drop")
  other <- filter(panel, species != !!species)
  classified <- pmap_dfr(loci, function(sacc, start, end, ...) {
    region <- tempfile(fileext = ".fna")
    system2("blastdbcmd", c("-db", shQuote(db), "-entry", sacc, "-range", paste0(max(1, start - 2000), "-", end + 2000),
                            "-out", region))
    x <- read_blast(system2("blastx", c("-query", region, "-db", shQuote(panel_db), "-evalue", "1e-5",
                                        "-outfmt", shQuote(paste("6", blastx_fields))), stdout = TRUE), blastx_fields) %>%
      inner_join(select(other, saccver = protein, label, class), by = "saccver") %>%
      group_by(label, class) %>%
      summarise(bits = summed_bits(qstart, qend, bitscore), identity = max(pident), .groups = "drop") %>%
      arrange(desc(bits))
    best <- x %>% group_by(class) %>% slice_max(bits, n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(desc(bits))
    tibble(class = best$class[1], best_panel = best$label[1], class_bits = best$bits[1],
           cyp27b1_bits = sum(best$bits[best$class == "CYP27B1"]),
           second_class = best$class[2], second_bits = best$bits[2])
  })
  ft <- read_features(label)
  chr_of <- distinct(ft, sacc = genomic_accession, chr = chromosome)
  genes <- ft %>% filter(feature == "gene") %>%
    transmute(sacc = genomic_accession, g_start = start, g_end = end,
              gene = paste0(coalesce(symbol, GeneID), " (", GeneID, ", ", class, ")"))
  overlap <- loci %>% select(sacc, locus, start, end) %>%
    inner_join(genes, by = "sacc", relationship = "many-to-many") %>%
    filter(g_start <= end, g_end >= start) %>%
    group_by(sacc, locus) %>% summarise(genes = paste(unique(gene), collapse = "; "), .groups = "drop")
  bind_cols(loci, classified) %>%
    left_join(overlap, by = c("sacc", "locus")) %>%
    left_join(chr_of, by = "sacc") %>%
    mutate(assembly = label, species = species, .before = 1) %>%
    select(assembly, species, chr, sacc, start, end, n_hits, queries, best_query, best_query_bits,
           best_query_identity, class, best_panel, class_bits, cyp27b1_bits, second_class, second_bits, genes)
}
loci <- map_dfr(c("fEsoLuc1", "Ssal_v3.1", "fSalSal1"), search_genome)
write_tsv(loci, file.path(out_dir, "loci.tsv"))
print(count(loci, assembly, class))
print(as.data.frame(select(filter(loci, class == "CYP27B1"), assembly, chr, start, end, best_query_bits, class_bits,
                           second_class, second_bits, genes)))
