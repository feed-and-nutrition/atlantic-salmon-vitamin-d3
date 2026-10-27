# Search for GC (vitamin D-binding protein) in the fish annotations and genomes (tracker item 36). Not part of
# run_all.R: it downloads about 3.5 GB (genome assemblies, feature tables and proteomes) from the NCBI FTP site, needs
# about 15 GB of disk for the BLAST databases and takes about an hour.
# The panel holds GC and albumin-family proteins (Table S4, plus the argentine protein annotated as gc, an albumin).
# A protein is assigned to the class (GC or albumin) of the panel protein with the highest bit score, leaving out the
# panel proteins of its own species.
# 1. Annotated genes: all ray-finned fish genes in NCBI Gene with the symbol gc or the name vitamin D-binding protein,
#    classified by BLASTp of their longest RefSeq protein, with their taxonomic group.
# 2. Genome loci: human, spotted gar and zebrafish GC and human, gar, pike and salmon albumins are searched with
#    tBLASTn (E <= 1e-5) against the genome assemblies of Atlantic salmon (Ssal_v3.1 of RS_2025_05, fSalSal1 of
#    RS_2026_09), northern pike (no Ss4R duplication), the greater argentine (Argentina silus), medaka (a
#    percomorph) and, as a positive control with a gc gene, Atlantic herring. A species is not searched with its own
#    proteins. Hits on the same sequence within 50 kb form a locus, which is searched with BLASTx against the panel
#    and classified by the summed bit score of non-overlapping hits.
# 3. Neighbourhood: the ten protein-coding genes on each side of gc in spotted gar, which keeps the gene order of
#    human 4q13.3 (SLC4A4, GC, NPFFR2, ADAMTS3), are searched with BLASTp (E <= 1e-10) against the proteomes of
#    zebrafish, herring, the argentine, pike, salmon (fSalSal1) and medaka. Genes within 80% of the best bit score
#    are kept (up to three, to catch Ss4R copies), with their position and their rank among the protein-coding genes
#    of their sequence, so that adjacent hits show the conserved order.
# Downloads are cached in annotation/ (delete a file to download it again). Needs tblastn, blastx, blastp,
# makeblastdb and blastdbcmd.
# Table S3 holds the BLASTp search of the Methods (gc_blast/blast_results.tsv) and the three searches.
# Input:  NCBI FTP (genomes, feature tables, proteomes), E-utilities (panel proteins, gc genes),
#         gc_blast/blast_results.tsv
# Output: vd3_genes/gc_search/{annotated_gc,loci,neighbourhood}.tsv, manuscript/v6/Table S3 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

ann_dir <- file.path(proj_dir, "annotation")
blast_dir <- file.path(ann_dir, "blastdb")
out_dir <- file.path(vd3_dir, "gc_search")
dir.create(blast_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
threads <- 8

ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"
asm <- tribble(~label, ~species, ~dir, ~path,
  "Ssal_v3.1", "Atlantic salmon", "RS_2025_05", "905/237/065/GCF_905237065.1_Ssal_v3.1",
  "fSalSal1", "Atlantic salmon", "RS_2026_09", "965/601/325/GCF_965601325.2_fSalSal1",
  "fEsoLuc1", "Northern pike", "gc_search", "011/004/845/GCF_011004845.1_fEsoLuc1.pri",
  "fArgSil1", "Greater argentine", "gc_search", "951/799/395/GCF_951799395.1_fArgSil1.1",
  "Olat3.1", "Medaka", "gc_search", "053/564/925/GCF_053564925.1_Olatipes_Hd-rR_3.1",
  "Ch_v2.0.2", "Atlantic herring", "gc_search", "900/700/415/GCF_900700415.2_Ch_v2.0.2",
  "fLepOcu1", "Spotted gar", "gc_search", "040/954/835/GCF_040954835.1_fLepOcu1.hap2",
  "GRCz13ab", "Zebrafish", "gc_search", "052/040/795/GCF_052040795.2_GRCz13ab")
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
           col_select = c(feature, class, chromosome, genomic_accession, start, end, strand,
                          product_accession, name, symbol, GeneID)) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
}
read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(protein = str_extract(lines[head], "^>(\\S+)", group = 1),
         organism = str_extract(lines[head], "\\[([^]]+)\\]$", group = 1),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}
read_blast <- function(out, fields) {
  read_tsv(I(paste(c(str_replace_all(fields, " ", "\t"), out), collapse = "\n")), show_col_types = FALSE)
}
write_fasta <- function(x, file) writeLines(paste0(">", x$protein, "\n", x$seq), file)

# --- Panel. query: used for tBLASTn
panel <- tribble(~label, ~species, ~class, ~protein, ~query,
  "Human_GC", "Human", "GC", "NP_000574.2", TRUE,
  "Chicken_GC", "Chicken", "GC", "NP_990213.1", FALSE,
  "Gar_gc", "Spotted gar", "GC", "XP_069043117.1", TRUE,
  "Zebrafish_gc", "Zebrafish", "GC", "NP_001002568.3", TRUE,
  "Human_ALB", "Human", "Albumin", "NP_000468.1", TRUE,
  "Human_AFP", "Human", "Albumin", "NP_001125.1", FALSE,
  "Human_AFM", "Human", "Albumin", "NP_001124.1", FALSE,
  "Chicken_ALB", "Chicken", "Albumin", "NP_990592.2", FALSE,
  "Gar_alb", "Spotted gar", "Albumin", "XP_015195668.1", TRUE,
  "Pike_alb.1", "Northern pike", "Albumin", "XP_012992300.1", FALSE,
  "Pike_alb.2", "Northern pike", "Albumin", "XP_034148109.1", TRUE,
  "Salmon_alb1", "Atlantic salmon", "Albumin", "NP_001117137.1", TRUE,
  "Salmon_alb2", "Atlantic salmon", "Albumin", "XP_013998969.2", FALSE,
  "Salmon_LOC100136922", "Atlantic salmon", "Albumin", "NP_001117164.1", FALSE,
  "Argentine_alb", "Greater argentine", "Albumin", "XP_083529340.1", FALSE)
panel_fa <- file.path(ann_dir, "gc_search", "panel.faa")
if (!file.exists(panel_fa)) {
  download.file(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=",
                       paste(panel$protein, collapse = ",")), panel_fa, quiet = TRUE)
}
panel <- inner_join(panel, read_fasta(panel_fa), by = "protein")
stopifnot(nrow(panel) == 15)
panel_db <- file.path(blast_dir, "gc_panel")
write_fasta(panel, fa <- tempfile(fileext = ".faa"))
stopifnot(system2("makeblastdb", c("-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", shQuote(panel_db)),
                  stdout = FALSE) == 0)

efetch_fasta <- function(ids, dest) {
  if (!file.exists(dest)) {
    for (chunk in split(ids, ceiling(seq_along(ids) / 150))) {
      download.file(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=",
                           paste(chunk, collapse = ",")), tmp <- tempfile(), quiet = TRUE)
      cat(readLines(tmp), file = dest, sep = "\n", append = TRUE)
      Sys.sleep(0.4)
    }
  }
  read_fasta(dest)
}
blastp_fields <- "qseqid saccver pident length qlen evalue bitscore"
# Best panel hit of each class for proteins of one species (scientific name)
classify <- function(proteins, organism) {
  write_fasta(proteins, fa <- tempfile(fileext = ".faa"))
  read_blast(system2("blastp", c("-query", fa, "-db", shQuote(panel_db), "-evalue", "1e-5", "-max_target_seqs", 50,
                                 "-outfmt", shQuote(paste("6", blastp_fields))), stdout = TRUE), blastp_fields) %>%
    inner_join(select(filter(panel, organism != !!organism), saccver = protein, class), by = "saccver") %>%
    group_by(protein = qseqid) %>%
    summarise(bits_GC = max(c(0, bitscore[class == "GC"])),
              identity_GC = pident[class == "GC"][which.max(bitscore[class == "GC"])][1],
              bits_Albumin = max(c(0, bitscore[class == "Albumin"])),
              identity_Albumin = pident[class == "Albumin"][which.max(bitscore[class == "Albumin"])][1]) %>%
    mutate(class = if_else(bits_GC > bits_Albumin, "GC", "Albumin"))
}

# --- 1. Annotated gc genes of ray-finned fishes
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
gene_file <- file.path(ann_dir, "gc_search", "annotated_gc_genes.json")
if (!file.exists(gene_file)) {
  ids <- jsonlite::fromJSON(paste0(eutils, "/esearch.fcgi?db=gene&retmax=500&retmode=json&term=",
                                   URLencode(paste("(gc[sym] OR \"vitamin D-binding protein\"[gene full name])",
                                                   "AND Actinopterygii[orgn] AND alive[prop]"), reserved = TRUE)))
  ids <- ids$esearchresult$idlist
  links <- jsonlite::fromJSON(paste0(eutils, "/elink.fcgi?dbfrom=gene&db=protein&linkname=gene_protein_refseq",
                                     "&retmode=json&", paste0("id=", ids, collapse = "&")), simplifyVector = FALSE)
  summ <- jsonlite::fromJSON(paste0(eutils, "/esummary.fcgi?db=gene&retmode=json&id=", paste(ids, collapse = ",")),
                             simplifyVector = FALSE)$result
  genes <- map(ids, function(id) {
    g <- summ[[id]]
    ls <- keep(links$linksets, ~ .x$ids[[1]] == id)[[1]]
    list(gene_id = id, symbol = g$name, description = g$description, species = g$organism$scientificname,
         taxid = g$organism$taxid, chr = g$chromosome,
         protein_uids = if (length(ls$linksetdbs)) unlist(ls$linksetdbs[[1]]$links) else character(0))
  })
  # Taxonomic group from the lineage
  taxids <- unique(map_chr(genes, ~ as.character(.x$taxid)))
  lineage <- xml2::read_xml(paste0(eutils, "/efetch.fcgi?db=taxonomy&id=", paste(taxids, collapse = ",")))
  taxa <- xml2::xml_find_all(lineage, "/TaxaSet/Taxon")
  groups <- set_names(map_chr(taxa, function(t) {
    l <- xml2::xml_text(xml2::xml_find_all(t, "./LineageEx/Taxon/ScientificName"))
    c(intersect(c("Euteleosteomorpha", "Otomorpha", "Osteoglossomorpha", "Elopomorpha"), l), "Non-teleost")[1]
  }), map_chr(taxa, ~ xml2::xml_text(xml2::xml_find_first(.x, "./TaxId"))))
  genes <- map(genes, ~ c(.x, group = unname(groups[as.character(.x$taxid)])))
  jsonlite::write_json(genes, gene_file, auto_unbox = TRUE)
}
genes_gc <- jsonlite::fromJSON(gene_file, simplifyVector = FALSE)
prot_gc <- efetch_fasta(unlist(map(genes_gc, "protein_uids")), file.path(ann_dir, "gc_search", "annotated_gc.faa"))
gene_of_protein <- map_dfr(genes_gc, function(g) {
  tibble(gene_id = g$gene_id, uid = as.character(unlist(g$protein_uids)))
})
# efetch returns the proteins in the order requested, so the uids map to the accessions in order (checked by species)
stopifnot(nrow(prot_gc) == nrow(gene_of_protein))
annotated <- bind_cols(gene_of_protein, prot_gc) %>%
  group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup() %>%
  left_join(map_dfr(genes_gc, ~ tibble(gene_id = .x$gene_id, symbol = .x$symbol, description = .x$description,
                                       species = .x$species, group = .x$group, chr = .x$chr)), by = "gene_id")
stopifnot(all(annotated$organism == annotated$species))
annotated <- annotated %>%
  group_by(species) %>% group_modify(~ left_join(.x, classify(.x, .y$species), by = "protein")) %>% ungroup() %>%
  transmute(group = factor(group, levels = c("Non-teleost", "Osteoglossomorpha", "Elopomorpha", "Otomorpha",
                                             "Euteleosteomorpha")),
            species, gene_id, symbol, description, chr, protein, length = nchar(seq), class,
            gc_bits = bits_GC, gc_identity = identity_GC, albumin_bits = bits_Albumin,
            albumin_identity = identity_Albumin) %>%
  arrange(group, class, species)
write_tsv(annotated, file.path(out_dir, "annotated_gc.tsv"))
print(count(annotated, group, class))

# --- 2. Genome loci
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
  hits_file <- file.path(ann_dir, "gc_search", paste0(label, "_tblastn.tsv"))
  if (!file.exists(hits_file)) {
    write_fasta(filter(panel, query, species != !!species), fa <- tempfile(fileext = ".faa"))
    writeLines(system2("tblastn", c("-query", fa, "-db", shQuote(db), "-evalue", "1e-5", "-num_threads", threads,
                                    "-outfmt", shQuote(paste("6", tblastn_fields))), stdout = TRUE), hits_file)
  }
  hits <- read_blast(readLines(hits_file), tblastn_fields) %>%
    mutate(sacc = saccver, from = pmin(sstart, send), to = pmax(sstart, send))
  if (!nrow(hits)) return(tibble(assembly = label, species = species, class = "No hit"))
  # Loci: hits on the same sequence within 50 kb
  loci <- hits %>% arrange(sacc, from) %>% group_by(sacc) %>%
    mutate(locus = cumsum(from > lag(cummax(to), default = -Inf) + 5e4)) %>%
    group_by(sacc, locus) %>%
    summarise(start = min(from), end = max(to), n_hits = n(), queries = n_distinct(qseqid),
              gc_query_hits = sum(qseqid %in% panel$protein[panel$class == "GC"]),
              best_query = qseqid[which.max(bitscore)], best_query_bits = max(bitscore),
              best_query_evalue = min(evalue), .groups = "drop")
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
    best_of <- function(cl) { y <- filter(x, class == cl); if (nrow(y)) y$bits[1] else 0 }
    tibble(class = if (nrow(x)) x$class[1] else NA, best_panel = if (nrow(x)) x$label[1] else NA,
           gc_bits = best_of("GC"), albumin_bits = best_of("Albumin"))
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
    select(assembly, species, chr, sacc, start, end, n_hits, queries, gc_query_hits, best_query, best_query_bits,
           best_query_evalue, class, best_panel, gc_bits, albumin_bits, genes)
}
loci <- map_dfr(c("Ch_v2.0.2", "fArgSil1", "fEsoLuc1", "Ssal_v3.1", "fSalSal1", "Olat3.1"), search_genome)
write_tsv(loci, file.path(out_dir, "loci.tsv"))
print(as.data.frame(select(loci, assembly, chr, start, end, best_query, best_query_bits, class, best_panel,
                           gc_bits, albumin_bits, genes)))

# --- 3. Neighbourhood of gc in spotted gar
coding_genes <- function(ft) {
  ft %>% filter(feature == "gene", class == "protein_coding") %>%
    transmute(gene_id = GeneID, symbol, chr = chromosome, seqid = genomic_accession, start, end) %>%
    arrange(seqid, start) %>% group_by(seqid) %>% mutate(rank = row_number()) %>% ungroup()
}
longest_protein <- function(label, ft) {
  cds <- ft %>% filter(feature == "CDS", !is.na(product_accession)) %>%
    distinct(protein = product_accession, gene_id = GeneID, product = name)
  read_fasta(get_file(label, "_protein.faa.gz")) %>% inner_join(cds, by = "protein") %>%
    group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup()
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

ft_gar <- read_features("fLepOcu1")
gar <- coding_genes(ft_gar)
gc_gar <- gar %>% filter(symbol == "gc")
flank <- gar %>% filter(seqid == gc_gar$seqid, abs(rank - gc_gar$rank) <= 10) %>%
  mutate(offset = rank - gc_gar$rank) %>%
  inner_join(select(longest_protein("fLepOcu1", ft_gar), gene_id, protein, seq, product), by = "gene_id")
write_fasta(flank, flank_fa <- tempfile(fileext = ".faa"))
neighbourhood <- map_dfr(c("GRCz13ab", "Ch_v2.0.2", "fArgSil1", "fEsoLuc1", "fSalSal1", "Olat3.1"), function(label) {
  ft <- read_features(label)
  cds <- ft %>% filter(feature == "CDS", !is.na(product_accession)) %>%
    distinct(saccver = product_accession, hit_gene_id = GeneID)
  read_blast(system2("blastp", c("-query", flank_fa, "-db", shQuote(proteome_db(label)), "-evalue", "1e-10",
                                 "-max_target_seqs", 50, "-num_threads", threads,
                                 "-outfmt", shQuote(paste("6", blastp_fields))), stdout = TRUE), blastp_fields) %>%
    inner_join(cds, by = "saccver") %>%
    group_by(qseqid, hit_gene_id) %>% summarise(bitscore = max(bitscore), identity = max(pident), .groups = "drop") %>%
    group_by(qseqid) %>% filter(bitscore >= 0.8 * max(bitscore)) %>%
    slice_max(bitscore, n = 3, with_ties = FALSE) %>% ungroup() %>%
    inner_join(rename_with(coding_genes(ft), ~ paste0("hit_", .x)), by = "hit_gene_id") %>%
    mutate(assembly = label, species = asm$species[asm$label == label])
}) %>%
  right_join(select(flank, offset, qseqid = protein, gar_gene_id = gene_id, gar_symbol = symbol, gar_product = product),
             by = "qseqid") %>%
  select(offset, gar_symbol, gar_gene_id, gar_product, assembly, species, hit_symbol, hit_gene_id, hit_chr,
         hit_seqid, hit_start, hit_rank, identity, bitscore) %>%
  arrange(offset, assembly, desc(bitscore))
write_tsv(neighbourhood, file.path(out_dir, "neighbourhood.tsv"))
print(as.data.frame(select(neighbourhood, offset, gar_symbol, assembly, hit_symbol, hit_chr, hit_start, hit_rank,
                           identity)))

# --- Table S3: the BLASTp search of the Methods (gc_blast/blast_results.tsv) and the three searches above
sheet_dir <- tempfile()
dir.create(sheet_dir)
sheet <- function(x, name, note) {
  file <- file.path(sheet_dir, paste0(name, ".tsv"))
  write_tsv(x, file, na = "")
  if (!is.null(note)) cat("\n", note, "\n", file = file, append = TRUE, sep = "")
  file
}
label_of <- set_names(str_replace(panel$label, "_", " "), panel$protein)
gene_label <- function(x) str_replace_all(x, "\\((\\d+), [a-z_]+\\)", "(Gene ID \\1)")
invisible(file.copy(file.path(proj_dir, "gc_blast", "blast_results.tsv"), file.path(sheet_dir, "BLAST results.tsv")))
sheets <- c(
  file.path(sheet_dir, "BLAST results.tsv"),
  sheet(annotated %>%
          transmute(Group = group, Species = species, `Gene ID` = gene_id, Symbol = symbol, Description = description,
                    Chromosome = chr, Protein = protein, `Length (aa)` = length, `Best match` = class,
                    `GC bit score` = gc_bits, `GC identity (%)` = round(gc_identity, 1),
                    `Albumin bit score` = albumin_bits, `Albumin identity (%)` = round(albumin_identity, 1)),
        "Annotated gc genes",
        paste("Ray-finned fish genes in NCBI Gene with the symbol gc or the name vitamin D-binding protein. The longest",
              "RefSeq protein of each gene was compared by BLASTp with the GC and albumin-family proteins of the gc",
              "tree (Supplementary Table S4) and the argentine protein annotated as gc, leaving out proteins of the",
              "same species. Best match: the class (GC or albumin family) with the higher bit score.")),
  sheet(loci %>%
          left_join(transmute(asm, assembly = label,
                              name = str_replace(basename(path), "^(GCF_[0-9.]+)_(.+)$", "\\2 (\\1)")), by = "assembly") %>%
          transmute(Species = species, Assembly = name, Chromosome = chr, Sequence = sacc, Start = start, End = end,
                    `tBLASTn hits` = n_hits, `Hits of the GC queries` = gc_query_hits, `Best query` = unname(label_of[best_query]),
                    `Best query bit score` = best_query_bits, `Best match` = class,
                    `Best BLASTx match` = str_replace(best_panel, "_", " "), `GC bit score` = round(gc_bits),
                    `Albumin bit score` = round(albumin_bits), `Annotated gene` = gene_label(genes)),
        "Genome search",
        paste("tBLASTn (E-value of 0.00001 or lower) of GC (human, spotted gar, zebrafish) and albumins (human, spotted",
              "gar, northern pike, Atlantic salmon) against each genome, without the proteins of the same species.",
              "Hits within 50 kb form a locus. Each locus was compared by BLASTx with the GC and albumin-family",
              "proteins and assigned to the class of the best match (summed bit scores of non-overlapping hits).",
              "Atlantic herring, which has a gc gene, is the positive control. The GC queries had no hit in the other",
              "genomes, and medaka had no hit at all.")),
  sheet(neighbourhood %>%
          filter(!is.na(assembly)) %>%
          mutate(hit = paste0(coalesce(hit_symbol, hit_gene_id), " (", hit_chr, ": ", hit_rank, ")"),
                 species = factor(species, levels = c("Zebrafish", "Atlantic herring", "Greater argentine",
                                                      "Northern pike", "Atlantic salmon", "Medaka"))) %>%
          group_by(offset, gar_symbol, gar_gene_id, species) %>%
          summarise(hit = paste(hit, collapse = "; "), .groups = "drop") %>%
          arrange(species) %>%
          pivot_wider(names_from = species, values_from = hit) %>%
          right_join(distinct(neighbourhood, offset, gar_symbol, gar_gene_id), by = c("offset", "gar_symbol", "gar_gene_id")) %>%
          arrange(offset) %>%
          rename(`Position from gc` = offset, `Spotted gar gene` = gar_symbol, `Gene ID` = gar_gene_id),
        "gc neighbourhood",
        paste("The ten protein-coding genes on each side of gc in spotted gar, and their best BLASTp matches (E-value of",
              "1e-10 or lower, within 80% of the best bit score, up to three genes) in each species. In parentheses:",
              "chromosome and rank among the protein-coding genes of that chromosome, so that consecutive ranks are",
              "neighbouring genes. The gc row lists the closest genes of each species, which are albumins in the species",
              "without gc, including the argentine gene annotated as gc.")))
table_file <- file.path(fig_dir, "Table S3 Search for GC in Atlantic salmon and other fishes.xlsx")
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(sheets),
                              shQuote(table_file))) == 0)
message("Wrote ", table_file)
