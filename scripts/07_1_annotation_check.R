# Check of the VD3-related genes (RS_2025_05 on Ssal_v3.1) against the newer RefSeq annotation RS_2026_09
# (GCF_965601325.2, assembly fSalSal1). Not part of run_all.R: it downloads about 95 MB from the NCBI FTP site
# and takes about 6 minutes.
# 1. Gene map: location in both annotations and identity of the manuscript protein (Table S4) to the protein
#    of the same GeneID in RS_2026_09. in_old: the cited protein version is a protein of the gene in RS_2025_05.
# 2. Paralog search: human and zebrafish proteins (Table S4) against the RS_2026_09 proteome (blastp, E <= 1e-10).
#    A hit gene is a candidate when its longest protein has its best human hit (reciprocal) in a query gene or,
#    for GC, in the albumin family (ALB, AFP, AFM, GC).
# 3. Synteny: for each paralog pair, the best hit (other chromosome) of the 10 genes on each side of one copy.
#    Ss4R paralogs sit in homeologous blocks, so most flanking genes should hit the chromosome of the other copy.
#    For a single-copy gene the same test names the chromosome where the second copy would be expected.
# Downloads are cached in annotation/ (delete a file to download it again). Needs blastp and makeblastdb.
# Input:  vd3_genes/salmon_vd3_orthologs_v1.tsv, NCBI FTP (RS_2025_05, RS_2026_09, human GRCh38.p14), E-utilities
# Output: vd3_genes/annotation_check/{gene_map,paralog_search,synteny}.tsv

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

ann_dir <- file.path(proj_dir, "annotation")
out_dir <- file.path(vd3_dir, "annotation_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"
releases <- list(
  old = c(dir = "RS_2025_05", url = file.path(ftp, "905/237/065/GCF_905237065.1_Ssal_v3.1/GCF_905237065.1_Ssal_v3.1")),
  new = c(dir = "RS_2026_09", url = file.path(ftp, "965/601/325/GCF_965601325.2_fSalSal1/GCF_965601325.2_fSalSal1")),
  human = c(dir = "human_GRCh38.p14", url = file.path(ftp, "000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14")))
get_file <- function(rel, suffix) {
  dest <- file.path(ann_dir, releases[[rel]][["dir"]], paste0(basename(releases[[rel]][["url"]]), suffix))
  if (!file.exists(dest)) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    options(timeout = max(600, getOption("timeout")))
    download.file(paste0(releases[[rel]][["url"]], suffix), dest, quiet = TRUE, mode = "wb")
  }
  dest
}

# --- Annotation tables
read_features <- function(rel) {
  file <- get_file(rel, "_feature_table.txt.gz")
  cols <- str_split(str_remove(readLines(file, n = 1), "^# "), "\t")[[1]]
  read_tsv(file, skip = 1, col_names = cols, col_types = cols(.default = "c"),
           col_select = c(feature, class, chromosome, genomic_accession, start, end, strand,
                          product_accession, name, symbol, GeneID))
}
genes_of <- function(ft) {
  ft %>% filter(feature == "gene", class == "protein_coding") %>%
    transmute(gene_id = GeneID, symbol, product = name, chr = chromosome, seqid = genomic_accession,
              start = as.numeric(start), end = as.numeric(end))
}
# Protein accession -> GeneID
cds_of <- function(ft) ft %>% filter(feature == "CDS", !is.na(product_accession)) %>% distinct(protein = product_accession, gene_id = GeneID)

read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(protein = str_extract(lines[head], "^>(\\S+)", group = 1),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}

ft <- map(c(old = "old", new = "new", human = "human"), read_features)
genes <- map(ft, genes_of)
prot <- map(c(old = "old", new = "new", human = "human"),
            ~ read_fasta(get_file(.x, "_protein.faa.gz")) %>% inner_join(cds_of(ft[[.x]]), by = "protein"))
# Longest protein of each gene
longest <- map(prot, ~ .x %>% group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup())

blast_dir <- file.path(ann_dir, "blastdb")
dir.create(blast_dir, showWarnings = FALSE)
make_db <- function(rel) {
  db <- file.path(blast_dir, rel)
  if (!file.exists(paste0(db, ".pin")) && !file.exists(paste0(db, ".00.pin"))) {
    fa <- tempfile(fileext = ".faa")
    writeLines(paste0(">", prot[[rel]]$protein, "\n", prot[[rel]]$seq), fa)
    stopifnot(system2("makeblastdb", c("-in", fa, "-dbtype", "prot", "-out", shQuote(db)), stdout = FALSE) == 0)
  }
  db
}
blastp <- function(query, rel, evalue = 1e-10, max_target = 500) {
  fa <- tempfile(fileext = ".faa")
  writeLines(paste0(">", query$protein, "\n", query$seq), fa)
  out <- system2("blastp", c("-query", fa, "-db", shQuote(make_db(rel)), "-evalue", evalue,
                             "-max_target_seqs", max_target, "-num_threads", 4,
                             "-outfmt", shQuote("6 qseqid sseqid pident length qlen slen evalue bitscore")),
                 stdout = TRUE)
  read_tsv(I(paste(c("qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore", out), collapse = "\n")),
           col_types = "ccdiiidd")
}
# Best hit per query, with the subject's gene
best_hit <- function(hits, rel) {
  hits %>% inner_join(cds_of(ft[[rel]]), by = c(sseqid = "protein")) %>%
    group_by(qseqid) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup()
}

# --- 1. Gene map
vd3 <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
  transmute(name, gene_id, protein_ms = str_trim(str_extract(`Protein refseq`, "^[^,]+")))
# Third albumin, annotated only in RS_2026_09 (msa/gc/gc_accessions.tsv). cyp27a1.b was added to the list after
# this check found it (the first run listed 106586130 as a candidate)
vd3 <- bind_rows(vd3, tibble(name = "LOC100136922", gene_id = "100136922", protein_ms = NA))

# Manuscript proteins at the versions cited (until v4 round 8 alb2 cited XP_013998969.1, older than RS_2025_05)
ms_fa <- file.path(ann_dir, "manuscript_proteins.faa")
if (!file.exists(ms_fa)) {
  download.file(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=",
                       paste(na.omit(vd3$protein_ms), collapse = ",")), ms_fa, quiet = TRUE)
}
ms_seq <- read_fasta(ms_fa) %>% select(protein_ms = protein, seq_ms = seq)
stopifnot(setequal(ms_seq$protein_ms, na.omit(vd3$protein_ms)))
identity_to_gene <- function(protein_ms, seq_ms, gene_id) {
  if (is.na(seq_ms)) return(tibble(protein_new = NA_character_, identity = NA_real_, coverage = NA_real_))
  hits <- blastp(tibble(protein = protein_ms, seq = seq_ms), "new", evalue = 1e-5, max_target = 5000) %>%
    filter(sseqid %in% prot$new$protein[prot$new$gene_id == gene_id]) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE)
  if (!nrow(hits)) return(tibble(protein_new = NA_character_, identity = 0, coverage = 0))
  tibble(protein_new = hits$sseqid, identity = hits$pident, coverage = round(100 * hits$length / hits$qlen, 1))
}
loc <- function(g) mutate(g, location = paste0(if_else(is.na(chr), seqid, chr), ":", start, "-", end)) %>%
  select(gene_id, symbol, product, chr, location)
gene_map <- vd3 %>%
  left_join(ms_seq, by = "protein_ms") %>%
  left_join(rename_with(loc(genes$old), ~ paste0(.x, "_old"), -gene_id), by = "gene_id") %>%
  left_join(rename_with(loc(genes$new), ~ paste0(.x, "_new"), -gene_id), by = "gene_id") %>%
  mutate(in_old = map2_lgl(seq_ms, gene_id, ~ .x %in% prot$old$seq[prot$old$gene_id == .y]),
         same_protein = protein_ms %in% prot$new$protein & seq_ms %in% prot$new$seq,
         map_dfr(pmap(list(protein_ms, seq_ms, gene_id), identity_to_gene), identity)) %>%
  select(-seq_ms)
write_tsv(gene_map, file.path(out_dir, "gene_map.tsv"))
print(as.data.frame(select(gene_map, name, gene_id, chr_old, chr_new, symbol_new, in_old, same_protein, identity, coverage)))

# --- 2. Paralog search
queries <- tribble(~query_gene, ~species, ~protein,
  "DHCR7", "human", "NP_001351.2", "DHCR7", "zebrafish", "NP_958487.2",
  "GC", "human", "NP_000574.2", "GC", "zebrafish", "NP_001002568.3",
  "CYP2R1", "human", "NP_078790.2", "CYP2R1", "zebrafish", "NP_001373291.1",
  "CYP27A1", "human", "NP_000775.1", "CYP27A1", "zebrafish", "XP_002663445.1",
  "CYP27B1", "human", "NP_000776.1", "CYP27B1", "zebrafish", "NP_001298720.1",
  "CYP24A1", "human", "NP_000773.2", "CYP24A1", "zebrafish", "NP_001082927.1",
  "VDR", "human", "NP_000367.1", "VDR", "zebrafish", "NP_570994.1", "VDR", "zebrafish", "XP_005166104.1")
query_fa <- file.path(ann_dir, "queries.faa")
if (!file.exists(query_fa)) {
  download.file(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=",
                       paste(queries$protein, collapse = ",")), query_fa, quiet = TRUE)
}
queries <- inner_join(queries, read_fasta(query_fa), by = "protein")
stopifnot(nrow(queries) == 15)
# Human genes that count as a reciprocal hit for each query gene
targets <- c(DHCR7 = "DHCR7", GC = "GC|ALB|AFP|AFM", CYP2R1 = "CYP2R1", CYP27A1 = "CYP27A1",
             CYP27B1 = "CYP27B1", CYP24A1 = "CYP24A1", VDR = "VDR")

fwd <- blastp(queries, "new") %>%
  inner_join(select(queries, qseqid = protein, query_gene), by = "qseqid") %>%
  inner_join(cds_of(ft$new), by = c(sseqid = "protein")) %>%
  group_by(query_gene, gene_id) %>%
  summarise(fwd_identity = max(pident), fwd_bits = max(bitscore), .groups = "drop")
candidates <- longest$new %>% filter(gene_id %in% fwd$gene_id)
rev <- best_hit(blastp(candidates, "human", max_target = 5), "human") %>%
  left_join(distinct(genes$human, gene_id, human_gene = symbol), by = "gene_id") %>%
  transmute(protein = qseqid, human_gene, rev_identity = pident)
paralog_search <- fwd %>%
  left_join(select(candidates, gene_id, protein), by = "gene_id") %>%
  left_join(rev, by = "protein") %>%
  filter(str_detect(human_gene, paste0("^(", targets[query_gene], ")$"))) %>%
  left_join(select(genes$new, gene_id, symbol, product, chr), by = "gene_id") %>%
  left_join(select(vd3, gene_id, name), by = "gene_id") %>%
  arrange(query_gene, desc(fwd_bits))
write_tsv(paralog_search, file.path(out_dir, "paralog_search.tsv"))
print(as.data.frame(select(paralog_search, query_gene, gene_id, name, symbol, chr, human_gene, fwd_identity, rev_identity)))
new_genes <- setdiff(paralog_search$gene_id, vd3$gene_id)
missed <- setdiff(vd3$gene_id, paralog_search$gene_id)
message("Candidates not in the gene list or LOC100136922: ", if (length(new_genes)) paste(new_genes, collapse = ", ") else "none")
message("Genes of the list without a reciprocal hit: ",
        if (length(missed)) paste(vd3$name[vd3$gene_id %in% missed], collapse = ", ") else "none")

# --- 3. Synteny of the paralog pairs in RS_2026_09
pairs <- tribble(~gene_a, ~gene_b,
                 "dhcr7.a", "dhcr7.b", "alb1", "LOC100136922", "alb1", "alb2", "cyp2r1.a", "cyp2r1.b",
                 "cyp27b1.a", "cyp27b1.b", "cyp24a1.a", "cyp24a1.b", "vdra.a", "vdra.b", "vdrb.a", "vdrb.b",
                 "cyp27a1.a", "cyp27a1.b", "cyp2r1.c", NA)
# A gene of the list by name, any other gene by GeneID
id_of <- function(x) if (x %in% vd3$name) vd3$gene_id[vd3$name == x] else x
coding_new <- genes$new %>% filter(gene_id %in% longest$new$gene_id) %>% arrange(seqid, start)
flanks <- function(name, n = 10) {
  g <- coding_new %>% filter(gene_id == id_of(name))
  coding_new %>% filter(seqid == g$seqid) %>% mutate(i = row_number()) %>%
    filter(abs(i - i[gene_id == g$gene_id]) %in% seq_len(n)) %>% pull(gene_id)
}
synteny_one <- function(a, b) {
  ga <- genes$new %>% filter(gene_id == id_of(a))
  hits <- blastp(filter(longest$new, gene_id %in% flanks(a)), "new", evalue = 1e-10, max_target = 50) %>%
    left_join(cds_of(ft$new), by = c(sseqid = "protein")) %>%
    left_join(select(genes$new, gene_id, hit_chr = chr, hit_start = start), by = "gene_id") %>%
    filter(hit_chr != ga$chr) %>%
    group_by(qseqid) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup()
  top <- names(sort(table(hits$hit_chr), decreasing = TRUE))[1]
  gb <- if (is.na(b)) NULL else genes$new %>% filter(gene_id == id_of(b))
  tibble(gene_a = a, chr_a = ga$chr, gene_b = b, chr_b = if (is.null(gb)) NA else gb$chr,
         flank_genes = 2 * 10, with_hit = nrow(hits), top_chr = top, on_top_chr = sum(hits$hit_chr == top),
         on_chr_b = if (is.null(gb)) NA else sum(hits$hit_chr == gb$chr),
         near_b_2mb = if (is.null(gb)) NA else sum(hits$hit_chr == gb$chr & abs(hits$hit_start - gb$start) < 2e6))
}
synteny <- bind_rows(pmap(pairs, function(gene_a, gene_b) {
  bind_rows(synteny_one(gene_a, gene_b), if (!is.na(gene_b)) synteny_one(gene_b, gene_a))
}))
write_tsv(synteny, file.path(out_dir, "synteny.tsv"))
print(as.data.frame(synteny))
