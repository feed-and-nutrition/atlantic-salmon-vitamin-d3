# Zebrafish as the teleost reference for the gene set (v11, reviewer 1 major comment 3): were any teleost genes
# missed by starting from the human genes? Not part of run_all.R. Uses the downloads of 07_1 (RS_2026_09, human
# GRCh38.p14) and 07_3 (zebrafish GRCz13ab) in annotation/ and needs blastp and makeblastdb.
# 1. Zebrafish genes: the human protein (longest isoform) of each canonical gene is searched with BLASTp
#    (E <= 1e-10) against the GRCz13ab proteome. A zebrafish gene is kept when its longest protein has its best human
#    hit in the query gene. This keeps teleost-specific copies (vdra/vdrb, the cyp27a group) and drops albumins.
# 2. Salmon genes: each zebrafish protein is searched against the RS_2026_09 proteome. A salmon gene is assigned to
#    the zebrafish gene that is the best zebrafish hit of its longest protein (so two Ss4R copies can share one
#    zebrafish gene), when its best human hit is also the query gene. Without this, the salmon albumins would be
#    assigned to zebrafish gc, their best zebrafish hit. For a zebrafish gene without a salmon ortholog, the best
#    salmon hit is given for reference.
# 3. Relationship per zebrafish gene: 1:0, 1:1 or 1:n salmon genes, with the manuscript names of the salmon genes.
# Input:  annotation/{RS_2026_09,human_GRCh38.p14,gc_search}/*_{protein.faa,feature_table.txt}.gz,
#         vd3_genes/salmon_vd3_orthologs_v1.tsv, msa/cyp27a1/cyp27a1_accessions.tsv
# Output: vd3_genes/zebrafish_orthologs/{zebrafish_genes,orthologs}.tsv

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

ann_dir <- file.path(proj_dir, "annotation")
out_dir <- file.path(vd3_dir, "zebrafish_orthologs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

query_genes <- c("DHCR7", "CYP2R1", "CYP27A1", "CYP27B1", "CYP24A1", "VDR", "GC")
prefix <- c(human = file.path(ann_dir, "human_GRCh38.p14", "GCF_000001405.40_GRCh38.p14"),
            zebrafish = file.path(ann_dir, "gc_search", "GCF_052040795.2_GRCz13ab"),
            salmon = file.path(ann_dir, "RS_2026_09", "GCF_965601325.2_fSalSal1"))
missing <- !file.exists(paste0(prefix, "_protein.faa.gz")) | !file.exists(paste0(prefix, "_feature_table.txt.gz"))
if (any(missing)) stop("missing downloads for ", paste(names(prefix)[missing], collapse = ", "), ": run 07_1 and 07_3 first")

read_features <- function(p) {
  file <- paste0(p, "_feature_table.txt.gz")
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
ft <- map(prefix, read_features)
genes <- map(ft, ~ .x %>% filter(feature == "gene") %>% distinct(gene_id = GeneID, .keep_all = TRUE) %>%
               select(gene_id, symbol, chr = chromosome))
cds <- map(ft, ~ .x %>% filter(feature == "CDS", !is.na(product_accession)) %>%
             distinct(protein = product_accession, gene_id = GeneID, product = name))
prot <- map(names(prefix), ~ read_fasta(paste0(prefix[[.x]], "_protein.faa.gz")) %>%
              inner_join(cds[[.x]], by = "protein")) %>% set_names(names(prefix))
longest <- map(prot, ~ .x %>% group_by(gene_id) %>% slice_max(nchar(seq), n = 1, with_ties = FALSE) %>% ungroup())

blast_dir <- file.path(ann_dir, "blastdb")
dir.create(blast_dir, showWarnings = FALSE)
make_db <- function(sp) {
  db <- file.path(blast_dir, paste0("orthologs_", sp))
  if (!file.exists(paste0(db, ".pin")) && !file.exists(paste0(db, ".00.pin"))) {
    fa <- tempfile(fileext = ".faa")
    writeLines(paste0(">", prot[[sp]]$protein, "\n", prot[[sp]]$seq), fa)
    stopifnot(system2("makeblastdb", c("-in", fa, "-dbtype", "prot", "-out", shQuote(db)), stdout = FALSE) == 0)
  }
  db
}
blastp <- function(query, sp, max_target = 500) {
  fa <- tempfile(fileext = ".faa")
  writeLines(paste0(">", query$protein, "\n", query$seq), fa)
  out <- system2("blastp", c("-query", fa, "-db", shQuote(make_db(sp)), "-evalue", "1e-10",
                             "-max_target_seqs", max_target, "-num_threads", 4,
                             "-outfmt", shQuote("6 qseqid sseqid pident length qlen slen evalue bitscore")),
                 stdout = TRUE)
  read_tsv(I(paste(c("qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore", out), collapse = "\n")),
           col_types = "ccdiiidd") %>%
    inner_join(select(cds[[sp]], sseqid = protein, s_gene = gene_id), by = "sseqid")
}
# Best subject gene per query protein
best <- function(hits) hits %>% group_by(qseqid) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup()
human_symbol <- function(ids) genes$human$symbol[match(ids, genes$human$gene_id)]

# --- 1. Zebrafish genes of the canonical pathway
human_query <- longest$human %>%
  inner_join(filter(genes$human, symbol %in% query_genes) %>% select(gene_id, query = symbol), by = "gene_id")
stopifnot(setequal(human_query$query, query_genes))
zf_candidates <- blastp(human_query, "zebrafish") %>% distinct(s_gene) %>% pull()
zf_back <- blastp(filter(longest$zebrafish, gene_id %in% zf_candidates), "human", 5) %>% best()
zf_genes <- longest$zebrafish %>%
  inner_join(transmute(zf_back, protein = qseqid, best_human = human_symbol(s_gene), pident_human = pident),
             by = "protein") %>%
  filter(best_human %in% query_genes) %>%
  left_join(select(genes$zebrafish, gene_id, symbol, chr), by = "gene_id") %>%
  transmute(query = best_human, zf_gene_id = gene_id, zf_symbol = symbol, zf_chr = chr, zf_protein = protein,
            zf_product = product, pident_human) %>%
  arrange(match(query, query_genes), zf_symbol)
write_tsv(zf_genes, file.path(out_dir, "zebrafish_genes.tsv"))

# --- 2. Salmon genes assigned to each zebrafish gene
zf_query <- filter(longest$zebrafish, gene_id %in% zf_genes$zf_gene_id)
sal_hits <- blastp(zf_query, "salmon")
sal_candidates <- unique(sal_hits$s_gene)
sal_query <- filter(longest$salmon, gene_id %in% sal_candidates)
sal_to_zf <- blastp(sal_query, "zebrafish", 5) %>% best() %>% select(protein = qseqid, zf_gene_id = s_gene, pident_zf = pident)
sal_to_human <- blastp(sal_query, "human", 5) %>% best() %>% transmute(protein = qseqid, sal_best_human = human_symbol(s_gene))
names_ms <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
  select(gene_id, name_ms = name)
names_cyp27a <- read_tsv(file.path(proj_dir, "msa", "cyp27a1", "cyp27a1_accessions.tsv"), col_names = FALSE,
                         col_types = cols(.default = "c")) %>%
  filter(str_starts(X1, "Salmon_LOC")) %>% transmute(gene_id = X3, name_ms = paste0(X1, " (Table S11)"))
salmon <- sal_query %>%
  inner_join(sal_to_zf, by = "protein") %>%
  inner_join(sal_to_human, by = "protein") %>%
  inner_join(select(zf_genes, zf_gene_id, query), by = "zf_gene_id") %>%
  filter(sal_best_human == query) %>%  # as in the main search, the best human hit must be the query gene
  left_join(select(genes$salmon, gene_id, sal_symbol = symbol, sal_chr = chr), by = "gene_id") %>%
  left_join(bind_rows(names_ms, names_cyp27a) %>% distinct(gene_id, .keep_all = TRUE), by = "gene_id") %>%
  transmute(zf_gene_id, sal_gene_id = gene_id, sal_symbol, sal_chr, name_ms, sal_best_human, pident_zf)

# Best salmon hit of zebrafish genes without a reciprocal salmon gene
nonrecip <- sal_hits %>% best() %>%
  inner_join(select(zf_query, qseqid = protein, zf_gene_id = gene_id), by = "qseqid") %>%
  left_join(select(genes$salmon, s_gene = gene_id, best_salmon_hit = symbol), by = "s_gene") %>%
  select(zf_gene_id, best_salmon_hit, best_salmon_pident = pident)

orthologs <- zf_genes %>%
  select(query, zf_gene_id, zf_symbol, zf_chr, pident_human) %>%
  left_join(salmon, by = "zf_gene_id", relationship = "one-to-many") %>%
  group_by(zf_gene_id) %>%
  mutate(n_salmon = sum(!is.na(sal_gene_id)),
         relationship = case_when(n_salmon == 0 ~ "1:0", n_salmon == 1 ~ "1:1", TRUE ~ paste0("1:", n_salmon))) %>%
  ungroup() %>%
  left_join(nonrecip, by = "zf_gene_id") %>%
  mutate(across(c(best_salmon_hit, best_salmon_pident), ~ if_else(n_salmon == 0, .x, NA))) %>%
  arrange(match(query, query_genes), zf_symbol, sal_chr)
write_tsv(orthologs, file.path(out_dir, "orthologs.tsv"))

print(as.data.frame(select(orthologs, query, zf_symbol, relationship, sal_symbol, sal_chr, name_ms, sal_best_human,
                           best_salmon_hit)), right = FALSE)
# Salmon genes of the manuscript set that no zebrafish gene claimed
unclaimed <- anti_join(names_ms, salmon, by = c(gene_id = "sal_gene_id"))
if (nrow(unclaimed)) {
  message("Manuscript genes without a zebrafish ortholog in this search:")
  print(as.data.frame(unclaimed))
}
