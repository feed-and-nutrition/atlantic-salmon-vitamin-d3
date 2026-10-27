# Trees built from an accession list: gc (Fig 2b, the albumin family with the two Atlantic salmon albumins),
# vdr_extended (Fig 2d, one protein per gene, with chicken, spotted gar, medaka, fugu, stickleback and northern pike
# added to resolve the teleost vdra and vdrb lineages) and cyp27a1 (Fig S4a, the CYP27A1 orthologs with the salmon
# cyp27a2 and cyp27a3 genes, rooted on CYP27B1).
# Sequences are listed in msa/<gene>/<gene>_accessions.tsv and downloaded from NCBI (E-utilities) on the first run.
# Delete msa/<gene>/<gene>.fasta to download them again.
# Tree with one isoform per salmon gene: dhcr7 (Fig 2a), from the isoform set in msa/dhcr7/dhcr7.fasta (Fig S4d).
# Each salmon gene keeps its longest protein (the first listed when several are equally long), written to
# msa/dhcr7_longest/dhcr7_longest.fasta under the gene name (dhcr7.b.1 -> dhcr7.b).
# The alignment and tree are rebuilt only when the fasta is newer than the tree (delete <gene>.afa.treefile to
# force a rebuild).
# Needs muscle5 and iqtree2 on the PATH, and internet access on the first run.
# Input:  msa/<gene>/<gene>_accessions.tsv
# Output: msa/<gene>/<gene>.fasta, <gene>.afa (Muscle5), <gene>.afa.* (IQ-TREE 2, -m MFP -B 1000), the same in
#         msa/<gene>_longest/

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")
suppressPackageStartupMessages(library(Biostrings))

build_tree <- function(gene) {
  gene_dir <- file.path(msa_dir, gene)
  fasta <- file.path(gene_dir, paste0(gene, ".fasta"))
  acc <- read_tsv(file.path(gene_dir, paste0(gene, "_accessions.tsv")), show_col_types = FALSE)
  stopifnot(!anyDuplicated(acc$label), !any(str_detect(acc$label, "[ (),:;]")))

  if (!file.exists(fasta)) {
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=",
                  paste(acc$accession, collapse = ","))
    lines <- readLines(url, warn = FALSE)
    lines <- lines[lines != ""]
    # Replace each NCBI header with the label of its accession
    hdr <- which(startsWith(lines, ">"))
    got <- str_extract(lines[hdr], "^>(\\S+)", group = 1)
    stopifnot(setequal(got, acc$accession))
    lines[hdr] <- paste0(">", acc$label[match(got, acc$accession)])
    writeLines(lines, fasta)
    message("Downloaded ", length(hdr), " sequences to ", fasta)
  }

  align_and_build(gene_dir, gene)
}

align_and_build <- function(gene_dir, name) {
  fasta <- file.path(gene_dir, paste0(name, ".fasta"))
  afa <- file.path(gene_dir, paste0(name, ".afa"))
  treefile <- paste0(afa, ".treefile")
  if (!file.exists(treefile) || file.mtime(fasta) > file.mtime(treefile)) {
    old <- list.files(gene_dir, pattern = paste0("^", name, "\\.afa"), full.names = TRUE)
    file.remove(old)
    run <- function(cmd, args) if (system2(cmd, args) != 0) stop(cmd, " failed")
    run("muscle5", c("-align", shQuote(fasta), "-output", shQuote(afa)))
    old_wd <- setwd(gene_dir)
    on.exit(setwd(old_wd))
    run("iqtree2", c("-s", basename(afa), "-m", "MFP", "-B", "1000", "--seed", "1", "-T", "AUTO", "--quiet"))
    message("Tree written to ", treefile)
  }
}

build_longest <- function(gene) {
  out_dir <- file.path(msa_dir, paste0(gene, "_longest"))
  fasta <- file.path(out_dir, paste0(gene, "_longest.fasta"))
  if (!file.exists(fasta)) {
    seqs <- readAAStringSet(file.path(msa_dir, gene, paste0(gene, ".fasta")))
    keep <- tibble(label = names(seqs), length = width(seqs), order = seq_along(seqs)) %>%
      mutate(gene = ifelse(startsWith(label, "Salmon_"), str_remove(label, "(?<=\\.[ab])\\.\\d+$"), label)) %>%
      arrange(desc(length), order) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(order)
    out <- seqs[keep$order]
    names(out) <- keep$gene
    dir.create(out_dir, showWarnings = FALSE)
    writeXStringSet(out, fasta)
    message("Kept ", paste(keep$label[keep$label != keep$gene], collapse = ", "), " in ", fasta)
  }
  align_and_build(out_dir, paste0(gene, "_longest"))
}

walk(c("gc", "vdr_extended", "cyp27a1"), build_tree)
build_longest("dhcr7")
