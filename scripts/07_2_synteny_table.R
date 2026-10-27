# Table S6: synteny of the paralog pairs in RS_2026_09 (fSalSal1), from the flanking-gene test of 07_1.
# For each copy, the proteins of the ten genes on each side were searched against the proteome, and the best match
# on another chromosome was counted. cyp2r1.c (excluded_genes) is left out. alb1 and alb2 are listed as a contrast:
# they are not Ss4R paralogs (Results), and the flanking genes of alb1 point to LOC100136922 instead.
# Input:  vd3_genes/annotation_check/synteny.tsv (07_1)
# Output: vd3_genes/annotation_check/table_s6.tsv, manuscript/v6/Table S6 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

check_dir <- file.path(vd3_dir, "annotation_check")
table_tsv <- file.path(check_dir, "table_s6.tsv")
table_file <- file.path(fig_dir, "Table S6 Synteny of VD3-related paralog pairs.xlsx")

# Pairs in the order of Table 1, the albumin pairs last
pair_order <- c("dhcr7", "cyp2r1", "cyp27a1", "cyp27b1", "cyp24a1", "vdra", "vdrb", "alb1/LOC100136922", "alb1/alb2")
note <- paste("Chromosomes (ssa) are those of the fSalSal1 assembly (RefSeq RS_2026_09). Flanking genes: the ten",
              "protein-coding genes on each side of the gene. For each, the best BLASTp match (E <= 1e-10) on another",
              "chromosome was counted. alb1 and alb2 are not Ss4R paralogs and are shown for comparison.")

synteny <- read_tsv(file.path(check_dir, "synteny.tsv"), col_types = cols(chr_a = "c", chr_b = "c", top_chr = "c")) %>%
  filter(!is.na(gene_b), !gene_a %in% excluded_genes, !gene_b %in% excluded_genes) %>%
  mutate(pair = case_when(gene_a == "alb1" ~ paste(gene_a, gene_b, sep = "/"),
                          gene_b == "alb1" ~ paste(gene_b, gene_a, sep = "/"),
                          .default = str_remove(gene_a, "\\.[ab]$")),
         pair = factor(pair, levels = pair_order)) %>%
  arrange(pair)
stopifnot(!anyNA(synteny$pair))

synteny %>%
  transmute(Gene = gene_a, Chromosome = paste0("ssa", chr_a), `Other copy` = gene_b,
            `Chromosome of other copy` = paste0("ssa", chr_b), `Flanking genes` = flank_genes,
            `With a match on another chromosome` = with_hit, `Most frequent chromosome` = paste0("ssa", top_chr),
            `Matches on the chromosome of the other copy` = on_chr_b,
            `Matches within 2 Mb of the other copy` = near_b_2mb) %>%
  mutate(across(everything(), as.character)) %>%
  add_row(Gene = c("", note)) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  write_tsv(table_tsv)
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(table_tsv),
                              shQuote(table_file))) == 0)
print(as.data.frame(select(synteny, gene_a, gene_b, on_chr_b, near_b_2mb)))
