# Table S11: sequence similarity (06_1), sequence divergence and expression profiles of the Ss4R paralog pairs.
# Divergence: the proteins of each pair (the transcripts with the longest CDS, cached by 06_1) are aligned with
# MUSCLE 5 and the CDS are aligned codon by codon on the protein alignment. dN, dS and dN/dS are estimated with yn00
# (PAML, Yang and Nielsen 2000). ORF: start and stop codon, no internal stop. Motifs: the heme-binding motif of
# cytochrome P450s (FxxGx(R/H/K)xCxG) and the two zinc fingers of the VDR DNA-binding domain (Cx2Cx13Cx2C,
# Cx5Cx9Cx2C), and the net charge of the first 30 residues (mitochondrial presequences are positive), each also for
# the pike ortholog (06_3). The yn00 results are cached in divergence.tsv (delete it to recompute).
# Rediploidization: pairs on one of the seven pairs of homeologous arms with delayed rediploidization (LORe,
# Robertson et al. 2017: 2p-5q, 2q-12qa, 3q-6p, 4p-8q, 7q-17qb, 11qa-26, 16qb-17qa) by their fSalSal1 chromosomes
# (07_1), the others AORe. Matched by chromosome only, which is enough for these pairs.
# Profiles: Pearson correlation of log2(TPM + 1) across tissues, between the two copies (atlas, 14 tissues, and the
# mean of the four PRJNA359635 pools per tissue, 10 tissues), and between each copy, or the sum of both, and the pike
# ortholog in the tissues shared with the pike atlas (06_3). Genes below 1 TPM in every tissue have no correlation.
# Needs muscle5 and yn00 on the PATH (or PAML_YN00 set to yn00) only when divergence.tsv is missing.
# Input:  vd3_genes/paralogs/ (06_1, 06_3), vd3_genes/annotation_check/synteny.tsv (07_1),
#         vd3_genes/tpm_vd3_genes.tsv, vd3_genes/tpm_vd3_genes_studies.tsv
# Output: vd3_genes/paralogs/divergence.tsv, profiles.tsv, table_s9_*.tsv (sheets), manuscript/v7/Table S11 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_dir <- file.path(vd3_dir, "paralogs")
seq_dir <- file.path(out_dir, "seq")
divergence_file <- file.path(out_dir, "divergence.tsv")
profiles_file <- file.path(out_dir, "profiles.tsv")
table_file <- file.path(fig_dir, paste("Table S11 Sequence similarity, divergence, and expression profiles of",
                                      "VD3-related paralog pairs.xlsx"))
pairs <- c("dhcr7", "cyp2r1", "cyp27a1", "cyp27b1", "cyp24a1", "vdra", "vdrb")

read_fasta <- function(file) {
  lines <- readLines(file)
  head <- grepl("^>", lines)
  tibble(header = sub("^>", "", lines[head]),
         seq = map_chr(split(lines[!head], cumsum(head)[!head]), paste, collapse = ""))
}
# One transcript per gene: the one with the longest CDS, as in 06_1
genes <- read_tsv(file.path(vd3_dir, "salmon_vd3_orthologs_v1.tsv"), col_types = cols(.default = "c")) %>%
  transmute(name, transcript = str_split(str_trim(`mRNA refseq`), ",")) %>%
  unnest(transcript)
coding <- function(file, col) {
  read_fasta(file) %>% transmute(transcript = str_match(header, "^lcl\\|(.+?)_(?:cds|prot)_")[, 2], !!col := seq)
}
rep_seq <- coding(file.path(seq_dir, "cds.fna"), "cds") %>%
  inner_join(coding(file.path(seq_dir, "protein.faa"), "protein"), by = "transcript") %>%
  inner_join(genes, by = "transcript") %>%
  group_by(name) %>%
  slice_max(nchar(cds), n = 1, with_ties = FALSE) %>%
  ungroup()
copies <- paste0(rep(pairs, each = 2), c(".a", ".b"))
stopifnot(all(copies %in% rep_seq$name))
seq_of <- function(gene, type) rep_seq[[type]][rep_seq$name == gene]

# --- Divergence (yn00)
codon_alignment <- function(a, b) {
  faa <- tempfile(fileext = ".faa")
  afa <- tempfile(fileext = ".afa")
  writeLines(paste0(">", c(a, b), "\n", c(seq_of(a, "protein"), seq_of(b, "protein"))), faa)
  stopifnot(system2("muscle5", c("-align", faa, "-output", afa), stdout = FALSE, stderr = FALSE) == 0)
  aln <- read_fasta(afa)
  map_chr(c(a, b), function(g) {
    aa <- strsplit(aln$seq[aln$header == g], "")[[1]]
    cds <- sub("(TAA|TAG|TGA)$", "", seq_of(g, "cds"))
    codons <- substring(cds, seq(1, nchar(cds), 3), seq(3, nchar(cds), 3))
    stopifnot(length(codons) == sum(aa != "-"))
    out <- rep("---", length(aa))
    out[aa != "-"] <- codons
    paste(out, collapse = "")
  })
}
yn00 <- function(a, b) {
  aln <- codon_alignment(a, b)
  wd <- tempfile()
  dir.create(wd)
  writeLines(c(paste(2, nchar(aln[1])), paste0("seq_a  ", aln[1]), paste0("seq_b  ", aln[2])), file.path(wd, "aln.phy"))
  writeLines(c("seqfile = aln.phy", "outfile = out.txt", "verbose = 0", "icode = 0", "weighting = 0",
               "commonf3x4 = 0"), file.path(wd, "yn00.ctl"))
  owd <- setwd(wd)
  status <- system2(Sys.getenv("PAML_YN00", "yn00"), "yn00.ctl", stdout = FALSE)
  setwd(owd)
  stopifnot(status == 0)
  # Yang and Nielsen (2000) block: seq. seq. S N t kappa omega dN +- SE dS +- SE
  out <- readLines(file.path(wd, "out.txt"))
  values <- scan(text = out[grep("^seq\\. seq\\.", out)[1] + 2], what = "", quiet = TRUE)
  values <- as.numeric(values[values != "+-"])
  tibble(dN = values[8], dN_se = values[9], dS = values[10], dS_se = values[11], dN_dS = values[7])
}
if (!file.exists(divergence_file)) {
  map_dfr(set_names(pairs), ~ yn00(paste0(.x, ".a"), paste0(.x, ".b")), .id = "pair") %>%
    write_tsv(divergence_file)
}
divergence <- read_tsv(divergence_file, show_col_types = FALSE)
lore_pairs <- c("2-5", "2-12", "3-6", "4-8", "7-17", "11-26", "16-17")
chromosomes <- read_tsv(file.path(vd3_dir, "annotation_check", "synteny.tsv"), col_types = cols(.default = "c")) %>%
  filter(gene_a %in% paste0(pairs, ".a"), gene_b == str_replace(gene_a, "\\.a$", ".b")) %>%
  transmute(pair = str_remove(gene_a, "\\.a$"), chromosomes = paste0("ssa", chr_a, " / ssa", chr_b),
            rediploidization = ifelse(map2_chr(as.integer(chr_a), as.integer(chr_b), ~ paste(sort(c(.x, .y)), collapse = "-"))
                                      %in% lore_pairs, "LORe", "AORe"))
stopifnot(setequal(chromosomes$pair, pairs))

# --- ORF and motifs, for the salmon copies and the pike ortholog
pike <- read_tsv(file.path(out_dir, "pike_orthologs.tsv"), col_types = cols(.default = "c"))
pike_faa <- file.path(proj_dir, "annotation", "gc_search", "GCF_011004845.1_fEsoLuc1.pri_protein.faa.gz")
pike_seq <- if (file.exists(pike_faa)) {
  read_fasta(gzfile(pike_faa)) %>%
    transmute(pike_protein = str_extract(header, "^\\S+"), protein = seq) %>%
    filter(pike_protein %in% pike$pike_protein)
}
motif_of <- function(pair) {
  if (str_detect(pair, "^cyp")) return(c(`heme-binding motif` = "F..G.[RHK].C.G"))
  if (str_detect(pair, "^vdr")) return(c(`zinc finger 1` = "C..C.{13}C..C", `zinc finger 2` = "C.{5}C.{9}C..C"))
  character()
}
net_charge <- function(p) sum(str_count(substr(p, 1, 30), c("[RK]", "[DE]")) * c(1, -1))
motifs <- function(protein, pair) {
  m <- motif_of(pair)
  if (!length(m)) return("not tested")
  paste(paste(names(m), ifelse(str_detect(protein, m), "present", "absent")), collapse = ", ")
}
orf <- tibble(gene = copies, pair = str_remove(gene, "\\.[ab]$")) %>%
  mutate(cds = map_chr(gene, seq_of, "cds"), protein = map_chr(gene, seq_of, "protein"),
         codons = map(cds, ~ substring(.x, seq(1, nchar(.x) - 2, 3), seq(3, nchar(.x), 3))),
         complete = map_lgl(codons, ~ .x[1] == "ATG" && tail(.x, 1) %in% c("TAA", "TAG", "TGA") &&
                              !any(head(.x, -1) %in% c("TAA", "TAG", "TGA"))) & nchar(cds) %% 3 == 0,
         length = nchar(protein), motifs = map2_chr(protein, pair, motifs),
         n_charge = map_int(protein, net_charge))
stopifnot(all(orf$complete))
pike_motifs <- pike %>%
  filter(str_remove(name, "\\.[ab]$") %in% pairs) %>%
  distinct(pair = str_remove(name, "\\.[ab]$"), pike_symbol, pike_protein) %>%
  left_join(pike_seq, by = "pike_protein") %>%
  mutate(pike_length = nchar(protein), pike_motifs = map2_chr(protein, pair, motifs),
         pike_n_charge = map_int(protein, net_charge)) %>%
  select(-protein)

# --- Expression profiles
atlas <- read_gene_tpm(file.path(vd3_dir, "tpm_vd3_genes.tsv"))$log_tpm
pool <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes_studies.tsv"), show_col_types = FALSE) %>%
  filter(study == "PRJNA359635", !excluded) %>%
  group_by(name, tissue = tissue_label(tissue)) %>%
  summarise(log_tpm = mean(log2(tpm + 1)), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = log_tpm) %>%
  column_to_rownames("name") %>%
  as.matrix()
pike_tpm <- read_tsv(file.path(out_dir, "pike_tpm.tsv"), col_types = cols(pike_gene_id = "c"))
pike_log <- log2(as.matrix(select(pike_tpm, -pike_gene_id, -pike_symbol)) + 1)
rownames(pike_log) <- pike_tpm$pike_gene_id
colnames(pike_log) <- tissue_label(colnames(pike_log))
# Correlation of two profiles over the given tissues, NA when either stays below 1 TPM
profile_r <- function(x, y) if (max(x) < 1 || max(y) < 1) NA_real_ else round(cor(x, y), 2)
# Table 2 pattern of a log2(TPM + 1) profile (expr_pattern, 00_config)
pattern_of <- function(x) expr_pattern(sd(x), tau_index(x), 2^max(x) - 1)
pike_id <- function(pair) unique(pike$pike_gene_id[pike$name %in% paste0(pair, c(".a", ".b"))])
profiles <- map_dfr(pairs, function(pair) {
  a <- paste0(pair, ".a")
  b <- paste0(pair, ".b")
  p <- pike_log[pike_id(pair), ]
  shared <- intersect(colnames(atlas), names(p))
  shared_pool <- intersect(colnames(pool), names(p))
  sum_log <- function(m, t) log2(2^m[a, t] + 2^m[b, t] - 1)
  tibble(pair, gene_a = a, gene_b = b,
         peak_a = round(2^max(atlas[a, ]) - 1, 1), peak_b = round(2^max(atlas[b, ]) - 1, 1),
         pike_peak = round(2^max(p) - 1, 1), pike_peak_tissue = names(p)[which.max(p)],
         pattern_a = pattern_of(atlas[a, ]), pattern_b = pattern_of(atlas[b, ]), pike_pattern = pattern_of(p),
         a_higher = sum(atlas[a, ] - atlas[b, ] >= 1 & atlas[a, ] > 1),
         b_higher = sum(atlas[b, ] - atlas[a, ] >= 1 & atlas[b, ] > 1),
         r_ab_atlas = profile_r(atlas[a, ], atlas[b, ]), r_ab_pool = profile_r(pool[a, ], pool[b, ]),
         n_pike = length(shared), r_a_pike = profile_r(atlas[a, shared], p[shared]),
         r_b_pike = profile_r(atlas[b, shared], p[shared]), r_sum_pike = profile_r(sum_log(atlas, shared), p[shared]),
         n_pike_pool = length(shared_pool), r_a_pike_pool = profile_r(pool[a, shared_pool], p[shared_pool]),
         r_b_pike_pool = profile_r(pool[b, shared_pool], p[shared_pool]),
         r_sum_pike_pool = profile_r(sum_log(pool, shared_pool), p[shared_pool]))
})
write_tsv(profiles, profiles_file)
print(as.data.frame(left_join(divergence, select(profiles, pair, starts_with("r_")), by = "pair")))

# --- Table S11: Similarity (06_1), Divergence and Profiles
divergence_tsv <- file.path(out_dir, "table_s9_divergence.tsv")
profiles_tsv <- file.path(out_dir, "table_s9_profiles.tsv")
by_copy <- function(col) set_names(orf[[col]], orf$gene)
divergence %>%
  left_join(pike_motifs, by = "pair") %>%
  left_join(chromosomes, by = "pair") %>%
  transmute(`Gene a` = paste0(pair, ".a"), `Gene b` = paste0(pair, ".b"), `Chromosomes, a / b` = chromosomes,
            Rediploidization = rediploidization,
            dN = sprintf("%.3f (%.3f)", dN, dN_se), dS = sprintf("%.3f (%.3f)", dS, dS_se),
            `dN/dS` = sprintf("%.2f", dN_dS),
            `Protein length, a / b / pike (residues)` = paste(by_copy("length")[`Gene a`], by_copy("length")[`Gene b`],
                                                               pike_length, sep = " / "),
            `Complete ORF, a and b` = ifelse(by_copy("complete")[`Gene a`] & by_copy("complete")[`Gene b`], "yes", "no"),
            `Motifs, gene a` = by_copy("motifs")[`Gene a`], `Motifs, gene b` = by_copy("motifs")[`Gene b`],
            `Motifs, pike` = pike_motifs,
            `Net charge of residues 1 to 30, a / b / pike` = paste(by_copy("n_charge")[`Gene a`],
                                                                    by_copy("n_charge")[`Gene b`], pike_n_charge,
                                                                    sep = " / "),
            `Pike ortholog` = paste0(pike_symbol, " (", pike_protein, ")")) %>%
  add_row(`Gene a` = c("", paste("dN and dS (standard error) by the method of Yang and Nielsen (2000), yn00 in PAML,",
                                 "on the codon alignment of the transcripts with the longest CDS (Similarity sheet).",
                                 "Heme-binding motif FxxGx(R/H/K)xCxG of cytochrome P450s, zinc fingers Cx2Cx13Cx2C",
                                 "and Cx5Cx9Cx2C of the VDR DNA-binding domain. Net charge: R and K minus D and E;",
                                 "mitochondrial presequences are positively charged. Pike ortholog: best BLASTp hit of",
                                 "both copies in the northern pike proteome (fEsoLuc1). Chromosomes of fSalSal1.",
                                 "Rediploidization: LORe if the chromosomes match one of the seven pairs of homeologous",
                                 "arms with delayed rediploidization in Robertson et al. (2017), AORe otherwise."))) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  write_tsv(divergence_tsv)
fmt <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))
profiles %>%
  transmute(`Gene a` = gene_a, `Gene b` = gene_b, `Atlas peak TPM, a` = peak_a, `Atlas peak TPM, b` = peak_b,
            `Atlas tissues with a at least 2-fold higher` = a_higher,
            `Atlas tissues with b at least 2-fold higher` = b_higher,
            `r, a vs b, atlas` = fmt(r_ab_atlas), `r, a vs b, PRJNA359635` = fmt(r_ab_pool),
            `Pattern, a / b` = paste(pattern_a, pattern_b, sep = " / "),
            `Pike peak TPM (tissue)` = paste0(pike_peak, " (", tissue_label(pike_peak_tissue), ")"),
            `Pike pattern` = pike_pattern,
            `r, a vs pike, atlas` = fmt(r_a_pike), `r, b vs pike, atlas` = fmt(r_b_pike),
            `r, a + b vs pike, atlas` = fmt(r_sum_pike),
            `r, a vs pike, PRJNA359635` = fmt(r_a_pike_pool), `r, b vs pike, PRJNA359635` = fmt(r_b_pike_pool),
            `r, a + b vs pike, PRJNA359635` = fmt(r_sum_pike_pool)) %>%
  mutate(across(everything(), as.character)) %>%
  add_row(`Gene a` = c("", paste0("r: Pearson correlation of log2(TPM + 1) across tissues, not given when a gene stays ",
                                  "below 1 TPM in all of them. Atlas: 14 tissues (a vs b) and the ", profiles$n_pike[1],
                                  " tissues shared with the pike samples of PRJNA221548. PRJNA359635: mean of the four ",
                                  "pools per tissue, 10 tissues (a vs b) and ", profiles$n_pike_pool[1], " shared with ",
                                  "pike. a + b: the summed TPM of both copies. 2-fold higher: a difference of at least 1 ",
                                  "in log2(TPM + 1), in tissues above 1 TPM. Pattern: the rule of Table 2 (restricted if ",
                                  "SD >= 1 or tau >= 0.8, not classified below 1 TPM), for pike over its 13 tissues."))) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  write_tsv(profiles_tsv)
pike_tsv <- file.path(out_dir, "table_s9_pike.tsv")
pike %>%
  mutate(pair = str_remove(name, "\\.[ab]$")) %>%
  filter(pair %in% pairs) %>%
  group_by(pair, pike_gene_id, pike_symbol) %>%
  summarise(`Salmon genes` = paste(name, collapse = ", "), .groups = "drop") %>%
  inner_join(pike_tpm, by = c("pike_gene_id", "pike_symbol")) %>%
  arrange(match(pair, pairs)) %>%
  select(-pair) %>%
  rename(`Pike gene` = pike_symbol, `Gene ID` = pike_gene_id) %>%
  relocate(`Salmon genes`) %>%
  rename_with(tissue_label, -c(`Salmon genes`, `Pike gene`, `Gene ID`)) %>%
  mutate(across(everything(), as.character)) %>%
  add_row(`Salmon genes` = c("", paste("TPM of the pike orthologs in the tissue RNA-seq of PRJNA221548 (one sample per",
                                       "tissue), from the first 4 million read pairs of each run, quantified with Salmon",
                                       "against the fEsoLuc1 RefSeq transcripts."))) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  write_tsv(pike_tsv)
sheet_dir <- tempfile()
dir.create(sheet_dir)
sheets <- file.path(sheet_dir, c("Similarity.tsv", "Divergence.tsv", "Profiles.tsv", "Pike TPM.tsv"))
invisible(file.copy(c(file.path(out_dir, "table_s9.tsv"), divergence_tsv, profiles_tsv, pike_tsv), sheets))
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(sheets),
                              shQuote(table_file))) == 0)
