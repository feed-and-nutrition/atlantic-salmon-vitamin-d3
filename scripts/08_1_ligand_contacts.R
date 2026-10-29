# Are the ligand-binding residues of the human proteins kept in the salmon genes?
# Contact residues: protein residues with a heavy atom within 4 A of the ligand in a human crystal structure
# (VDR with calcitriol 1DB1, CYP2R1 with vitamin D3 3C6G, DBP/GC with 25(OH)D3 1J78). A contact counts as a hydrogen
# bond when a protein N or O atom lies within 3.5 A of a ligand O atom. The PDB residues are mapped onto the human
# sequence of each phylogeny alignment in msa/ by the numbering offset that matches the modelled chain best (1J78 is
# numbered from the mature protein, 16 residues on), and the salmon residues are read from the
# alignment column (all isoforms of a gene, which agree at every contact). PDB files are cached in
# msa/ligand_contacts/pdb/ (internet access on the first run).
# Input:  msa/vdr/vdr.afa, msa/cyp2r1/cyp2r1.afa, msa/gc/gc.afa
# Output: msa/ligand_contacts/contacts.tsv (one row per contact residue), summary.tsv (per sequence),
#         <protein>.tsv (one sheet each) and manuscript/v7/Table S7 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")
suppressPackageStartupMessages(library(Biostrings))

out_dir <- file.path(msa_dir, "ligand_contacts")
pdb_dir <- file.path(out_dir, "pdb")
dir.create(pdb_dir, recursive = TRUE, showWarnings = FALSE)

structures <- tribble(~protein, ~pdb, ~chain, ~ligand, ~msa, ~human,
                      "VDR", "1DB1", "A", "VDX", "vdr", "Human_VDR",
                      "CYP2R1", "3C6G", "A", "VD3", "cyp2r1", "Human_CYP2R1",
                      "GC", "1J78", "B", "VDY", "gc", "Human_GC")
table_file <- file.path(fig_dir, "Table S7 Ligand-contact residues of VDR, CYP2R1, and GC.xlsx")
# Sequences shown after the human one, salmon last: the other orthologs, and for GC the true fish gc genes, the
# salmonid proteins annotated as GC (albumins, Figure 2b) and human ALB. LOC106586122 is cyp2r1.c (excluded_genes)
shown <- list(VDR = c("Mouse_Vdr", "Zebrafish_vdra", "Zebrafish_vdrb",
                      "Salmon_vdra.a", "Salmon_vdra.b", "Salmon_vdrb.a", "Salmon_vdrb.b"),
              CYP2R1 = c("Mouse_Cyp2r1", "Zebrafish_cyp2r1", "Salmon_cyp2r1.a", "Salmon_cyp2r1.b", "Salmon_cyp2r1.c"),
              GC = c("Mouse_Gc", "Chicken_GC", "Gar_gc", "Zebrafish_gc", "Trout_gc", "Chum_gc", "Pink_gc", "Human_ALB",
                     "Salmon_alb1", "Salmon_alb2", "Salmon_LOC100136922"))
species <- c(Human = "Human", Mouse = "Mouse", Chicken = "Chicken", Gar = "Spotted gar", Zebrafish = "Zebrafish",
             Trout = "Rainbow trout", Chum = "Chum salmon", Pink = "Pink salmon")
display_name <- function(x) {
  sp <- str_extract(x, "^[^_]+")
  gene <- str_remove(x, "^[^_]+_")
  case_when(x == "Salmon_cyp2r1.c" ~ "LOC106586122", sp == "Salmon" ~ gene, .default = paste(species[sp], gene))
}

read_pdb <- function(id) {
  file <- file.path(pdb_dir, paste0(id, ".pdb"))
  if (!file.exists(file)) download.file(paste0("https://files.rcsb.org/download/", id, ".pdb"), file, quiet = TRUE)
  lines <- readLines(file)
  lines <- lines[startsWith(lines, "ATOM  ") | startsWith(lines, "HETATM")]
  tibble(record = str_trim(substr(lines, 1, 6)), atom = str_trim(substr(lines, 13, 16)),
         alt = substr(lines, 17, 17), resn = str_trim(substr(lines, 18, 20)), chain = substr(lines, 22, 22),
         resi = as.integer(substr(lines, 23, 26)), x = as.numeric(substr(lines, 31, 38)),
         y = as.numeric(substr(lines, 39, 46)), z = as.numeric(substr(lines, 47, 54)),
         element = str_trim(substr(lines, 77, 78))) %>%
    filter(alt %in% c(" ", "A"), element != "H")
}

read_afa <- function(gene) {
  aln <- readAAStringSet(file.path(msa_dir, gene, paste0(gene, ".afa")))
  set_names(as.character(aln), names(aln))
}

contacts_of <- function(s) {
  atoms <- read_pdb(s$pdb) %>% filter(chain == s$chain)
  lig <- filter(atoms, record == "HETATM", resn == s$ligand)
  prot <- filter(atoms, record == "ATOM")
  d <- sqrt(outer(prot$x, lig$x, "-")^2 + outer(prot$y, lig$y, "-")^2 + outer(prot$z, lig$z, "-")^2)
  polar <- outer(prot$element %in% c("N", "O"), lig$element == "O", "&")
  prot %>%
    mutate(min_dist = apply(d, 1, min), hbond = rowSums(polar & d <= 3.5) > 0) %>%
    group_by(resi, resn) %>%
    summarise(min_dist = min(min_dist), hbond = any(hbond), .groups = "drop") %>%
    filter(min_dist <= 4) %>%
    mutate(aa = names(AMINO_ACID_CODE)[match(str_to_title(resn), AMINO_ACID_CODE)])
}

# One-letter sequence of the chain as modelled, with the residue number of each letter
chain_seq <- function(s) {
  read_pdb(s$pdb) %>%
    filter(record == "ATOM", chain == s$chain) %>%
    distinct(resi, resn) %>%
    mutate(aa = names(AMINO_ACID_CODE)[match(str_to_title(resn), AMINO_ACID_CODE)])
}

map_contacts <- function(s) {
  aln <- read_afa(s$msa)
  human_gapped <- aln[[s$human]]
  human <- gsub("-", "", human_gapped)
  modelled <- chain_seq(s)
  # Offset from PDB to human numbering: the one under which most modelled residues match the human sequence
  human_aa <- strsplit(human, "")[[1]]
  match_rate <- function(o) {
    idx <- modelled$resi + o
    ok <- idx >= 1 & idx <= length(human_aa)
    mean(human_aa[idx[ok]] == modelled$aa[ok])
  }
  offsets <- -100:100
  rates <- map_dbl(offsets, match_rate)
  offset <- offsets[which.max(rates)]
  stopifnot(max(rates) > 0.95)
  modelled <- mutate(modelled, human_pos = resi + offset)
  # Human position -> alignment column
  col_of <- which(strsplit(human_gapped, "")[[1]] != "-")

  # Salmon isoforms (vdra.a.1, vdra.a.2, ...) are read separately and merged below
  seqs <- names(aln)[str_remove(names(aln), "(?<=\\.[ab])\\.\\d+$") %in% c(s$human, shown[[s$protein]])]
  contacts_of(s) %>%
    left_join(select(modelled, resi, human_pos), by = "resi") %>%
    mutate(column = col_of[human_pos]) %>%
    { stopifnot(!anyNA(.$column), all(substring(human, .$human_pos, .$human_pos) == .$aa)); . } %>%
    crossing(sequence = seqs) %>%
    mutate(residue = substring(aln[sequence], column, column)) %>%
    transmute(protein = s$protein, pdb = s$pdb, human_residue = paste0(aa, resi), human_pos,
              hbond, min_dist = round(min_dist, 2), sequence, residue)
}

contacts <- structures %>% split(seq_len(nrow(.))) %>% map(map_contacts) %>% list_rbind()

# Salmon isoforms -> gene (vdra.a.1 -> vdra.a). Isoforms must agree at every contact residue
contacts <- contacts %>%
  mutate(sequence = str_remove(sequence, "(?<=\\.[ab])\\.\\d+$")) %>%
  distinct()
stopifnot(!anyDuplicated(select(contacts, protein, human_residue, sequence)))
write_tsv(contacts, file.path(out_dir, "contacts.tsv"))

summary <- contacts %>%
  group_by(protein, sequence) %>%
  summarise(contacts = n(), identical = sum(residue == str_sub(human_residue, 1, 1)),
            hbond_contacts = sum(hbond), hbond_identical = sum(hbond & residue == str_sub(human_residue, 1, 1)),
            changes = paste0(human_residue[residue != str_sub(human_residue, 1, 1)], residue[residue != str_sub(human_residue, 1, 1)],
                             collapse = ", "),
            .groups = "drop") %>%
  mutate(order = match(sequence, c(structures$human, unlist(shown)))) %>%
  arrange(order) %>%
  select(-order)
write_tsv(summary, file.path(out_dir, "summary.tsv"))
print(as.data.frame(summary))

# Table: one sheet per protein, a column per sequence, then the count of residues identical to human
notes <- c(VDR = "human VDR ligand-binding domain with 1alpha,25-dihydroxyvitamin D3 (PDB 1DB1)",
           CYP2R1 = "human CYP2R1 with vitamin D3 (PDB 3C6G, chain A)",
           GC = "human vitamin D-binding protein (GC) with 25-hydroxyvitamin D3 (PDB 1J78, chain B)")
sheets <- map_chr(structures$protein, function(p) {
  x <- filter(contacts, protein == p)
  seq_order <- c(structures$human[structures$protein == p], shown[[p]])
  wide <- x %>%
    mutate(sequence = factor(sequence, levels = seq_order)) %>%
    arrange(human_pos, sequence) %>%
    select(human_residue, hbond, min_dist, sequence, residue) %>%
    pivot_wider(names_from = sequence, values_from = residue)
  counts <- x %>%
    group_by(sequence) %>%
    summarise(n = paste0(sum(residue == str_sub(human_residue, 1, 1)), "/", n())) %>%
    { set_names(.$n, .$sequence)[seq_order] }
  note <- paste0("Residues with a heavy atom within 4 A of the ligand in the structure of the ", notes[[p]],
                 ", numbered as in the PDB entry. Hydrogen bond: a protein N or O atom within 3.5 A of a ligand O ",
                 "atom. Residues of the other proteins are taken from the alignments of the phylogenetic analysis.")
  out <- wide %>%
    transmute(Residue = human_residue, `Hydrogen bond` = ifelse(hbond, "yes", ""),
              `Closest distance (A)` = sprintf("%.2f", min_dist), !!!select(wide, all_of(seq_order))) %>%
    rename_with(display_name, all_of(seq_order)) %>%
    add_row(Residue = "Identical to human") %>%
    mutate(across(everything(), ~ replace_na(as.character(.x), "")))
  out[nrow(out), display_name(seq_order)] <- as.list(counts)
  out <- add_row(out, Residue = c("", note)) %>% mutate(across(everything(), ~ replace_na(.x, "")))
  file <- file.path(out_dir, paste0(p, ".tsv"))
  write_tsv(out, file)
  file
})
stopifnot(system2("python3", c(shQuote(file.path(proj_dir, "scripts", "tsv_to_xlsx.py")), shQuote(sheets),
                              shQuote(table_file))) == 0)
