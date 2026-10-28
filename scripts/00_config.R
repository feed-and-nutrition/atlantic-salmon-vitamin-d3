# Shared paths and helpers. Every analysis script sources this first.
# Override locations with environment variables when running elsewhere:
#   VD3_ROOT     project root (default: working directory, or its parent when run from scripts/)
#   VD3_TXI_DIR  directory holding txi_gene.rds of the atlas from tximport (only 00_1, and only without
#                salmon_quant/PRJNA72713/tpm_gene.tsv.gz from 05_1)
#   VD3_FIG_DIR  where run_all.R copies the figures and the scripts write the supplementary tables
#                (default: manuscript/v6)
# Scripts start with:
#   source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# Map R's generic "sans"/"serif" to the fonts used for the v2 figures (Arimo, Nimbus Roman),
# so text widths and layout do not depend on the machine's fontconfig defaults
register_generic_font <- function(generic, family) {
  path <- function(weight, italic) systemfonts::match_fonts(family, weight = weight, italic = italic)$path
  if (!grepl(family, basename(path("normal", FALSE)), ignore.case = TRUE) &&
      !grepl(str_remove_all(family, " "), path("normal", FALSE), ignore.case = TRUE)) {
    warning("Font '", family, "' not found. Figures will use the system default ", generic, " font.")
    return(invisible(FALSE))
  }
  systemfonts::register_font(generic, plain = path("normal", FALSE), bold = path("bold", FALSE),
                             italic = path("normal", TRUE), bolditalic = path("bold", TRUE))
}
register_generic_font("sans", "Arimo")
register_generic_font("serif", "Nimbus Roman")
register_generic_font("", "Arimo")  # ggplot2 themes default to family "", which ragg would resolve to the fontconfig default

proj_dir <- Sys.getenv("VD3_ROOT", unset = if (file.exists("00_config.R")) normalizePath("..") else getwd())
stopifnot(file.exists(file.path(proj_dir, "scripts", "00_config.R")))

msa_dir   <- file.path(proj_dir, "msa")
hk_dir    <- file.path(proj_dir, "hk_genes")
vd3_dir   <- file.path(proj_dir, "vd3_genes")
quant_dir <- file.path(proj_dir, "salmon_quant", "PRJNA72713")
fig_dir   <- Sys.getenv("VD3_FIG_DIR", unset = file.path(proj_dir, "manuscript", "v6"))  # current manuscript version
txi_dir   <- Sys.getenv("VD3_TXI_DIR", unset = NA)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Genome-wide gene-level TPM (one column per tissue), exported by 00_1 from txi_gene.rds
tpm_all_file <- file.path(quant_dir, "tpm_atlantic_salmon_15_tissues.tsv.gz")

expected_tissues <- c("brain", "eye", "gill", "gut", "head_kidney", "heart", "kidney", "liver",
                      "muscle", "nose", "ovary", "pyloric_caecum", "skin", "spleen", "testis")

# The PRJNA72713 heart sample is mostly liver tissue (Fig S2, 05_4), so expression analyses use 14 tissues.
# Raw exports (00_1) and the QC figure (Fig S1) keep all 15.
excluded_tissues <- "heart"
analysis_tissues <- setdiff(expected_tissues, excluded_tissues)
# Samples whose profile does not match their label (05_4, salmon_quant/tissue_identity.tsv), as study:sample
excluded_samples <- c("PRJNA72713:heart" = "mostly liver",
                      "PRJNA543940:heart_h59" = "pyloric caecum profile",
                      "PRJNA543940:heart_l63" = "pyloric caecum profile")
# Text-search candidates that are not VD3 genes, left out of the expression analyses (00_1, 05_2). cyp2r1.c
# (LOC106586122) is named "vitamin D(3) 25-hydroxylase" in RefSeq but is the ortholog of zebrafish cyp2x12
# (si:dkey-286h2.3), a fish CYP2 subfamily without a mammalian ortholog. It stays in the cyp2r1 tree (Fig 2c).
excluded_genes <- c(cyp2r1.c = "LOC106586122")

# Replicate studies (Fig 4, Fig S5, Table S10), studies with the most tissues first. Not the PRJNA72713 atlas (Fig 3)
# or PRJNA260929 (one sample per tissue, and its heart sample does not match heart in other studies).
replicate_studies <- c("PRJNA359635", "PRJEB47409", "ROS_1", "PRJNA647285", "PRJNA543940", "PRJNA380580", "PRJNA626461")

# Tissue order of Figs 3 and 4, grouped by organ system (atlas and replicate-study tissues)
tissue_order <- c("brain", "eye", "nose", "gill", "skin", "muscle", "bone", "blood", "heart", "liver", "gut",
                  "distal_intestine", "pyloric_caecum", "spleen", "head_kidney", "kidney", "ovary", "testis")

# Gene order of Figs 3 and 4: gene families in pathway order (first appearance), Ss4R copies .a before .b
pathway_order <- function(genes) {
  family <- str_remove(genes, "\\.[ab]$")
  genes[order(match(family, unique(family)), genes)]
}

# "head_kidney" / "Head_kidney" -> "Head kidney"
tissue_label <- function(x) str_to_sentence(str_replace_all(x, "_", " "))

# Stop early if a TPM matrix does not hold exactly the expected tissues
check_tissues <- function(tissues) {
  got <- str_to_lower(tissues)
  missing <- setdiff(expected_tissues, got)
  extra <- setdiff(got, expected_tissues)
  if (length(missing) || length(extra)) {
    stop("Unexpected tissue set. Missing: ", paste(missing, collapse = ", "),
         ". Extra: ", paste(extra, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Read a per-gene TPM table (Name, Gene_id, ..., one column per tissue) and derive
# log2(TPM + 1), gene-wise Z-scores, and the SD of log2(TPM + 1) per gene. Excluded tissues are dropped.
read_gene_tpm <- function(file) {
  df <- read_tsv(file, col_types = cols(Gene_id = col_character()))
  tissue_cols <- setdiff(colnames(df), c("Name", "Gene_id", "Type", "Gene_symbol", "Gene_name"))
  check_tissues(tissue_cols)
  tissue_cols <- tissue_cols[!(str_to_lower(tissue_cols) %in% excluded_tissues)]
  tpm <- as.matrix(df[, tissue_cols])
  rownames(tpm) <- df$Name
  colnames(tpm) <- tissue_label(tissue_cols)
  log_tpm <- log2(tpm + 1)
  list(info = df[, setdiff(colnames(df), tissue_cols)],
       tpm = tpm,
       log_tpm = log_tpm,
       z = t(scale(t(log_tpm))),
       sd = apply(log_tpm, 1, sd))
}

# Long format with the "expressed" flag used for the white x markers (TPM > 1)
to_long <- function(log_tpm) {
  as_tibble(log_tpm, rownames = "gene") %>%
    pivot_longer(-gene, names_to = "tissue", values_to = "logTPM") %>%
    mutate(expressed = logTPM > log2(1 + 1))
}

stability_annotation <- function(sd_values) {
  data.frame(SD = sd_values, Stability = ifelse(sd_values < 1, "Stable", "Variable"),
             row.names = names(sd_values))
}
stability_colors <- list(Stability = c(Stable = "#1B9E77", Variable = "#D95F02"))

# Tissue specificity index tau (Yanai et al. 2005) of log2(TPM + 1): 0 = even, 1 = one tissue
tau_index <- function(x) if (max(x) == 0) NA_real_ else sum(1 - x / max(x)) / (length(x) - 1)

# Table 2 pattern: restricted if the level varies strongly (SD >= 1) or is concentrated in few
# tissues (tau >= 0.8, Kryuchkova-Mostacci and Robinson-Rechavi 2017). The SD alone misses weakly
# expressed tissue-specific genes such as cyp27b1. Not "tissue enriched", which has a stricter meaning in the
# Human Protein Atlas (at least four-fold above all other tissues). Genes below 1 TPM in every tissue are not
# classified, because SD and tau are unstable near background (cyp2r1.b).
expr_pattern <- function(sd, tau, peak_tpm = Inf, sd_cut = 1, tau_cut = 0.8) {
  ifelse(peak_tpm < 1, "Not classified",
         ifelse(sd >= sd_cut | (!is.na(tau) & tau >= tau_cut), "Restricted", "Broad"))
}
# Fig 3 and Fig 4 panels: expression level by the atlas maximum of log2(TPM + 1), high > 7, medium > 3.5,
# low otherwise. Returns a list of gene names per level.
expr_level_groups <- function(log_tpm) {
  max_log <- apply(log_tpm, 1, max)
  list(High = names(max_log)[max_log > 7],
       Medium = names(max_log)[max_log > 3.5 & max_log <= 7],
       Low = names(max_log)[max_log <= 3.5])
}
# Fig S6 sidebars: the two criteria and the resulting pattern, as fixed categories so every panel
# uses the same colours. 04_2 lists "Not classified" only in the panel where it occurs (it does not fit the others)
pattern_annotation <- function(sd_values, tau_values, peak_tpm) {
  pattern <- expr_pattern(sd_values, tau_values, peak_tpm)
  data.frame(SD = factor(ifelse(sd_values >= 1, "\u2265 1", "< 1"), levels = c("\u2265 1", "< 1")),
             tau = factor(ifelse(!is.na(tau_values) & tau_values >= 0.8, "\u2265 0.8", "< 0.8"),
                          levels = c("\u2265 0.8", "< 0.8")),
             Pattern = factor(pattern, levels = c("Restricted", "Broad", intersect("Not classified", pattern))),
             row.names = names(sd_values))
}
pattern_colors <- list(SD = c("#4D4D4D", "#E0E0E0"),
                       tau = c("#7570B3", "#DADAEB"),
                       Pattern = c(Restricted = "#D95F02", Broad = "#1B9E77", `Not classified` = "#BDBDBD"))
names(pattern_colors$SD) <- c("\u2265 1", "< 1")
names(pattern_colors$tau) <- c("\u2265 0.8", "< 0.8")

# Z-score heatmap with sidebars: SD/Stability by default (Fig S3b), SD/tau/Pattern for Fig S6
z_heatmap <- function(z, sd_values, title, annotation = stability_annotation(sd_values[rownames(z)]),
                      annotation_colors = stability_colors, ...) {
  pheatmap::pheatmap(
    z,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_col = 0,
    annotation_row = annotation,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks = seq(-3, 3, length.out = 101),
    angle_col = "45",
    fontsize = 8,
    fontsize_row = 9,
    fontsize_col = 9,
    main = title,
    silent = TRUE,
    ...
  )
}

# Row and column order of a clustered pheatmap, so the absolute-expression panels
# use the same layout as the Z-score panels
ph_order <- function(ph, z) {
  list(genes = rownames(z)[ph$tree_row$order], tissues = colnames(z)[ph$tree_col$order])
}

# Expression heatmap (log2(TPM + 1)) with white x for TPM <= 1
expr_tile_plot <- function(long, genes, tissues, title, key_height = NULL) {
  long <- long %>%
    filter(gene %in% genes) %>%
    mutate(tissue = factor(tissue, levels = tissues),
           gene = fct_rev(factor(gene, levels = genes)))
  p <- ggplot(long, aes(tissue, gene, fill = logTPM)) +
    geom_tile() +
    geom_point(data = filter(long, !expressed), color = "white", shape = 4, size = 1.5) +
    scale_fill_viridis_c(name = "log2(TPM + 1)") +
    labs(title = title, x = "", y = "") +
    scale_y_discrete(position = "right") +
    theme_grey(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(size = 11, hjust = 0.5),
          plot.title.position = "plot",  # centre over the whole panel so long titles clear the a/b/c labels
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          plot.margin = margin(5.5, 5.5, 5.5, 45))  # room for a long first rotated tissue label
  if (!is.null(key_height)) p <- p + theme(legend.key.height = unit(key_height, "cm"))
  p
}

# Save a figure as SVG and as 600 dpi PNG (180 mm wide = 4252 px)
save_fig <- function(plot, path_no_ext, width, height) {
  ggsave2(paste0(path_no_ext, ".svg"), plot = plot, width = width, height = height,
          units = "mm", fix_text_size = FALSE)
  ggsave2(paste0(path_no_ext, ".png"), plot = plot, width = width, height = height,
          units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
  message("Saved ", path_no_ext, ".{svg,png}")
}
