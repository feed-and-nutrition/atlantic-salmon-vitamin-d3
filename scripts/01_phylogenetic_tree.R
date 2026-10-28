# Fig 2 (cyp27a1 with cyp27a2/a3, gc, cyp2r1, vdr, one protein per gene) and Fig S4 (cyp27b1, cyp24a1, and the
# dhcr7 and vdr trees with all salmon isoforms)
# Input:  msa/<tree>/<tree>.afa.treefile (IQ-TREE 2, -m MFP -B 1000). The gc, vdr_extended and cyp27a1 trees are built
#         by 01_0_build_trees.R
# Output: msa/tree2.{svg,png} (Fig 2), msa/tree1.{svg,png} (Fig S4)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
})

if (!exists("is.waive")) {
  is.waive <- function(x) inherits(x, "waiver")
}

# Root on the first mammalian tip (human or mouse), or on the outgroup clade given in outgroups. edgelabel = TRUE
# keeps each support value on its branch (without it, root() leaves the labels on the old nodes)
# (cyp27a1: the CYP27A family, with CYP27B1 as outgroup)
outgroups <- list(gc = c("Human GC", "Mouse Gc", "Chicken GC", "Gar gc", "Zebrafish gc"), cyp27a1 = "Human CYP27B1",
                  vdr_extended = c("Human VDR", "Mouse Vdr", "Chicken VDR"))
# Tips shown under another name: cyp2r1.c is not in the gene set (excluded_genes in 00_config.R) and its
# zebrafish ortholog si:dkey-286h2.3 (":" is not allowed in Newick) is cyp2x12
tip_names <- c("Salmon cyp2r1.c" = "Salmon LOC106586122", "Zebrafish si_dkey-286h2.3" = "Zebrafish cyp2x12")
load_tree <- function(gene) {
  tree <- read.tree(file.path(msa_dir, gene, paste0(gene, ".afa.treefile")))
  tree$tip.label <- str_replace(tree$tip.label, "_", " ")
  tree$tip.label <- coalesce(tip_names[tree$tip.label], tree$tip.label)
  if (!is.null(outgroups[[gene]])) {
    stopifnot(all(outgroups[[gene]] %in% tree$tip.label), is.monophyletic(tree, outgroups[[gene]]))
    return(root(tree, outgroup = outgroups[[gene]], resolve.root = TRUE, edgelabel = TRUE))
  }
  mammal_id <- grep("Homo|Human|Mus|Mouse", tree$tip.label, value = TRUE)
  if (length(mammal_id) > 0) root(tree, outgroup = mammal_id[1], resolve.root = TRUE, edgelabel = TRUE) else tree
}

# gc (Fig 2b) is drawn as a cladogram: its salmonid branches are too short to show next to the
# GC/albumin split, and the panel is about topology. The other trees are phylograms.
cladograms <- "gc"
panel_titles <- c(gc = "gc and albumins", cyp27a1 = "cyp27a1, cyp27a2 and cyp27a3",
                  vdr_extended = "vdr", dhcr7 = "dhcr7, all salmon isoforms", vdr = "vdr, all salmon isoforms")

create_tree <- function(gene, xlim) {
  clado <- gene %in% cladograms
  p <- ggtree(load_tree(gene), branch.length = if (clado) "none" else "branch.length") +
    geom_tree(size = 0.8) +
    geom_tiplab(size = 3, fontface = "italic", family = "serif",
                align = TRUE, linesize = 0.5, linetype = "dotted", offset = if (clado) 0.3 else 0.05) +
    # UFBoot support, shown only above 70
    geom_nodelab(aes(x = branch, label = label,
                     subset = !isTip & (suppressWarnings(as.numeric(label)) > 70)),
                 geom = "label", alpha = 0.8, size = 2.5, family = "sans") +
    ggtitle(if (gene %in% names(panel_titles)) panel_titles[[gene]] else gene) +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    xlim(0, xlim)
  if (clado) p + theme_tree() + theme(plot.title = element_text(hjust = 0.5, size = 11))
  else p + geom_treescale(x = 0, y = 0, width = 0.1, fontsize = 2.5) + theme_tree2() +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
}

# x-axis limits chosen so tip labels are not cut off (gc: number of nodes, as a cladogram)
xlims <- c(cyp27a1 = 1.9, gc = 15, cyp2r1 = 3, vdr_extended = 1.45, cyp27b1 = 2, cyp24a1 = 1, dhcr7 = 0.8, vdr = 0.75)
p <- imap(xlims, ~ create_tree(.y, .x))

# cyp27b1, cyp24a1 (5 tips each) and the dhcr7 isoform tree (7 tips) are stacked on the left, the vdr isoform tree
# (17 tips) takes the right column
p_tree1 <- plot_grid(plot_grid(p$cyp27b1, p$cyp24a1, p$dhcr7, ncol = 1, rel_heights = c(1, 1, 1.3),
                               labels = c("a", "b", "c"), label_size = 11),
                     plot_grid(p$vdr, labels = "d", label_size = 11),
                     ncol = 2)
# Columns are stacked separately so each tree gets a height that suits its number of tips
p_tree2 <- plot_grid(plot_grid(p$cyp27a1, p$cyp2r1, ncol = 1, rel_heights = c(1.5, 1), labels = c("a", "c"),
                               label_size = 11),
                     plot_grid(p$gc, p$vdr_extended, ncol = 1, rel_heights = c(1.45, 1.1), labels = c("b", "d"),
                               label_size = 11),
                     ncol = 2)

save_fig(p_tree1, file.path(msa_dir, "tree1"), 180, 170)
save_fig(p_tree2, file.path(msa_dir, "tree2"), 180, 190)
