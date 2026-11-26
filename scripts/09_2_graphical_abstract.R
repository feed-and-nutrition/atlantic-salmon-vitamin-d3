# Graphical abstract (required by CBP Part D): the main results on a drawing of an Atlantic salmon, with callouts in
# the colours of Figure 1 (09_1). Drawn in code because Elsevier does not allow generative AI images in graphical
# abstracts. Canvas 150 x 60 mm, y downwards. Elsevier asks for a 5:2 image of at least 1328 x 531 px and
# 300 dpi that stays legible when shown 200 px high, so no text is smaller than 10 pt. The heatmap icon is the
# atlas log2(TPM + 1) of Figure 3.
# The salmon outline is traced from a public-domain USFWS drawing by 09_0.
# Input:  vd3_genes/tpm_vd3_genes.tsv (00_1), d3pathway/salmon_outline.tsv (09_0)
# Output: d3pathway/graphical_abstract.{svg,png}, manuscript/v11/graphical_abstract.tiff (600 dpi, LZW)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

W <- 150
H <- 60
lw <- 0.3
col <- c(gene = "#d5e8d4", hypothesis = "#fff2cc", region = "#ea6b66", outline = "grey30", body = "#e3eaf1",
         back = "#a9b8c8", spot = "grey15",
         fin = "#9fb2c4", fin_light = "#c9d4de", ray = "#6f8599",
         brain = "#f2b5b0", liver = "#a65a3c", kidney = "#7d2424", absent = "grey40", cross = "#d62728")
pt <- function(size) size / .pt  # font size in points -> geom_text size
text_mm <- function(x, size = 10, ...) systemfonts::string_width(x, family = "sans", size = size, res = 72, ...) * 25.4 / 72

# As in 09_1
rounded_rect <- function(x0, y0, x1, y1, r) {
  arc <- function(cx, cy, from) {
    a <- seq(from, from + pi / 2, length.out = 8)
    tibble(x = cx + r * cos(a), y = cy + r * sin(a))
  }
  bind_rows(arc(x1 - r, y1 - r, 0), arc(x0 + r, y1 - r, pi / 2), arc(x0 + r, y0 + r, pi), arc(x1 - r, y0 + r, 3 * pi / 2))
}
# Atlantic salmon traced from a public-domain drawing (09_0), in body lengths (snout at 0, tip of the caudal fin at
# 1, y up), head to the left. The fins lie behind the body except the pectoral fin. All colours are flat.
fish <- read_tsv(file.path(proj_dir, "d3pathway", "salmon_outline.tsv"), show_col_types = FALSE)
shape <- function(p) filter(fish, part == p)
fx <- 48.5   # snout position and body length on the canvas
fy <- 35
L <- 67
place <- function(d) mutate(d, x = fx + x * L, y = fy - y * L)
at <- function(x, y) c(fx + x * L, fy - y * L)
fish <- place(fish)
fins <- shape("fin") %>%
  mutate(dark = group %in% c("dorsal", "caudal", "adipose"))
eye <- filter(fish, str_starts(part, "eye_")) %>%
  mutate(part = str_remove(part, "eye_"))
# Callout anchors: the brain and liver centres, and the trunk kidney at 60% of the body length
brain_at <- colMeans(select(shape("brain"), x, y))
liver_at <- colMeans(select(shape("liver"), x, y))
kidney_at <- shape("kidney") %>%
  filter(abs(x - at(0.6, 0)[1]) < 1.5) %>%
  summarise(x = mean(x), y = mean(y)) %>%
  unlist()

# Left column icons: a tree with two pairs of tips, and the atlas heatmap (genes x tissues, as in Figure 3)
tree <- tribble(~x, ~y, ~xend, ~yend,
  3, 14, 6, 14,  6, 11, 6, 17,  6, 11, 9, 11,  6, 17, 9, 17,  9, 9.5, 9, 12.5,  9, 15.5, 9, 18.5,
  9, 9.5, 12, 9.5,  9, 12.5, 12, 12.5,  9, 15.5, 12, 15.5,  9, 18.5, 12, 18.5)
atlas <- read_tsv(file.path(vd3_dir, "tpm_vd3_genes.tsv"), show_col_types = FALSE) %>%
  pivot_longer(-c(Name, Gene_id, Type, Gene_symbol, Gene_name), names_to = "tissue", values_to = "tpm") %>%
  mutate(tissue = str_to_lower(tissue)) %>%
  filter(!tissue %in% excluded_tissues) %>%
  mutate(gene = match(Name, pathway_order(unique(Name))), tissue = match(tissue, intersect(tissue_order, tissue)),
         x = 2.6 + tissue * 0.75, y = 24.6 + gene * 0.75, value = log2(tpm + 1))

# Gene boxes, sized to their labels: the brain and kidney above the fish, the liver row below it
pad <- 1.5
gene_box <- function(x0, label, row) {
  tibble(x0, x1 = x0 + text_mm(label, italic = TRUE) + 2 * pad, label, row,
         y0 = if_else(row == "top", 7, 50.5), y1 = y0 + 5.5)
}
brain_x <- brain_at
kidney_x <- kidney_at
liver_x <- liver_at
genes <- bind_rows(
  gene_box(brain_x[1] - 10, "cyp27b1.a", "top"),
  gene_box(kidney_x[1] - 2, "cyp27b1", "top"),
  gene_box(32, "cyp27b1.b", "bottom"))
genes <- bind_rows(genes, gene_box(last(genes$x1) + 1.5, "alb1, alb2", "bottom")) %>%
  mutate(present = label != "cyp27b1")
gc_box <- tibble(x0 = last(genes$x1) + 1.5, x1 = x0 + text_mm("gc", italic = TRUE) + 2 * pad + 4, y0 = 50.5, y1 = 56)
box_poly <- function(d) pmap_dfr(d, function(x0, x1, y0, y1, ...) mutate(rounded_rect(x0, y0, x1, y1, 1.2), id = paste(x0, y0)))
organs <- tribble(~x, ~y, ~label,
  genes$x0[1], 4.7, "Brain",
  genes$x0[2], 4.7, "Kidney",
  32, 48.2, "Liver")
alb <- filter(genes, label == "alb1, alb2")
leaders <- tribble(~x, ~y, ~xend, ~yend,
  brain_x[1], 12.5, brain_x[1], brain_x[2],
  kidney_x[1], 16.8, kidney_x[1], kidney_x[2],
  (alb$x0 + alb$x1) / 2, 50.5, liver_x[1], liver_x[2])
no_gene_x <- gc_box$x1 + 1.5
hypotheses <- tribble(~x0, ~y0, ~x1, ~y1, ~label,
  121, 4, 148, 24, "Vitamin D₃\nactivation\noutside the\nkidney?",
  121, 40, 148, 56, "Carried by\nalbumins or\nlipoproteins?")
arrows <- tribble(~x, ~y, ~xend, ~yend,
  genes$x1[2] + 1.5, 9.75, 120.3, 9.75,
  no_gene_x + text_mm("no gene") + 1.5, 53.25, 120.3, 53.25)

p <- ggplot() +
  # Left column: what was studied
  geom_text(aes(2, 2, label = "Atlantic salmon"), hjust = 0, vjust = 1, fontface = "bold", size = pt(11),
            family = "sans") +
  geom_segment(data = tree, aes(x, y, xend = xend, yend = yend), linewidth = lw * 1.5, lineend = "square") +
  geom_text(aes(15, 11.5, label = "16 genes"), hjust = 0, fontface = "bold", size = pt(11), family = "sans") +
  geom_text(aes(15, 16.5, label = "7 paralog pairs"), hjust = 0, size = pt(10), family = "sans") +
  geom_tile(data = atlas, aes(x, y, fill = value), width = 0.75, height = 0.75) +
  scale_fill_viridis_c(guide = "none") +
  geom_text(aes(15, 28.5, label = "14-tissue atlas"), hjust = 0, size = pt(10), family = "sans") +
  geom_text(aes(15, 33.5, label = "7 RNA-seq studies"), hjust = 0, size = pt(10), family = "sans") +
  # Fish with organs
  geom_polygon(data = filter(fins, dark), aes(x, y, group = group), fill = col[["fin"]], colour = col[["outline"]],
               linewidth = lw) +
  geom_polygon(data = filter(fins, !dark), aes(x, y, group = group), fill = col[["fin_light"]],
               colour = col[["outline"]], linewidth = lw) +
  geom_path(data = shape("ray"), aes(x, y, group = group), colour = col[["ray"]], linewidth = lw * 0.4) +
  geom_polygon(data = shape("body"), aes(x, y), fill = col[["body"]]) +
  geom_polygon(data = shape("back"), aes(x, y), fill = col[["back"]]) +
  geom_polygon(data = shape("spot"), aes(x, y, group = group), fill = col[["spot"]]) +
  geom_path(data = shape("outline"), aes(x, y), colour = col[["outline"]], linewidth = lw * 1.5, linejoin = "round") +
  geom_path(data = shape("lateral_line"), aes(x, y), colour = "grey50", linewidth = lw * 0.6) +
  geom_path(data = shape("operculum"), aes(x, y), colour = col[["outline"]], linewidth = lw) +
  geom_path(data = shape("preopercle"), aes(x, y), colour = col[["outline"]], linewidth = lw * 0.6) +
  geom_path(data = shape("maxilla"), aes(x, y), colour = col[["outline"]], linewidth = lw) +
  geom_polygon(data = shape("pectoral"), aes(x, y), fill = col[["fin_light"]], colour = col[["outline"]],
               linewidth = lw) +
  geom_path(data = shape("pectoral_ray"), aes(x, y, group = group), colour = col[["ray"]], linewidth = lw * 0.4) +
  geom_polygon(data = shape("brain"), aes(x, y), fill = col[["brain"]], colour = col[["outline"]],
               linewidth = lw * 0.6) +
  geom_polygon(data = shape("liver"), aes(x, y), fill = col[["liver"]], colour = col[["outline"]],
               linewidth = lw * 0.6) +
  geom_polygon(data = shape("kidney"), aes(x, y), fill = col[["kidney"]], colour = col[["outline"]],
               linewidth = lw * 0.6) +
  geom_polygon(data = filter(eye, part == "iris"), aes(x, y), fill = "#d9cda4") +
  geom_polygon(data = filter(eye, part == "pupil"), aes(x, y), fill = "grey10") +
  geom_polygon(data = filter(eye, part == "glint"), aes(x, y), fill = "white") +
  # Callouts: organ, genes expressed (green), barely expressed (dashed), absent (crossed)
  geom_segment(data = leaders, aes(x, y, xend = xend, yend = yend), linewidth = lw) +
  geom_point(data = leaders, aes(xend, yend), size = 0.9) +
  geom_text(data = organs, aes(x, y, label = label), hjust = 0, colour = col[["region"]], fontface = "bold",
            size = pt(10.5), family = "sans") +
  geom_polygon(data = box_poly(filter(genes, present)), aes(x, y, group = id), fill = col[["gene"]],
               colour = "black", linewidth = lw) +
  geom_polygon(data = box_poly(bind_rows(filter(genes, !present), gc_box)), aes(x, y, group = id), fill = "white",
               colour = col[["absent"]], linewidth = lw, linetype = "22") +
  geom_text(data = genes, aes((x0 + x1) / 2, (y0 + y1) / 2, label = label, colour = if_else(present, "black", "grey25")),
            fontface = "italic", size = pt(10), family = "sans") +
  scale_colour_identity() +
  geom_text(aes(gc_box$x0 + pad, 53.25, label = "gc"), hjust = 0, colour = "grey25", fontface = "italic",
            size = pt(10), family = "sans") +
  geom_segment(aes(x = gc_box$x1 - pad - 3, y = 51.75, xend = gc_box$x1 - pad, yend = 54.75), colour = col[["cross"]],
               linewidth = lw * 2) +
  geom_segment(aes(x = gc_box$x1 - pad - 3, y = 54.75, xend = gc_box$x1 - pad, yend = 51.75), colour = col[["cross"]],
               linewidth = lw * 2) +
  geom_text(aes(no_gene_x, 53.25, label = "no gene"), hjust = 0, size = pt(10), family = "sans") +
  geom_text(aes(genes$x0[2], 15, label = "barely expressed"), hjust = 0, size = pt(10), family = "sans") +
  # Hypotheses
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), linewidth = lw * 1.5,
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_polygon(data = box_poly(hypotheses), aes(x, y, group = id), fill = col[["hypothesis"]], colour = "black",
               linewidth = lw) +
  geom_text(data = hypotheses, aes((x0 + x1) / 2, (y0 + y1) / 2, label = label), fontface = "bold",
            size = pt(10), lineheight = 0.95, family = "sans") +
  scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
  scale_y_reverse(limits = c(H, 0), expand = c(0, 0)) +
  coord_fixed() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA))

save_fig(p, file.path(proj_dir, "d3pathway", "graphical_abstract"), W, H)
ggsave(file.path(fig_dir, "graphical_abstract.tiff"), p, width = W, height = H, units = "mm", dpi = 600,
       bg = "white", device = ragg::agg_tiff, compression = "lzw")
message("Wrote ", file.path(fig_dir, "graphical_abstract.tiff"))
