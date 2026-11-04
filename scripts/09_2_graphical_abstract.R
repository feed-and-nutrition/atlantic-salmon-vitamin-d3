# Graphical abstract (required by CBP Part D): the main results on a drawing of an Atlantic salmon, with callouts in
# the colours of Figure 1 (09_1). Drawn in code because Elsevier does not allow generative AI images in graphical
# abstracts. Canvas 150 x 60 mm, y downwards. Elsevier asks for a 5:2 image of at least 1328 x 531 px and
# 300 dpi that stays legible when shown 200 px high, so no text is smaller than 10 pt. The heatmap icon is the
# atlas log2(TPM + 1) of Figure 3.
# Input:  vd3_genes/tpm_vd3_genes.tsv (00_1)
# Output: d3pathway/graphical_abstract.{svg,png}, manuscript/v10/graphical_abstract.tiff (600 dpi, LZW)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

W <- 150
H <- 60
lw <- 0.3
col <- c(gene = "#d5e8d4", hypothesis = "#fff2cc", region = "#ea6b66", outline = "grey30", body = "#e3eaf1",
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
cubic <- function(s, n = 16) {
  t <- seq(0, 1, length.out = n)
  tibble(x = (1 - t)^3 * s[1] + 3 * (1 - t)^2 * t * s[3] + 3 * (1 - t) * t^2 * s[5] + t^3 * s[7],
         y = (1 - t)^3 * s[2] + 3 * (1 - t)^2 * t * s[4] + 3 * (1 - t) * t^2 * s[6] + t^3 * s[8])
}
ellipse <- function(cx, cy, rx, ry, n = 40) {
  a <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + rx * cos(a), y = cy + ry * sin(a))
}

# Atlantic salmon in body lengths (snout at 0, tip of the caudal fin at 1, y up), head to the left. The fins lie
# behind the body except the pectoral fin. All colours are flat.
seg <- function(...) map_dfr(list(...), cubic)
back <- seg(
  c(0.000, 0.000, 0.004, 0.030, 0.060, 0.075, 0.180, 0.088),      # snout to nape
  c(0.180, 0.088, 0.260, 0.100, 0.340, 0.114, 0.420, 0.116),      # back
  c(0.420, 0.116, 0.460, 0.117, 0.500, 0.113, 0.530, 0.108),
  c(0.530, 0.108, 0.600, 0.100, 0.700, 0.082, 0.760, 0.070),
  c(0.760, 0.070, 0.765, 0.085, 0.778, 0.093, 0.792, 0.092),      # adipose fin
  c(0.792, 0.092, 0.803, 0.091, 0.808, 0.080, 0.812, 0.062),
  c(0.812, 0.062, 0.835, 0.052, 0.865, 0.045, 0.885, 0.042))      # caudal peduncle
belly <- seg(
  c(0.885, -0.042, 0.860, -0.045, 0.835, -0.050, 0.810, -0.056),
  c(0.810, -0.056, 0.785, -0.062, 0.760, -0.070, 0.735, -0.077),
  c(0.735, -0.077, 0.680, -0.090, 0.620, -0.097, 0.565, -0.100),  # belly
  c(0.565, -0.100, 0.450, -0.108, 0.300, -0.100, 0.200, -0.083),
  c(0.200, -0.083, 0.120, -0.066, 0.020, -0.035, 0.000, 0.000))   # throat and chin
# The body is closed across the base of the caudal fin, which is left without an outline
body <- bind_rows(back, seg(c(0.885, 0.042, 0.893, 0.020, 0.893, -0.020, 0.885, -0.042)), belly)
body_outline <- bind_rows(belly, back)
# Fins as closed outlines (base inside the body), each with the margin that its rays end on
dorsal_margin <- seg(c(0.458, 0.215, 0.480, 0.185, 0.520, 0.150, 0.552, 0.130))
caudal_margin <- seg(c(1.000, 0.132, 0.985, 0.090, 0.968, 0.040, 0.962, 0.000),
                     c(0.962, 0.000, 0.968, -0.040, 0.985, -0.090, 1.000, -0.132))
anal_margin <- seg(c(0.760, -0.130, 0.785, -0.115, 0.808, -0.095, 0.822, -0.083))
fins <- bind_rows(
  dorsal = bind_rows(seg(c(0.415, 0.105, 0.420, 0.140, 0.440, 0.190, 0.458, 0.215)), dorsal_margin,
                     seg(c(0.552, 0.130, 0.548, 0.120, 0.544, 0.110, 0.540, 0.104))),
  caudal = bind_rows(seg(c(0.870, 0.040, 0.910, 0.060, 0.960, 0.100, 1.000, 0.132)), caudal_margin,
                     seg(c(1.000, -0.132, 0.960, -0.100, 0.910, -0.060, 0.870, -0.040))),
  anal = bind_rows(seg(c(0.738, -0.070, 0.742, -0.095, 0.750, -0.115, 0.760, -0.130)), anal_margin,
                   seg(c(0.822, -0.083, 0.818, -0.072, 0.815, -0.062, 0.812, -0.052))),
  pelvic = seg(c(0.505, -0.098, 0.530, -0.120, 0.560, -0.138, 0.588, -0.146),
               c(0.588, -0.146, 0.575, -0.125, 0.568, -0.110, 0.560, -0.096)),
  .id = "fin") %>%
  mutate(dark = fin %in% c("dorsal", "caudal"))
pectoral_margin <- seg(c(0.312, -0.094, 0.296, -0.102, 0.262, -0.098, 0.224, -0.080))
pectoral <- bind_rows(seg(c(0.205, -0.046, 0.240, -0.056, 0.280, -0.075, 0.312, -0.094)), pectoral_margin)
# Rays from evenly spaced points on the fin base to evenly spaced points on its margin
rays <- function(base0, base1, margin, n) {
  m <- margin[round(seq(1, nrow(margin), length.out = n)), ]
  tibble(x = seq(base0[1], base1[1], length.out = n), y = seq(base0[2], base1[2], length.out = n),
         xend = m$x, yend = m$y)
}
fin_rays <- bind_rows(rays(c(0.420, 0.106), c(0.540, 0.104), dorsal_margin, 9),
                      rays(c(0.872, 0.034), c(0.872, -0.034), caudal_margin, 13),
                      rays(c(0.742, -0.072), c(0.812, -0.054), anal_margin, 7))
pectoral_rays <- rays(c(0.207, -0.050), c(0.214, -0.072), pectoral_margin, 6)
fx <- 48.5   # snout position and body length on the canvas
fy <- 35
L <- 67
place <- function(d) mutate(d, across(any_of(c("x", "xend")), ~ fx + .x * L), across(any_of(c("y", "yend")), ~ fy - .x * L))
at <- function(x, y) c(fx + x * L, fy - y * L)
body <- place(body)
body_outline <- place(body_outline)
fins <- place(fins)
pectoral <- place(pectoral)
fin_rays <- place(fin_rays)
pectoral_rays <- place(pectoral_rays)
lateral_line <- place(cubic(c(0.205, 0.022, 0.450, 0.018, 0.700, 0.010, 0.885, 0.004)))
operculum <- place(cubic(c(0.185, 0.072, 0.205, 0.035, 0.208, -0.015, 0.185, -0.072)))
mouth <- place(cubic(c(0.003, -0.006, 0.030, -0.012, 0.055, -0.016, 0.075, -0.020)))
eye <- place(bind_rows(iris = ellipse(0.055, 0.030, 0.013, 0.013), pupil = ellipse(0.055, 0.030, 0.0085, 0.0085),
                       glint = ellipse(0.059, 0.034, 0.0028, 0.0028), .id = "part"))
brain <- place(ellipse(0.100, 0.050, 0.027, 0.016))
liver <- place(ellipse(0.275, -0.040, 0.062, 0.024))
kidney <- place(seg(c(0.170, 0.039, 0.170, 0.048, 0.176, 0.055, 0.185, 0.055),   # head kidney
                    c(0.185, 0.055, 0.210, 0.056, 0.235, 0.057, 0.250, 0.056),
                    c(0.250, 0.056, 0.270, 0.055, 0.290, 0.048, 0.320, 0.047),
                    c(0.320, 0.047, 0.450, 0.047, 0.620, 0.042, 0.690, 0.036),   # trunk kidney
                    c(0.690, 0.036, 0.700, 0.033, 0.680, 0.027, 0.660, 0.027),
                    c(0.660, 0.027, 0.550, 0.027, 0.400, 0.031, 0.320, 0.031),
                    c(0.320, 0.031, 0.290, 0.031, 0.270, 0.023, 0.250, 0.023),
                    c(0.250, 0.023, 0.230, 0.022, 0.205, 0.023, 0.185, 0.024),
                    c(0.185, 0.024, 0.176, 0.025, 0.170, 0.031, 0.170, 0.039)))

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
brain_x <- at(0.100, 0.050)
kidney_x <- at(0.600, 0.034)
liver_x <- at(0.275, -0.040)
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
  geom_polygon(data = filter(fins, dark), aes(x, y, group = fin), fill = col[["fin"]], colour = col[["outline"]],
               linewidth = lw) +
  geom_polygon(data = filter(fins, !dark), aes(x, y, group = fin), fill = col[["fin_light"]],
               colour = col[["outline"]], linewidth = lw) +
  geom_segment(data = fin_rays, aes(x, y, xend = xend, yend = yend), colour = col[["ray"]], linewidth = lw * 0.4) +
  geom_polygon(data = body, aes(x, y), fill = col[["body"]]) +
  geom_path(data = body_outline, aes(x, y), colour = col[["outline"]], linewidth = lw * 1.5, linejoin = "round") +
  geom_path(data = lateral_line, aes(x, y), colour = "grey50", linewidth = lw * 0.6) +
  geom_path(data = operculum, aes(x, y), colour = col[["outline"]], linewidth = lw) +
  geom_path(data = mouth, aes(x, y), colour = col[["outline"]], linewidth = lw) +
  geom_polygon(data = pectoral, aes(x, y), fill = col[["fin_light"]], colour = col[["outline"]], linewidth = lw) +
  geom_segment(data = pectoral_rays, aes(x, y, xend = xend, yend = yend), colour = col[["ray"]],
               linewidth = lw * 0.4) +
  geom_polygon(data = brain, aes(x, y), fill = col[["brain"]], colour = col[["outline"]], linewidth = lw * 0.6) +
  geom_polygon(data = liver, aes(x, y), fill = col[["liver"]], colour = col[["outline"]], linewidth = lw * 0.6) +
  geom_polygon(data = kidney, aes(x, y), fill = col[["kidney"]], colour = col[["outline"]], linewidth = lw * 0.6) +
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
