# Fig 1: the canonical vitamin D3 pathway (human model), redrawn in the style of the draw.io figure of v2 to v5
# (d3pathway/vd3_pathwa_colors.svg). Previtamin D3 and its thermal isomerization to vitamin D3 lie in the skin,
# dietary vitamin D is absorbed with lipids (not carried by DBP from the gut), and CYP24A1 acts on both 25(OH)D3 and
# 1,25(OH)2D3. Coordinates are those of a 640 x 376.5 canvas (170 x 100 mm), y downwards, 10 units = 10 px font.
# Output: d3pathway/pathway.{svg,png} (Fig 1)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

W <- 640
H <- W / 1.7
u <- 170 / W                     # mm per canvas unit
lw <- 0.36                       # linewidth of a 1-unit line
col <- c(metabolite = "#fff2cc", enzyme = "#d5e8d4", sky = "#d4e1f5", region = "#ea6b66", region_fill = "#fffafa")

rounded_rect <- function(x0, y0, x1, y1, r) {
  arc <- function(cx, cy, from) {
    a <- seq(from, from + pi / 2, length.out = 8)
    tibble(x = cx + r * cos(a), y = cy + r * sin(a))
  }
  if (r == 0) return(tibble(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1)))
  bind_rows(arc(x1 - r, y1 - r, 0), arc(x0 + r, y1 - r, pi / 2), arc(x0 + r, y0 + r, pi), arc(x1 - r, y0 + r, 3 * pi / 2))
}
cubic <- function(p0, p1, p2, p3, n = 12) {
  t <- seq(0, 1, length.out = n)
  tibble(x = (1 - t)^3 * p0[1] + 3 * (1 - t)^2 * t * p1[1] + 3 * (1 - t) * t^2 * p2[1] + t^3 * p3[1],
         y = (1 - t)^3 * p0[2] + 3 * (1 - t)^2 * t * p1[2] + 3 * (1 - t) * t^2 * p2[2] + t^3 * p3[2])
}
# Cloud and lightning outlines of the draw.io figure, scaled into a box
fit <- function(d, x0, y0, x1, y1) {
  mutate(d, x = x0 + (x - min(x)) / diff(range(x)) * (x1 - x0), y = y0 + (y - min(y)) / diff(range(y)) * (y1 - y0))
}
cloud_seg <- list(c(33.5, 75.5, 7.5, 75.5, 1, 93, 21.8, 96.5), c(21.8, 96.5, 1, 104.2, 24.4, 121, 41.3, 114),
                  c(41.3, 114, 53, 128, 92, 128, 105, 114), c(105, 114, 131, 114, 131, 100, 114.75, 93),
                  c(114.75, 93, 131, 79, 105, 65, 82.25, 72), c(82.25, 72, 66, 61.5, 40, 61.5, 33.5, 75.5))
cloud <- map_dfr(cloud_seg, ~ cubic(.x[1:2], .x[3:4], .x[5:6], .x[7:8]))
lightning <- tibble(x = c(275, 288.33, 277, 288.33, 277, 297.67, 315, 299, 308.33, 294.33, 303),
                    y = c(88, 71.25, 59.52, 49.47, 34.39, 8, 29.36, 41.09, 52.82, 62.87, 72.08))

# Tissue regions (double outline) and their labels
regions <- tribble(~name, ~x0, ~y0, ~x1, ~y1, ~lx, ~ly, ~hjust,
  "Skin",   318, 10,  628, 196, 312, 20,  1,
  "Liver",  12,  214, 256, 310, 14,  206, 0,
  "Kidney", 340, 214, 628, 366, 342, 206, 0)
# Metabolites (square corners) and enzymes or binding proteins (rounded, bold)
boxes <- tribble(~id, ~type, ~x0, ~y0, ~x1, ~y1, ~label,
  "chol",   "metabolite", 476, 20,  616, 48,  "Cholesterol",
  "dhcr7",  "enzyme",     540, 62,  616, 88,  "DHCR7",
  "dhc",    "metabolite", 476, 100, 616, 130, "7-Dehydrocholesterol\n(Provitamin D₃)",
  "pre",    "metabolite", 334, 100, 434, 130, "Previtamin D₃",
  "d3",     "metabolite", 334, 156, 434, 186, "Cholecalciferol\n(Vitamin D₃)",
  "dbp1",   "enzyme",     225, 180, 291, 204, "DBP / GC",
  "cyp2r1", "enzyme",     20,  228, 92,  254, "CYP2R1",
  "cyp27a1","enzyme",     20,  268, 92,  294, "CYP27A1",
  "oh_l",   "metabolite", 126, 266, 246, 296, "Calcidiol\n(25(OH)D₃)",
  "dbp2",   "enzyme",     262, 290, 328, 314, "DBP / GC",
  "oh_k",   "metabolite", 354, 230, 464, 260, "Calcidiol\n(25(OH)D₃)",
  "cyp27b1","enzyme",     384, 280, 452, 306, "CYP27B1",
  "cyp24a1","enzyme",     471, 280, 539, 306, "CYP24A1",
  "oh2",    "metabolite", 354, 326, 464, 356, "Calcitriol\n(1α,25(OH)₂D₃)",
  "oh24",   "metabolite", 520, 230, 620, 260, "24,25(OH)₂D₃\n(inactive)",
  "oh124",  "metabolite", 520, 326, 620, 356, "1,24,25(OH)₃D₃\n(inactive)",
  "vdr",    "enzyme",     100, 330, 170, 356, "VDR",
  "dbp3",   "enzyme",     230, 350, 296, 374, "DBP / GC")
# Edges as polylines; end: block = conversion, oval = catalysis, open = binding, none = translocation
edges <- list(
  list(c(510, 100, 510, 48), "block"),                       # 7-DHC -> cholesterol
  list(c(540, 75, 510, 75), "oval"),                         # DHCR7
  list(c(476, 115, 434, 115), "block"),                      # 7-DHC -> previtamin D3
  list(c(362, 64, 455, 64, 455, 115), "oval"),               # UV B
  list(c(384, 130, 384, 156), "block"),                      # thermal isomerization
  list(c(186, 93, 186, 171), "none"),                        # diet, absorbed with lipids
  list(c(334, 171, 186, 171), "none"),                       # vitamin D3 from the skin
  list(c(258, 180, 258, 171), "open"),                       # DBP binds vitamin D3
  list(c(186, 171, 186, 266), "block"),                      # 25-hydroxylation
  list(c(92, 241, 186, 241), "oval"),                        # CYP2R1
  list(c(92, 281, 112, 281, 112, 241), "none"),              # CYP27A1
  list(c(246, 276, 334, 276, 334, 245, 354, 245), "none"),   # 25(OH)D3 to the kidney
  list(c(295, 290, 295, 276), "open"),                       # DBP binds 25(OH)D3
  list(c(372, 260, 372, 326), "block"),                      # 1-alpha-hydroxylation
  list(c(384, 293, 372, 293), "oval"),                       # CYP27B1
  list(c(464, 245, 520, 245), "block"),                      # 24-hydroxylation of 25(OH)D3
  list(c(464, 341, 520, 341), "block"),                      # 24-hydroxylation of 1,25(OH)2D3
  list(c(505, 280, 505, 245), "oval"),                       # CYP24A1
  list(c(505, 306, 505, 341), "oval"),
  list(c(354, 341, 170, 341), "open"),                       # 1,25(OH)2D3 binds VDR in target cells
  list(c(263, 350, 263, 341), "open"))                       # DBP binds 1,25(OH)2D3
legend <- tribble(~y, ~end, ~label,
  46, "block", "Conversion",
  64, "oval",  "Catalysis",
  82, "open",  "Binding",
  100, "none", "Translocation")
edges <- c(edges, pmap(legend, function(y, end, label) list(c(20, y, 46, y), end)))

# Arrow heads at the last point of each edge, pointing along the last segment
head_geom <- function(pts, end) {
  n <- length(pts) / 2
  p <- matrix(pts, ncol = 2, byrow = TRUE)
  tip <- p[n, ]
  d <- (p[n, ] - p[n - 1, ]) / sqrt(sum((p[n, ] - p[n - 1, ])^2))
  nrm <- c(-d[2], d[1])
  switch(end,
    block = tibble(x = c(tip[1], tip[1] - 7 * d[1] + 3.5 * nrm[1], tip[1] - 7 * d[1] - 3.5 * nrm[1]),
                   y = c(tip[2], tip[2] - 7 * d[2] + 3.5 * nrm[2], tip[2] - 7 * d[2] - 3.5 * nrm[2])),
    open = tibble(x = c(tip[1] - 6 * d[1] + 3.5 * nrm[1], tip[1], tip[1] - 6 * d[1] - 3.5 * nrm[1]),
                  y = c(tip[2] - 6 * d[2] + 3.5 * nrm[2], tip[2], tip[2] - 6 * d[2] - 3.5 * nrm[2])),
    oval = tibble(x = tip[1] - 3 * d[1] + 3 * cos(seq(0, 2 * pi, length.out = 25)),
                  y = tip[2] - 3 * d[2] + 3 * sin(seq(0, 2 * pi, length.out = 25))),
    NULL)
}
# The line stops at the base of a filled head or at the circle
edge_lines <- imap_dfr(edges, function(e, i) {
  p <- matrix(e[[1]], ncol = 2, byrow = TRUE)
  n <- nrow(p)
  d <- (p[n, ] - p[n - 1, ]) / sqrt(sum((p[n, ] - p[n - 1, ])^2))
  p[n, ] <- p[n, ] - d * switch(e[[2]], block = 6, oval = 6, 0)
  tibble(edge = i, x = p[, 1], y = p[, 2])
})
heads <- imap_dfr(edges, ~ if (.x[[2]] == "none") NULL else mutate(head_geom(.x[[1]], .x[[2]]), edge = .y, end = .x[[2]]))

poly_df <- function(d, fn) pmap_dfr(d, function(...) {
  r <- list(...)
  mutate(fn(r), id = paste(r$x0, r$y0, r$x1, r$y1))
})
region_outer <- poly_df(regions, function(r) rounded_rect(r$x0, r$y0, r$x1, r$y1, 16))
region_inner <- poly_df(regions, function(r) rounded_rect(r$x0 + 2, r$y0 + 2, r$x1 - 2, r$y1 - 2, 14.5))
box_poly <- poly_df(boxes, function(r) rounded_rect(r$x0, r$y0, r$x1, r$y1, if (r$type == "enzyme") 5 else 0)) %>%
  left_join(transmute(boxes, id = paste(x0, y0, x1, y1), type), by = "id")
txt <- function(size) size * u     # font size in canvas units -> geom_text size (mm)

p <- ggplot() +
  geom_polygon(data = region_outer, aes(x, y, group = id), fill = col[["region_fill"]], colour = col[["region"]],
               linewidth = lw) +
  geom_polygon(data = region_inner, aes(x, y, group = id), fill = NA, colour = col[["region"]], linewidth = lw) +
  geom_text(data = regions, aes(lx, ly, label = name, hjust = hjust), colour = col[["region"]], fontface = "bold",
            size = txt(12), family = "sans") +
  geom_rect(aes(xmin = 12, xmax = 122, ymin = 32, ymax = 112), fill = NA, colour = "black", linewidth = lw,
            linetype = "22") +
  geom_text(data = legend, aes(54, y, label = label), hjust = 0, size = txt(10), family = "sans") +
  geom_text(aes(12, 16, label = "Vitamin D metabolism"), hjust = 0, fontface = "bold", size = txt(13), family = "sans") +
  geom_polygon(data = fit(cloud, 126, 31, 246, 93), aes(x, y), fill = col[["sky"]], colour = "black", linewidth = lw) +
  geom_text(aes(186, 63, label = "Vitamin D\nfrom diet"), size = txt(10), lineheight = 0.95, family = "sans") +
  geom_text(aes(179, 130, label = "Absorbed with\ndietary lipids"), hjust = 1, size = txt(9), lineheight = 0.95,
            family = "sans") +
  geom_polygon(data = fit(lightning, 335, 16, 362, 76), aes(x, y), fill = col[["sky"]], colour = "black",
               linewidth = lw) +
  geom_text(aes(408, 40, label = "UV B radiation\nfrom sunlight"), size = txt(10), lineheight = 0.95, family = "sans") +
  geom_text(aes(392, 143, label = "Thermal\nisomerization"), hjust = 0, size = txt(9), lineheight = 0.95,
            family = "sans") +
  geom_path(data = edge_lines, aes(x, y, group = edge), linewidth = lw, linejoin = "mitre") +
  geom_polygon(data = filter(heads, end == "block"), aes(x, y, group = edge), fill = "black", colour = "black",
               linewidth = lw * 0.5) +
  geom_path(data = filter(heads, end == "open"), aes(x, y, group = edge), linewidth = lw) +
  geom_polygon(data = filter(heads, end == "oval"), aes(x, y, group = edge), fill = "white", colour = "black",
               linewidth = lw) +
  geom_polygon(data = box_poly, aes(x, y, group = id, fill = type), colour = "black", linewidth = lw) +
  geom_text(data = boxes, aes((x0 + x1) / 2, (y0 + y1) / 2, label = label,
                              fontface = if_else(type == "enzyme", "bold", "plain")),
            size = txt(10), lineheight = 0.95, family = "sans") +
  scale_fill_manual(values = col, guide = "none") +
  scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
  scale_y_reverse(limits = c(H, 0), expand = c(0, 0)) +
  coord_fixed() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA))

save_fig(p, file.path(proj_dir, "d3pathway", "pathway"), 170, 100)
