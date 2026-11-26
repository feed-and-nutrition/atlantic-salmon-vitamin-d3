# Outline of an Atlantic salmon for the graphical abstract (09_2), traced from a public-domain drawing.
# Source: Timothy Knepp, U.S. Fish and Wildlife Service, "Salmo salar" (Wikimedia Commons, File:Salmo_salar.jpg,
# public domain as a work of the US federal government), reduced to 1236 x 404 px. Only the outline is used.
# The silhouette is the largest contour of the non-white pixels. Each fin that sticks out of it (dorsal, adipose,
# caudal, anal, pelvic) is cut off at its base, where the body outline is bridged by a smooth curve. The open mouth
# is kept as traced. The pectoral fin, which lies on the body, the eye, the gill cover and jaw lines, the lateral
# line and the organs (brain, kidney, liver) were read off the drawing by hand (pixel coordinates below). The dark back
# (above the lateral line) and the black spots follow the colouring of Atlantic salmon, with spots placed at random.
# Coordinates are written in body lengths: snout at (0, 0), tip of the caudal fin at x = 1, y up.
# Not in run_all (the outline only changes if this script does).
# Input:  d3pathway/salmo_salar_knepp.png
# Output: d3pathway/salmon_outline.tsv (part, group, x, y)

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

img <- png::readPNG(file.path(proj_dir, "d3pathway", "salmo_salar_knepp.png"))
fish <- (pmin(img[, , 1], img[, , 2], img[, , 3]) < 0.92) * 1
iso <- isoband::isolines(seq_len(ncol(fish)), seq_len(nrow(fish)), fish, 0.5)[[1]]  # y = pixel row, down
rings <- split(tibble(x = iso$x, y = iso$y), iso$id)
ring <- rings[[which.max(map_int(rings, nrow))]]
# Smooth the pixel steps (circular moving average over 11 points) and keep every fourth point
k <- 5
smooth_ring <- function(v) {
  n <- length(v)
  as.numeric(stats::filter(c(v[(n - k + 1):n], v, v[1:k]), rep(1 / (2 * k + 1), 2 * k + 1)))[(k + 1):(n + k)]
}
ring <- tibble(x = smooth_ring(ring$x), y = smooth_ring(ring$y))[seq(1, nrow(ring), by = 4), ]
n <- nrow(ring)
near <- function(x, y) which.min((ring$x - x)^2 + (ring$y - y)^2)  # ring point closest to a pixel position
from_to <- function(a, b) if (a <= b) a:b else c(a:n, 1:b)        # ring indices from a to b, wrapping round

# Smooth bridge (cubic Hermite) from ring point a to ring point b, following the outline's direction at both ends
bridge <- function(a, b, m = 24) {
  p0 <- unlist(ring[a, ])
  p1 <- unlist(ring[b, ])
  d <- sqrt(sum((p1 - p0)^2))
  tangent <- function(i, j) { v <- unlist(ring[j, ]) - unlist(ring[i, ]); v / sqrt(sum(v^2)) * d }
  t0 <- tangent((a - 4 - 1) %% n + 1, a)
  t1 <- tangent(b, b %% n + 4)
  s <- seq(0, 1, length.out = m)
  h <- cbind(2 * s^3 - 3 * s^2 + 1, s^3 - 2 * s^2 + s, -2 * s^3 + 3 * s^2, s^3 - s^2)
  tibble(x = h %*% c(p0[1], t0[1], p1[1], t1[1]) %>% as.numeric(),
         y = h %*% c(p0[2], t0[2], p1[2], t1[2]) %>% as.numeric())
}

# Fin bases on the outline (pixel positions), and the two ends of the margin that the rays run to
fins <- tribble(~fin, ~a_x, ~a_y, ~b_x, ~b_y, ~m0_x, ~m0_y, ~m1_x, ~m1_y, ~rays,
  "dorsal", 432, 103, 590, 107, 535, 24, 620, 96, 10,
  "adipose", 860, 130, 900, 143, NA, NA, NA, NA, 0,
  "caudal", 1005, 164, 1006, 252, 1189, 76, 1187, 324, 15,
  "anal", 903, 261, 770, 302, 927, 268, 882, 337, 8,
  "pelvic", 597, 334, 525, 336, 600, 338, 616, 371, 7) %>%
  mutate(a = map2_int(a_x, a_y, near), b = map2_int(b_x, b_y, near))

# Body: the outline with every fin bridged, from the lower to the upper end of the caudal base.
# The caudal base is the only edge without an outline.
cuts <- select(fins, a, b)
caudal <- filter(fins, fin == "caudal")
body <- tibble()
i <- caudal$b
repeat {
  cut <- filter(cuts, a == i)
  if (nrow(cut) && cut$a != caudal$a) {
    body <- bind_rows(body, bridge(cut$a, cut$b))
    i <- cut$b %% n + 1
  } else {
    body <- bind_rows(body, ring[i, ])
    if (i == caudal$a) break
    i <- i %% n + 1
  }
}

# Smooth the body outline behind the head once more (moving average over 9 points), which removes the unevenness
# from the shading of the light belly. The head and jaws (x < 130 px) keep the traced outline.
smooth_path <- function(v, w = 4) map_dbl(seq_along(v), \(j) mean(v[max(1, j - w):min(length(v), j + w)]))
body <- mutate(body, x = if_else(x < 130, x, smooth_path(x)), y = if_else(x < 130, y, smooth_path(y)))

# Fins: the part of the outline outside the body, closed along the bridge (hidden behind the body)
fin_poly <- pmap_dfr(fins, function(fin, a, b, ...) {
  # The caudal fin reaches 15 px into the body, so no gap opens at its base
  closing <- if (fin == "caudal") mutate(ring[c(b, b, a, a), ], x = x - c(0, 15, 15, 0)) else bridge(a, b)[24:1, ]
  mutate(bind_rows(ring[from_to(a, b), ], closing), group = fin)
})
# Rays from evenly spaced points on the base to evenly spaced points on the margin
along <- function(d, m) {  # m points evenly spaced along a path
  s <- c(0, cumsum(sqrt(diff(d$x)^2 + diff(d$y)^2)))
  tibble(x = approx(s, d$x, seq(0, max(s), length.out = m))$y, y = approx(s, d$y, seq(0, max(s), length.out = m))$y)
}
fin_rays <- pmap_dfr(filter(fins, rays > 0), function(fin, a, b, m0_x, m0_y, m1_x, m1_y, rays, ...) {
  base <- if (fin == "caudal") ring[c(a, b), ] else bridge(a, b)
  base <- along(base, rays + 2)[2:(rays + 1), ]
  margin <- along(ring[from_to(near(m0_x, m0_y), near(m1_x, m1_y)), ], rays)
  tibble(group = paste(fin, seq_len(rays)), x = base$x, y = base$y, xend = margin$x, yend = margin$y) %>%
    pivot_longer(c(x, xend), values_to = "px") %>%
    mutate(py = if_else(name == "x", y, yend)) %>%
    select(group, x = px, y = py)
})

# Read off the drawing (pixels)
pts <- function(...) { v <- c(...); tibble(x = v[c(TRUE, FALSE)], y = v[c(FALSE, TRUE)]) }
open_spline <- function(d, m = 40) {  # smooth curve through the points
  t <- seq_len(nrow(d))
  tibble(x = spline(t, d$x, n = m, method = "natural")$y, y = spline(t, d$y, n = m, method = "natural")$y)
}
closed_spline <- function(d, m = 80) {  # smooth closed curve through the points
  d <- bind_rows(d, d[1, ])
  t <- seq_len(nrow(d))
  tibble(x = spline(t, d$x, n = m, method = "periodic")$y, y = spline(t, d$y, n = m, method = "periodic")$y)
}
ellipse <- function(cx, cy, rx, ry, m = 40) {
  a <- seq(0, 2 * pi, length.out = m)
  tibble(x = cx + rx * cos(a), y = cy + ry * sin(a))
}
pectoral <- closed_spline(pts(252, 283, 300, 272, 345, 266, 372, 266, 363, 280, 335, 295, 295, 303, 264, 298), m = 60)
pectoral_rays <- {
  base <- along(pts(254, 283, 264, 297), 7)
  margin <- along(pts(370, 267, 362, 280, 335, 294, 297, 301), 7)
  tibble(group = paste("pectoral", 1:7), x = base$x, y = base$y, xend = margin$x, yend = margin$y) %>%
    pivot_longer(c(x, xend), values_to = "px") %>%
    mutate(py = if_else(name == "x", y, yend)) %>%
    select(group, x = px, y = py)
}
lateral_line <- along(pts(245, 200, 600, 201, 1005, 207), 30)
# Dark back: the body above the lateral line, closed in the head along a line above the eye to the upper jaw tip
snout <- which.min((body$x - 39)^2 + (body$y - 217)^2)
upper <- body[snout:nrow(body), ]                     # upper jaw tip, back, upper end of the caudal base
back <- bind_rows(upper, lateral_line[nrow(lateral_line):1, ], pts(200, 196, 150, 195, 110, 196, 70, 207, 45, 214))
# Black spots on the back, mostly above the lateral line as in Atlantic salmon, plus two on the gill cover. Spots
# keep clear of the dorsal fin base edge, the brain and the kidney, and of the kidney callout line (x near 734 px).
set.seed(1)
top_at <- function(x) approx(upper$x, upper$y, x, ties = mean)$y
spot_xy <- tibble(x = rep(seq(270, 960, by = 46), each = 2) + runif(32, -12, 12)) %>%
  mutate(y = top_at(x) + runif(n(), 0.25, 0.85) * (183 - 6 - top_at(x)) + 6) %>%
  filter(abs(x - 734) > 14) %>%
  bind_rows(pts(205, 212, 222, 228))
spots <- pmap_dfr(mutate(spot_xy, k = row_number()), \(x, y, k) mutate(ellipse(x, y, 4.2, 3.6, 16), group = paste("spot", k)))

other <- bind_rows(
  outline = body,
  back = back,
  pectoral = pectoral,
  operculum = open_spline(pts(236, 195, 243, 233, 230, 263, 207, 282, 177, 290)),
  preopercle = open_spline(pts(173, 203, 177, 240, 172, 262, 158, 275)),
  maxilla = open_spline(pts(62, 224, 90, 236, 118, 256)),
  lateral_line = lateral_line,
  eye_iris = ellipse(120, 212, 13.5, 13.5),
  eye_pupil = ellipse(120, 212, 8.5, 8.5),
  eye_glint = ellipse(123.5, 208.5, 2.8, 2.8),
  brain = ellipse(160, 197, 20, 11),
  liver = ellipse(320, 250, 60, 18),
  kidney = closed_spline(pts(244, 204, 252, 186, 272, 180, 300, 182, 330, 190, 380, 193, 500, 194, 700, 196,
                             790, 199, 812, 204, 790, 210, 700, 211, 500, 212, 380, 212, 330, 214, 300, 221,
                             272, 223, 252, 219), m = 120),
  .id = "part") %>%
  mutate(group = part)

# To body lengths: snout at (0, 0) on the mouth line, tip of the caudal fin at x = 1, y up
x0 <- min(body$x)
y0 <- 233
len <- max(ring$x) - x0
out <- bind_rows(mutate(body, part = "body", group = "body"),
                 mutate(fin_poly, part = "fin"),
                 mutate(fin_rays, part = "ray"),
                 mutate(pectoral_rays, part = "pectoral_ray"),
                 mutate(spots, part = "spot"),
                 other) %>%
  transmute(part, group, x = round((x - x0) / len, 4), y = round((y0 - y) / len, 4))
write_tsv(out, file.path(proj_dir, "d3pathway", "salmon_outline.tsv"))
message("Wrote d3pathway/salmon_outline.tsv: ", n_distinct(out$group), " shapes, ", nrow(out), " points")
