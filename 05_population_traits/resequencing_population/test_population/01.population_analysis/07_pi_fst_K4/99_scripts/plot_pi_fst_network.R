args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript plot_pi_fst_network_K4_no_outgroup.R <analysis_dir>")
}

base_dir <- normalizePath(args[1], mustWork = TRUE)
pi_file <- file.path(base_dir, "pi_summary.tsv")
fst_file <- file.path(base_dir, "fst_summary.tsv")
plot_dir <- file.path(base_dir, "04_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pi_df <- read.delim(pi_file, check.names = FALSE)
fst_df <- read.delim(fst_file, check.names = FALSE)

groups <- c("C1", "C2", "C3", "C4")
pi_df <- pi_df[match(groups, pi_df$group), ]

coords <- data.frame(
  group = groups,
  x = c(-1, 0, 1, 0),
  y = c(0, 1, 0, -1),
  color = c("#ef6a6a", "#405a93", "#159c89", "#b786bd"),
  stringsAsFactors = FALSE
)
coords$pi <- pi_df$pi_snp_mean_valid_sites

pi_range <- range(coords$pi, finite = TRUE)
if (diff(pi_range) == 0) {
  coords$radius <- 0.31
} else {
  coords$radius <- 0.25 + 0.08 * (coords$pi - pi_range[1]) / diff(pi_range)
}

edge_key <- function(a, b) paste(sort(c(a, b)), collapse = "-")
fst_df$key <- mapply(edge_key, fst_df$group1, fst_df$group2)
fst_lookup <- setNames(fst_df$weighted_fst, fst_df$key)

edges <- data.frame(
  group1 = c("C1", "C2", "C1", "C1", "C2", "C3"),
  group2 = c("C2", "C3", "C3", "C4", "C4", "C4"),
  label_x = c(-0.72, 0.72, -0.18, -0.72, 0.40, 0.72),
  label_y = c(0.68, 0.68, 0.23, -0.68, -0.02, -0.68),
  stringsAsFactors = FALSE
)
edges$key <- mapply(edge_key, edges$group1, edges$group2)
edges$fst <- fst_lookup[edges$key]

draw_label <- function(x, y, label, cex = 1.35) {
  pad_x <- 0.045
  pad_y <- 0.035
  w <- strwidth(label, cex = cex)
  h <- strheight(label, cex = cex)
  rect(x - w / 2 - pad_x, y - h / 2 - pad_y,
       x + w / 2 + pad_x, y + h / 2 + pad_y,
       col = "white", border = NA)
  text(x, y, label, cex = cex, col = "#111111")
}

draw_panel <- function() {
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = c(-1.45, 1.45), ylim = c(-1.35, 1.35), asp = 1,
       axes = FALSE, xlab = "", ylab = "", bty = "n")

  for (i in seq_len(nrow(edges))) {
    p1 <- coords[coords$group == edges$group1[i], ]
    p2 <- coords[coords$group == edges$group2[i], ]
    segments(p1$x, p1$y, p2$x, p2$y, lty = 2, lwd = 1.4, col = "#222222")
  }

  for (i in seq_len(nrow(edges))) {
    draw_label(edges$label_x[i], edges$label_y[i], sprintf("%.3f", edges$fst[i]))
  }

  for (i in seq_len(nrow(coords))) {
    symbols(coords$x[i], coords$y[i], circles = coords$radius[i], inches = FALSE,
            add = TRUE, bg = coords$color[i], fg = coords$color[i], lwd = 1.5)
    text(coords$x[i], coords$y[i] + 0.055, coords$group[i],
         col = "white", cex = 1.25, font = 2)
    text(coords$x[i], coords$y[i] - 0.085, sprintf("%.3f", coords$pi[i]),
         col = "white", cex = 1.12, font = 2)
  }
}

png(file.path(plot_dir, "pi_fst_K4_network.png"), width = 1200, height = 1000, res = 220, type = "cairo")
draw_panel()
dev.off()

pdf(file.path(plot_dir, "pi_fst_K4_network.pdf"), width = 5.5, height = 4.8, useDingbats = FALSE)
draw_panel()
dev.off()

writeLines(c(
  file.path(plot_dir, "pi_fst_K4_network.png"),
  file.path(plot_dir, "pi_fst_K4_network.pdf")
))
