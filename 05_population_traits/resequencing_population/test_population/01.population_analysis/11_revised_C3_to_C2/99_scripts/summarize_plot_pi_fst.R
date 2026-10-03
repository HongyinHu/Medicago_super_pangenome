args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript summarize_plot_pi_fst.R <out_dir>")
}

out_dir <- normalizePath(args[1], mustWork = TRUE)
groups <- c("C1", "C2", "C3", "C4")
counts <- read.delim(file.path(out_dir, "00_groups", "revised_group_counts.tsv"), stringsAsFactors = FALSE)
counts <- counts[counts$group %in% groups, ]

pi_rows <- lapply(groups, function(g) {
  f <- file.path(out_dir, "03_pi_fst", "02_pi", paste0(g, ".sites.pi"))
  x <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  vals <- suppressWarnings(as.numeric(x$PI))
  data.frame(
    group = g,
    n_samples = counts$n[match(g, counts$group)],
    n_sites = length(vals),
    n_valid_sites = sum(is.finite(vals)),
    n_nan_sites = sum(!is.finite(vals)),
    pi_snp_mean_valid_sites = mean(vals[is.finite(vals)]),
    stringsAsFactors = FALSE
  )
})
pi_summary <- do.call(rbind, pi_rows)
write.table(pi_summary, file.path(out_dir, "03_pi_fst", "pi_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parse_fst_log <- function(log_file) {
  lines <- readLines(log_file, warn = FALSE)
  mean_line <- grep("Weir and Cockerham mean Fst estimate", lines, value = TRUE)
  weighted_line <- grep("Weir and Cockerham weighted Fst estimate", lines, value = TRUE)
  mean_fst <- as.numeric(sub(".*: *", "", mean_line[length(mean_line)]))
  weighted_fst <- as.numeric(sub(".*: *", "", weighted_line[length(weighted_line)]))
  c(mean_fst = mean_fst, weighted_fst = weighted_fst)
}

pairs <- list(c("C1", "C2"), c("C1", "C3"), c("C1", "C4"), c("C2", "C3"), c("C2", "C4"), c("C3", "C4"))
fst_rows <- lapply(pairs, function(p) {
  log_file <- file.path(out_dir, "03_pi_fst", "03_fst", paste0(p[1], "_vs_", p[2], ".stdout.log"))
  vals <- parse_fst_log(log_file)
  data.frame(group1 = p[1], group2 = p[2], mean_fst = vals["mean_fst"], weighted_fst = vals["weighted_fst"], stringsAsFactors = FALSE)
})
fst_summary <- do.call(rbind, fst_rows)
write.table(fst_summary, file.path(out_dir, "03_pi_fst", "fst_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_dir <- file.path(out_dir, "04_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

coords <- data.frame(
  group = groups,
  x = c(-1, 0, 1, 0),
  y = c(0, 1, 0, -1),
  color = c("#ef6a6a", "#405a93", "#159c89", "#b786bd"),
  stringsAsFactors = FALSE
)
coords$pi <- pi_summary$pi_snp_mean_valid_sites[match(coords$group, pi_summary$group)]
pi_range <- range(coords$pi, finite = TRUE)
coords$radius <- if (diff(pi_range) == 0) 0.31 else 0.25 + 0.08 * (coords$pi - pi_range[1]) / diff(pi_range)

edge_key <- function(a, b) paste(sort(c(a, b)), collapse = "-")
fst_summary$key <- mapply(edge_key, fst_summary$group1, fst_summary$group2)
fst_lookup <- setNames(fst_summary$weighted_fst, fst_summary$key)
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
  rect(x - w / 2 - pad_x, y - h / 2 - pad_y, x + w / 2 + pad_x, y + h / 2 + pad_y, col = "white", border = NA)
  text(x, y, label, cex = cex, col = "#111111")
}

draw_panel <- function() {
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = c(-1.45, 1.45), ylim = c(-1.35, 1.35), asp = 1, axes = FALSE, xlab = "", ylab = "", bty = "n")
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
    text(coords$x[i], coords$y[i] + 0.055, coords$group[i], col = "white", cex = 1.25, font = 2)
    text(coords$x[i], coords$y[i] - 0.085, sprintf("%.3f", coords$pi[i]), col = "white", cex = 1.12, font = 2)
  }
}

png(file.path(plot_dir, "pi_fst_revised_C3_to_C2_network.png"), width = 1200, height = 1000, res = 220, type = "cairo")
draw_panel()
dev.off()
pdf(file.path(plot_dir, "pi_fst_revised_C3_to_C2_network.pdf"), width = 5.5, height = 4.8, useDingbats = FALSE)
draw_panel()
dev.off()

print(pi_summary)
print(fst_summary)
