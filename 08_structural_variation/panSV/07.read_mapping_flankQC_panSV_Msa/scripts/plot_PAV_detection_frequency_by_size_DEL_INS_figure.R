#!/usr/bin/env Rscript
# Extended Data Fig. 6b (revised) figure, drawn from the count table written by
# plot_PAV_detection_frequency_by_size_DEL_INS.R (DEL/INS only; 50-499, 500-4,999, >=5,000 bp).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_PAV_detection_frequency_by_size_DEL_INS_figure.R <count_table.tsv> <out_prefix>")
}
count_tab <- read.delim(args[1], sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
out_prefix <- args[2]

size_levels <- c("50-499", "500-4999", ">=5000")
size_labels <- c("50–499 bp", "500–4,999 bp", "≥5,000 bp")
freq_levels <- sort(unique(count_tab$detection_frequency))
mat <- 100 * xtabs(percent_of_size_bin ~ size_bin + detection_frequency, data = count_tab)
mat <- mat[size_levels, as.character(freq_levels), drop = FALSE]
low_freq_pct <- 100 * sum(count_tab$event_count[count_tab$detection_frequency <= 2]) / sum(count_tab$event_count)
ann_text <- sprintf("%.1f%% in one or two genomes", low_freq_pct)
cols <- c("#5075AF", "#CE954C", "#6BA862")
x_labels <- ifelse(freq_levels <= 6 | freq_levels %% 2 == 0, freq_levels, NA)

plot_one <- function(device, file, width, height, res = 600) {
  if (device == "pdf") {
    cairo_pdf(file, width = width, height = height, family = "Times New Roman")
  } else if (device == "png") {
    png(file, width = width, height = height, units = "in", res = res, type = "cairo", family = "Times New Roman")
  } else if (device == "svg") {
    svg(file, width = width, height = height, family = "Times New Roman")
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(3.2, 4.0, 0.8, 0.6), mgp = c(2.0, 0.45, 0), tcl = -0.3, xaxs = "i", yaxs = "i")
  ylim <- c(0, 90)
  space <- c(0, 0.7)
  bp <- barplot(mat, beside = TRUE, space = space, plot = FALSE)
  xlim <- c(min(bp) - 1.2, max(bp) + 1.2)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
  shade_right <- (max(bp[, 2]) + min(bp[, 3])) / 2
  rect(xlim[1], 0, shade_right, ylim[2], col = "#ECECEC", border = NA)
  barplot(mat, beside = TRUE, space = space, col = cols, border = NA, add = TRUE,
          axes = FALSE, names.arg = rep("", ncol(mat)), xlim = xlim, ylim = ylim)
  axis(2, at = seq(0, 90, 30), labels = paste0(seq(0, 90, 30), "%"), las = 1, lwd = 0.8, cex.axis = 0.9)
  centers <- colMeans(bp)
  show <- !is.na(x_labels)
  axis(1, at = centers[show], labels = x_labels[show], lwd = 0, lwd.ticks = 0.8, cex.axis = 0.9)
  segments(xlim[1], 0, xlim[2], 0, lwd = 0.8)
  segments(xlim[1], 0, xlim[1], ylim[2], lwd = 0.8)
  title(xlab = "Detection frequency", ylab = "PAVs", cex.lab = 1.05)
  text(shade_right + 0.8, 84, ann_text, adj = c(0, 0.5), cex = 0.85)
  legend(x = centers[9], y = 72, legend = size_labels, fill = cols, border = NA,
         bty = "n", cex = 0.72, y.intersp = 1.25, x.intersp = 0.6)
}
plot_one("pdf", paste0(out_prefix, ".pdf"), 4.3, 3.8)
plot_one("png", paste0(out_prefix, ".png"), 4.3, 3.8)
plot_one("svg", paste0(out_prefix, ".svg"), 4.3, 3.8)
cat("DONE", out_prefix, "\n")
