#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) != 1 || hit == length(args)) stop(paste("Missing", flag))
  args[[hit + 1]]
}

overview_path <- get_arg("--overview")
stats_path <- get_arg("--pansv-stats")
out_prefix <- get_arg("--out-prefix")

svtypes <- c("INV", "DUP", "TRA", "INS", "DEL")
legend_order <- c("DEL", "INS", "TRA", "DUP", "INV")
labels <- c(DEL="Deletion", INS="Insertion", TRA="Translocation", DUP="Duplication", INV="Inversion")
cols <- c(DEL="#5AB6C4", INS="#42A98F", TRA="#596782", DUP="#D98955", INV="#D85E5C")
small_svs <- c("INV", "DUP", "TRA")

overview <- read.delim(overview_path, check.names = FALSE, stringsAsFactors = FALSE)
for (sv in names(labels)) {
  if (!sv %in% colnames(overview)) overview[[sv]] <- 0L
}
overview$total_for_plot <- rowSums(overview[, names(labels), drop = FALSE])
overview <- overview[order(overview$total_for_plot, overview$species), ]

stats <- read.delim(stats_path, check.names = FALSE, stringsAsFactors = FALSE)
sv_stats <- stats[stats$section == "SVTYPE", c("key", "count")]
pie_counts <- setNames(rep(0L, length(names(labels))), names(labels))
for (i in seq_len(nrow(sv_stats))) {
  key <- sv_stats$key[i]
  if (key %in% names(pie_counts)) pie_counts[[key]] <- as.integer(sv_stats$count[i])
}

y_fmt <- function(x) {
  ifelse(x == 0, "0", ifelse(x %% 1000 == 0, paste0(as.integer(x / 1000), "k"), paste0(format(round(x / 1000, 1), nsmall = 1), "k")))
}

bar_segments <- function(x_centers, mat, width = 0.64) {
  bottom <- rep(0, ncol(mat))
  for (sv in rownames(mat)) {
    vals <- as.numeric(mat[sv, ])
    rect(x_centers - width/2, bottom, x_centers + width/2, bottom + vals,
         col = cols[sv], border = NA)
    bottom <- bottom + vals
  }
}

draw_data_pie <- function(cx, cy, ry, values, colors, title = NULL) {
  usr <- par("usr")
  pin <- par("pin")
  rx <- ry * diff(usr[1:2]) / diff(usr[3:4]) * pin[2] / pin[1]
  vals <- as.numeric(values)
  fr <- vals / sum(vals)
  start <- pi / 2
  for (i in seq_along(vals)) {
    end <- start - 2 * pi * fr[i]
    theta <- seq(start, end, length.out = max(20, ceiling(120 * fr[i])))
    polygon(c(cx, cx + rx * cos(theta)), c(cy, cy + ry * sin(theta)),
            col = colors[i], border = "white", lwd = 0.7)
    start <- end
  }
  if (!is.null(title)) text(cx, cy + ry * 1.28, title, cex = 0.72, font = 2)
  invisible(rx)
}

add_pie_labels <- function(cx, cy, ry, values, pie_svs) {
  usr <- par("usr")
  pin <- par("pin")
  rx <- ry * diff(usr[1:2]) / diff(usr[3:4]) * pin[2] / pin[1]
  fr <- values / sum(values)
  starts <- pi / 2 - c(0, cumsum(fr[-length(fr)]) * 2 * pi)
  ends <- pi / 2 - cumsum(fr) * 2 * pi
  mids <- (starts + ends) / 2
  fmt <- function(v) format(v, big.mark = ",", scientific = FALSE)
  for (sv in c("DEL", "INS")) {
    i <- match(sv, pie_svs)
    a <- mids[i]
    x1 <- cx + rx * 0.92 * cos(a); y1 <- cy + ry * 0.92 * sin(a)
    x2 <- cx + rx * 1.18 * cos(a); y2 <- cy + ry * 1.18 * sin(a)
    segments(x1, y1, x2, y2, col = "grey30", lwd = 0.6)
    text(x2 + ifelse(cos(a) > 0, 0.25, -0.25), y2, fmt(values[i]),
         adj = ifelse(cos(a) > 0, 0, 1), cex = 0.62)
  }
  small <- c("TRA", "DUP", "INV")
  small_y <- cy + ry * c(1.16, 1.02, 0.88)
  small_x <- cx + rx * 1.50
  for (j in seq_along(small)) {
    sv <- small[j]
    i <- match(sv, pie_svs)
    a <- mids[i]
    x1 <- cx + rx * 0.95 * cos(a); y1 <- cy + ry * 0.95 * sin(a)
    segments(x1, y1, small_x - 0.08, small_y[j], col = "grey30", lwd = 0.55)
    text(small_x, small_y[j], paste0(sv, " ", fmt(values[i])), adj = 0, cex = 0.55)
  }
}

plot_one <- function() {
  mat <- as.matrix(t(overview[, svtypes, drop = FALSE]))
  colnames(mat) <- overview$species
  totals <- colSums(mat)
  x <- seq_len(ncol(mat))
  xlim <- c(0.35, ncol(mat) + 0.65)

  lower_step <- 500       # 0.5k
  upper_step <- 40000     # 40k
  small_max <- max(colSums(mat[small_svs, , drop = FALSE]))
  lower_max <- ceiling((small_max * 1.18) / lower_step) * lower_step
  lower_max <- max(lower_max, 3200)
  lower_ticks <- seq(0, lower_max, by = lower_step)

  upper_min <- floor(min(totals) * 0.55 / upper_step) * upper_step
  upper_min <- max(upper_min, 40000)
  upper_max <- ceiling(max(totals) * 1.10 / upper_step) * upper_step
  upper_ticks <- seq(upper_min, upper_max, by = upper_step)

  layout(matrix(c(1, 2), nrow = 2), heights = c(2.45, 1.45))

  par(family = "sans", mar = c(0.15, 4.95, 0.7, 0.9), xaxs = "i", yaxs = "i", cex = 0.9)
  plot(NA, xlim = xlim, ylim = c(upper_min, upper_max), axes = FALSE, xlab = "", ylab = "")
  bar_segments(x, mat)
  axis(2, at = upper_ticks, labels = y_fmt(upper_ticks), las = 1, tck = -0.015)
  abline(h = upper_ticks, col = adjustcolor("grey70", 0.18), lwd = 0.6)
  box(bty = "l", lwd = 1.05)
  mtext("SV number", side = 2, line = 3.45, cex = 0.98)
  legend("topleft", legend = labels[legend_order], fill = cols[legend_order], border = NA,
         bty = "n", cex = 0.78, y.intersp = 0.80, x.intersp = 0.48, inset = c(0.01, 0.02))
  text(par("usr")[1] - 0.10 * diff(par("usr")[1:2]), upper_max * 0.995, "A", font = 2, cex = 1.35)

  pie_svs <- c("DEL", "INS", "TRA", "DUP", "INV")
  pie_vals <- as.numeric(pie_counts[pie_svs])
  pie_cx <- mean(xlim) + 0.8
  pie_cy <- upper_min + 0.63 * (upper_max - upper_min)
  pie_ry <- 0.19 * (upper_max - upper_min)
  draw_data_pie(pie_cx, pie_cy, pie_ry, pie_vals, cols[pie_svs], title = "Non-redundant panSV")
  add_pie_labels(pie_cx, pie_cy, pie_ry, pie_vals, pie_svs)

  par(family = "sans", mar = c(5.55, 4.95, 0.05, 0.9), xaxs = "i", yaxs = "i", cex = 0.9)
  plot(NA, xlim = xlim, ylim = c(0, lower_max), axes = FALSE, xlab = "", ylab = "")
  bar_segments(x, mat)
  axis(2, at = lower_ticks, labels = y_fmt(lower_ticks), las = 1, tck = -0.015, cex.axis = 0.63)
  axis(1, at = x, labels = FALSE, tck = -0.015)
  abline(h = lower_ticks, col = adjustcolor("grey70", 0.16), lwd = 0.45)
  box(bty = "l", lwd = 1.05)
  text(x = x, y = par("usr")[3] - 0.065 * lower_max,
       labels = overview$species, srt = 45, adj = 1, cex = 0.68)
  mtext("Species", side = 1, line = 4.45, cex = 0.95)

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  segments(0.088, 0.366, 0.106, 0.389, xpd = NA, lwd = 1.0)
  segments(0.088, 0.384, 0.106, 0.407, xpd = NA, lwd = 1.0)
  text(0.09, 0.018,
       "Bars: per-species merged SVs after read-mapping flank QC; rare SV classes are stacked at the base and zoomed below. Pie: non-redundant panSV events.",
       adj = c(0, 0), cex = 0.55, col = "grey35")
}

out_dir <- dirname(out_prefix)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- paste0(out_prefix, ".R.split.v4.pdf")
svg_path <- paste0(out_prefix, ".R.split.v4.svg")
png_path <- paste0(out_prefix, ".R.split.v4.png")
tsv_path <- paste0(out_prefix, ".R.split.v4.sorted_species_counts.tsv")

pdf(pdf_path, width = 7.4, height = 4.9, useDingbats = FALSE)
plot_one()
dev.off()

svg(svg_path, width = 7.4, height = 4.9, pointsize = 10)
plot_one()
dev.off()

png(png_path, width = 7.4, height = 4.9, units = "in", res = 500, type = "cairo-png")
plot_one()
dev.off()

sorted_out <- overview[, c("species", "total_for_plot", names(labels))]
colnames(sorted_out)[colnames(sorted_out) == "total_for_plot"] <- "total"
write.table(sorted_out, file = tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("WROTE\t", pdf_path, "\n", sep = "")
cat("WROTE\t", svg_path, "\n", sep = "")
cat("WROTE\t", png_path, "\n", sep = "")
cat("WROTE\t", tsv_path, "\n", sep = "")
cat("Y_TICKS\tlower=0.5k\tupper=40k\n")
cat("STACK_ORDER_BOTTOM_TO_TOP\t", paste(svtypes, collapse = ","), "\n", sep = "")
cat("PIE_COUNTS\t", paste(paste0(names(labels), "=", pie_counts[names(labels)]), collapse = "\t"), "\n", sep = "")