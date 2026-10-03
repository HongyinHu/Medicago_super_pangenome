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

svtypes <- c("DEL", "INS", "TRA", "DUP", "INV")
labels <- c(DEL="Deletion", INS="Insertion", TRA="Translocation", DUP="Duplication", INV="Inversion")
cols <- c(DEL="#56B4C3", INS="#45A68C", TRA="#5A6682", DUP="#D88A56", INV="#D9635F")

overview <- read.delim(overview_path, check.names = FALSE, stringsAsFactors = FALSE)
for (sv in svtypes) {
  if (!sv %in% colnames(overview)) overview[[sv]] <- 0L
}
overview$total_for_plot <- rowSums(overview[, svtypes, drop = FALSE])
overview <- overview[order(overview$total_for_plot, overview$species), ]

stats <- read.delim(stats_path, check.names = FALSE, stringsAsFactors = FALSE)
sv_stats <- stats[stats$section == "SVTYPE", c("key", "count")]
pie_counts <- setNames(rep(0L, length(svtypes)), svtypes)
for (i in seq_len(nrow(sv_stats))) {
  key <- sv_stats$key[i]
  if (key %in% names(pie_counts)) pie_counts[[key]] <- as.integer(sv_stats$count[i])
}

y_fmt <- function(x) {
  ifelse(x >= 1000, paste0(round(x / 1000), "k"), as.character(x))
}

plot_one <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(family = "sans", mar = c(5.4, 4.7, 1.1, 1.0), xpd = NA, cex = 0.9)
  mat <- as.matrix(t(overview[, svtypes, drop = FALSE]))
  ymax <- max(colSums(mat)) * 1.18
  bp <- barplot(mat,
                col = cols[svtypes], border = NA, space = 0.55,
                ylim = c(0, ymax), axes = FALSE, names.arg = rep("", ncol(mat)),
                ylab = "SV number", xlab = "")
  axis(2, at = pretty(c(0, ymax)), labels = y_fmt(pretty(c(0, ymax))), las = 1, tck = -0.015)
  axis(1, at = bp, labels = FALSE, tck = -0.015)
  abline(h = pretty(c(0, ymax)), col = adjustcolor("grey70", 0.18), lwd = 0.6)
  box(bty = "l", lwd = 1.1)
  mtext("Species", side = 1, line = 4.2, cex = 0.95)
  text(x = bp, y = par("usr")[3] - 0.035 * ymax,
       labels = overview$species, srt = 45, adj = 1, cex = 0.72)
  legend("topleft", legend = labels[svtypes], fill = cols[svtypes], border = NA,
         bty = "n", cex = 0.82, y.intersp = 0.85, x.intersp = 0.55)
  text(par("usr")[1] - 0.12 * diff(par("usr")[1:2]), ymax * 1.02, "A", font = 2, cex = 1.45)

  # Inset pie chart: non-redundant panSV SVTYPE counts.
  par(fig = c(0.37, 0.70, 0.55, 0.93), new = TRUE, mar = c(0, 0, 1.0, 0))
  pie_vals <- as.numeric(pie_counts[svtypes])
  pie_labels <- ifelse(pie_vals > 0, format(pie_vals, big.mark = ",", scientific = FALSE), "")
  pie(pie_vals,
      labels = pie_labels,
      col = cols[svtypes], border = "white", lwd = 0.7,
      clockwise = TRUE, init.angle = 90,
      cex = 0.72)
  title("Non-redundant panSV", cex.main = 0.78, line = -0.2)

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  text(0.09, 0.035,
       "Bars: per-species merged SVs after read-mapping flank QC; Pie: non-redundant panSV events.",
       adj = c(0, 0), cex = 0.58, col = "grey35")
}

out_dir <- dirname(out_prefix)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- paste0(out_prefix, ".R.pdf")
svg_path <- paste0(out_prefix, ".R.svg")
png_path <- paste0(out_prefix, ".R.png")
tsv_path <- paste0(out_prefix, ".R.sorted_species_counts.tsv")

pdf(pdf_path, width = 7.2, height = 4.6, useDingbats = FALSE)
plot_one()
dev.off()

svg(svg_path, width = 7.2, height = 4.6, pointsize = 10)
plot_one()
dev.off()

png(png_path, width = 7.2, height = 4.6, units = "in", res = 450, type = "cairo-png")
plot_one()
dev.off()

sorted_out <- overview[, c("species", "total_for_plot", svtypes)]
colnames(sorted_out)[colnames(sorted_out) == "total_for_plot"] <- "total"
write.table(sorted_out, file = tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("WROTE\t", pdf_path, "\n", sep = "")
cat("WROTE\t", svg_path, "\n", sep = "")
cat("WROTE\t", png_path, "\n", sep = "")
cat("WROTE\t", tsv_path, "\n", sep = "")
cat("PIE_COUNTS\t", paste(paste0(svtypes, "=", pie_counts[svtypes]), collapse = "\t"), "\n", sep = "")