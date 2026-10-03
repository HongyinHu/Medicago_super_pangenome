#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) != 1 || hit == length(args)) stop(paste("Missing", flag))
  args[[hit + 1]]
}
in_tsv <- get_arg("--input")
out_prefix <- get_arg("--out-prefix")

df <- read.delim(in_tsv, check.names = FALSE, stringsAsFactors = FALSE)
sv_order <- c("DEL", "INS")
bin_order <- c(as.character(1:10), ">10")
cols <- c(DEL = "#5AB6C4", INS = "#42A98F")
legend_labels <- c(DEL = "Deletion", INS = "Insertion")

df$SVTYPE <- factor(df$SVTYPE, levels = sv_order)
df$bin_label <- factor(df$bin_label, levels = bin_order)

fmt_y <- function(x) {
  out <- ifelse(x >= 1000, paste0(round(x / 1000), "k"), as.character(x))
  out[x == 0] <- "0"
  out
}

plot_one <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)

  y_max <- 300000
  y_ticks <- seq(0, y_max, by = 100000)

  # Equal layout heights + identical per-panel margins make the physical y-scale identical.
  layout(matrix(c(1, 2), ncol = 1), heights = c(1, 1))
  par(oma = c(3.0, 0.0, 0.0, 0.0))

  for (sv in sv_order) {
    sub <- df[df$SVTYPE == sv, ]
    vals <- sub$event_count[match(bin_order, as.character(sub$bin_label))]
    bottom_panel <- sv == "INS"

    par(
      family = "serif",
      mar = c(1.15, 4.65, 0.45, 0.80),
      mgp = c(1.75, 0.45, 0),
      xaxs = "i",
      yaxs = "i",
      tcl = -0.22,
      cex.axis = 0.9
    )

    bp <- barplot(
      vals,
      col = cols[sv],
      border = NA,
      ylim = c(0, y_max),
      axes = FALSE,
      names.arg = rep("", length(vals)),
      space = 0.25
    )

    axis(2, at = y_ticks, labels = fmt_y(y_ticks), las = 1, lwd = 1.05, lwd.ticks = 1.05)
    box(bty = "l", lwd = 1.05)

    legend(
      "topright",
      legend = legend_labels[sv],
      fill = cols[sv],
      border = NA,
      bty = "n",
      cex = 1.02,
      inset = c(0.09, 0.05)
    )

    if (bottom_panel) {
      axis(1, at = bp, labels = bin_order, lwd = 1.05, lwd.ticks = 1.05, cex.axis = 0.9)
    }
  }

  mtext("PAV Length (kb)", side = 1, outer = TRUE, line = 1.7, cex = 1.05, family = "serif")

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  text(0.030, 0.535, "PAV number", srt = 90, cex = 1.05, family = "serif")
}

out_dir <- dirname(out_prefix)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (ext in c("pdf", "svg", "png")) {
  out <- paste0(out_prefix, ".", ext)
  if (file.exists(out)) unlink(out)
  if (ext == "pdf") pdf(out, width = 4.2, height = 5.4, useDingbats = FALSE)
  if (ext == "svg") svg(out, width = 4.2, height = 5.4, pointsize = 10)
  if (ext == "png") png(out, width = 4.2, height = 5.4, units = "in", res = 500, type = "cairo-png")
  plot_one()
  dev.off()
  cat("WROTE\t", out, "\n", sep = "")
}
