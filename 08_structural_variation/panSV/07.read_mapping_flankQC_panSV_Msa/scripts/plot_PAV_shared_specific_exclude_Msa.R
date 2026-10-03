#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_PAV_shared_specific_exclude_Msa.R <per_species.tsv> <out_prefix>")
}
infile <- args[1]
out_prefix <- args[2]

df <- read.table(infile, header=TRUE, sep="\t", quote="", check.names=FALSE)
df <- df[order(-df$specific_count, -df$shared_count, df$species), ]

draw_mirror <- function(file, mode=c("countK", "percent"), device=c("pdf","png","svg")) {
  mode <- match.arg(mode)
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width=4.4, height=5.4, family="Times")
  if (device == "png") png(file, width=1320, height=1620, res=300, type="cairo")
  if (device == "svg") svg(file, width=4.4, height=5.4, family="Times")
  op <- par(no.readonly=TRUE)
  on.exit({par(op); dev.off()}, add=TRUE)

  if (mode == "countK") {
    left <- df$shared_K
    right <- df$specific_K
    left_title <- "Shared (K)"
    right_title <- "Specific (K)"
    label_fmt <- function(x) sprintf("%.1f", x)
    left_ticks <- pretty(c(0, max(left) * 1.08), n=4)
    right_ticks <- pretty(c(0, max(right) * 1.08), n=4)
  } else {
    left <- df$shared_percent
    right <- df$specific_percent
    left_title <- "Shared (%)"
    right_title <- "Specific (%)"
    label_fmt <- function(x) sprintf("%.1f", x)
    left_ticks <- seq(0, 100, by=25)
    right_ticks <- seq(0, 100, by=25)
  }

  n <- nrow(df)
  y <- seq(n, 1)
  bar_h <- 0.72
  left_col <- "#d97932"
  right_col <- "#cbd6ea"

  layout(matrix(c(1,2,3), nrow=1), widths=c(1.45, 0.82, 1.45))

  par(mar=c(1.2, 3.0, 2.5, 0.1), family="Times", xaxs="i", yaxs="i")
  plot(NA, xlim=c(max(left_ticks), 0), ylim=c(0.4, n + 0.6), axes=FALSE, xlab="", ylab="")
  rect(0, y - bar_h/2, left, y + bar_h/2, col=left_col, border=NA)
  axis(3, at=left_ticks, labels=left_ticks, cex.axis=0.82)
  mtext(left_title, side=3, line=1.2, cex=0.95)
  text(left * 0.5, y, label_fmt(left), cex=0.72)
  box(bty="n")

  par(mar=c(1.2, 0.0, 2.5, 0.0), family="Times", xaxs="i", yaxs="i")
  plot(NA, xlim=c(0, 1), ylim=c(0.4, n + 0.6), axes=FALSE, xlab="", ylab="")
  text(0.5, y, df$species, cex=0.82)

  par(mar=c(1.2, 0.1, 2.5, 3.0), family="Times", xaxs="i", yaxs="i")
  plot(NA, xlim=c(0, max(right_ticks)), ylim=c(0.4, n + 0.6), axes=FALSE, xlab="", ylab="")
  rect(0, y - bar_h/2, right, y + bar_h/2, col=right_col, border=NA)
  axis(3, at=right_ticks, labels=right_ticks, cex.axis=0.82)
  mtext(right_title, side=3, line=1.2, cex=0.95)
  text(pmax(right * 0.08, max(right_ticks) * 0.025), y, label_fmt(right), cex=0.72, adj=0)
  box(bty="n")
}

draw_mirror(paste0(out_prefix, ".countK.pdf"), "countK", "pdf")
draw_mirror(paste0(out_prefix, ".countK.png"), "countK", "png")
draw_mirror(paste0(out_prefix, ".countK.svg"), "countK", "svg")
draw_mirror(paste0(out_prefix, ".percent.pdf"), "percent", "pdf")
draw_mirror(paste0(out_prefix, ".percent.png"), "percent", "png")
draw_mirror(paste0(out_prefix, ".percent.svg"), "percent", "svg")
