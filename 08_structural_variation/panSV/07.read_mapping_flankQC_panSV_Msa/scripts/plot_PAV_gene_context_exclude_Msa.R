#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_PAV_gene_context_exclude_Msa.R <prefix> <out_prefix>")
}
prefix <- args[1]
out_prefix <- args[2]

dist <- read.table(paste0(prefix, ".distance_bins.tsv"), header=TRUE, sep="\t", quote="", check.names=FALSE)
catdf <- read.table(paste0(prefix, ".category_counts.tsv"), header=TRUE, sep="\t", quote="", check.names=FALSE)

cols_cat <- c("Upstream"="#ef6a63", "Downstream"="#7aa63f", "Gene body"="#2cb7b8", "Intergenic"="#2f9f8c")
bar_col <- "#5a9cad"

plot_one <- function(file, device=c("pdf","png","svg")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width=7.6, height=5.2, family="Times")
  if (device == "png") png(file, width=2200, height=1500, res=300, type="cairo")
  if (device == "svg") svg(file, width=7.6, height=5.2, family="Times")

  op <- par(no.readonly=TRUE)
  on.exit({par(op); dev.off()}, add=TRUE)

  counts <- dist$event_count
  n <- length(counts)
  ymax <- max(counts)
  lower_max <- 90000
  upper_min <- 160000
  upper_max <- max(180000, ceiling(ymax / 10000) * 10000)
  if (upper_max <= ymax) upper_max <- upper_max + 10000

  layout(matrix(c(1,2), nrow=2), heights=c(1.1, 2.6))
  par(mar=c(0.2, 4.4, 1.1, 0.8), xaxs="i", yaxs="i", family="Times")
  barplot(counts, col=bar_col, border=NA, ylim=c(upper_min, upper_max), axes=FALSE, space=0.15)
  top_ticks <- seq(upper_min, upper_max, by=10000)
  axis(2, at=top_ticks, labels=paste0(top_ticks/1000, "k"), las=1)
  box(bty="l")
  usr <- par("usr")
  segments(usr[1]-0.35, upper_min*0.99, usr[1]-0.15, upper_min*1.03, xpd=NA)
  segments(usr[1]-0.15, upper_min*0.99, usr[1]+0.05, upper_min*1.03, xpd=NA)

  par(fig=c(0.58, 0.93, 0.58, 0.96), new=TRUE, mar=c(0,0,0,0), family="Times")
  pie_counts <- setNames(catdf$count, catdf$category)
  pie_counts <- pie_counts[names(cols_cat)]
  pie_labels <- paste0(round(100 * pie_counts / sum(pie_counts)), "%")
  pie(pie_counts, labels=pie_labels, col=cols_cat[names(pie_counts)], border="white", cex=0.95)
  legend("bottom", legend=names(cols_cat), fill=cols_cat, bty="n", cex=0.72, ncol=2, inset=-0.18, xpd=NA)

  par(fig=c(0,1,0,0.66), new=TRUE)
  par(mar=c(4.2, 4.4, 0.5, 0.8), xaxs="i", yaxs="i", family="Times")
  bp <- barplot(counts, col=bar_col, border=NA, ylim=c(0, lower_max), axes=FALSE, space=0.15)
  yt <- seq(0, lower_max, by=10000)
  axis(2, at=yt, labels=ifelse(yt >= 1000, paste0(yt/1000, "k"), yt), las=1)
  major_labels <- rep("", n)
  major_labels[1] <- "<-5k"
  major_labels[which(dist$bin_label == "Gene")] <- "Gene"
  major_labels[n] <- ">5k"
  centers <- dist$bin_center_bp
  for (i in seq_along(centers)) {
    if (abs(centers[i]) < 5000 && centers[i] != 0 && abs(centers[i]) %% 2000 == 250) {
      major_labels[i] <- paste0(ifelse(centers[i] < 0, "-", ""), abs(round(centers[i]/1000)), "k")
    }
  }
  axis(1, at=bp, labels=major_labels, las=1, tick=TRUE)
  mtext("Gene", side=1, line=2.7, cex=1.2)
  mtext("PAV number", side=2, line=2.7, cex=1.2)
  box(bty="l")
  usr <- par("usr")
  segments(usr[1]-0.35, lower_max*1.00, usr[1]-0.15, lower_max*1.05, xpd=NA)
  segments(usr[1]-0.15, lower_max*1.00, usr[1]+0.05, lower_max*1.05, xpd=NA)
}

plot_one(paste0(out_prefix, ".pdf"), "pdf")
plot_one(paste0(out_prefix, ".png"), "png")
plot_one(paste0(out_prefix, ".svg"), "svg")
