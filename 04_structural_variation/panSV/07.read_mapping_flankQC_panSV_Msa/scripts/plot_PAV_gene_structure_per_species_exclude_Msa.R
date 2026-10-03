#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_PAV_gene_structure_per_species_exclude_Msa.R <per_species.tsv> <out_prefix>")
}
infile <- args[1]
out_prefix <- args[2]

df <- read.table(infile, header=TRUE, sep="\t", quote="", check.names=FALSE)
levels_cat <- c("Intergenic", "2kb Upstream", "2kb downstream", "Intron", "Exon")
df$category <- factor(df$category, levels=levels_cat)
df$proportion_percent <- as.numeric(df$proportion_percent)
cols <- c(
  "Intergenic"="#ef7670",
  "2kb Upstream"="#b3b21e",
  "2kb downstream"="#1ab88a",
  "Intron"="#00a8e0",
  "Exon"="#d957d9"
)

plot_one <- function(file, device=c("pdf","png","svg")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width=3.55, height=4.6, family="Times")
  if (device == "png") png(file, width=1065, height=1380, res=300, type="cairo")
  if (device == "svg") svg(file, width=3.55, height=4.6, family="Times")
  op <- par(no.readonly=TRUE)
  on.exit({par(op); dev.off()}, add=TRUE)

  set.seed(1)
  x <- seq_along(levels_cat)
  vals <- lapply(levels_cat, function(cat) df$proportion_percent[df$category == cat])
  names(vals) <- levels_cat

  par(mar=c(6.0, 4.8, 0.8, 0.6), xaxs="i", yaxs="i", family="Times")
  plot(NA, xlim=c(0.35, 5.65), ylim=c(0, 42), axes=FALSE, xlab="", ylab="")
  for (i in x) {
    y <- vals[[i]]
    boxplot(y, at=i, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.34, lwd=1.0)
    stripchart(y, at=i, add=TRUE, vertical=TRUE, method="jitter",
               jitter=0.08, pch=16, cex=0.58, col=adjustcolor(cols[i], alpha.f=0.85))
  }
  axis(2, at=seq(0, 40, by=5), las=1, cex.axis=0.95)
  axis(1, at=x, labels=FALSE, tick=FALSE)
  text(x, par("usr")[3] - 2.2, labels=levels_cat, srt=35, adj=1, xpd=NA, cex=0.95)
  mtext("Proportion (%)", side=2, line=3.1, cex=1.1)
  box(bty="l", lwd=1.1)
}

plot_one(paste0(out_prefix, ".pdf"), "pdf")
plot_one(paste0(out_prefix, ".png"), "png")
plot_one(paste0(out_prefix, ".svg"), "svg")
