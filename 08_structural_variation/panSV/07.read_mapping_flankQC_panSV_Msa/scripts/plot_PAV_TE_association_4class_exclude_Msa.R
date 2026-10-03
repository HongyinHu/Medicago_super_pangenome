#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_PAV_TE_association_4class_exclude_Msa.R <plot4class.tsv> <out_prefix>")
}
infile <- args[1]
out_prefix <- args[2]

df <- read.table(infile, header=TRUE, sep="\t", quote="", check.names=FALSE, stringsAsFactors=FALSE)
df$SVTYPE <- factor(df$SVTYPE, levels=c("DEL", "INS"), labels=c("Deletion", "Insertion"))
df$TE_class <- factor(df$TE_class, levels=c("DNA-TE", "LTR-TE", "Other TE", "Non-TE"))
cols <- c("DNA-TE"="#40507E", "LTR-TE"="#55BECA", "Other TE"="#E69A35", "Non-TE"="#029E73")

plot_one <- function(file, device=c("pdf","png","svg")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width=3.8, height=3.4, family="Times")
  if (device == "png") png(file, width=1140, height=1020, res=300, type="cairo")
  if (device == "svg") svg(file, width=3.8, height=3.4, family="Times")
  op <- par(no.readonly=TRUE)
  on.exit({par(op); dev.off()}, add=TRUE)
  par(mar=c(3.0, 3.9, 0.6, 4.35), xaxs="i", yaxs="i", family="Times")

  mat <- xtabs(percent_of_svtype ~ TE_class + SVTYPE, df)
  mat <- mat[c("DNA-TE", "LTR-TE", "Other TE", "Non-TE"), c("Deletion", "Insertion")]
  bp0 <- barplot(mat, plot=FALSE, width=0.55, space=c(0.55,0.65))
  bp <- barplot(
    mat, col=cols[rownames(mat)], border=NA, ylim=c(0,100), axes=FALSE,
    width=0.55, space=c(0.55,0.65), names.arg=rep("", ncol(mat)),
    xlim=c(0, max(bp0) + 1.55)
  )
  axis(2, at=seq(0,100,25), las=1, cex.axis=0.9)
  axis(1, at=bp, labels=colnames(mat), las=1, tick=FALSE, cex.axis=0.95, line=0.15)
  mtext("Percent (%)", side=2, line=2.5, cex=1.0)
  box(bty="l", lwd=1)
  legend(max(bp) + 0.36, 71, legend=rownames(mat), fill=cols[rownames(mat)],
         title="Type of PAVs", bty="n", cex=0.78, title.cex=0.82, xpd=NA)
}

for (ext in c("pdf", "png", "svg")) {
  out <- paste0(out_prefix, ".", ext)
  if (file.exists(out)) unlink(out)
  plot_one(out, ext)
  cat("WROTE\t", out, "\n", sep="")
}
