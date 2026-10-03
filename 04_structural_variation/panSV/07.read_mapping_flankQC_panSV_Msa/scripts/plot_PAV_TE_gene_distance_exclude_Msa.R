#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_PAV_TE_gene_distance_exclude_Msa.R <events.tsv> <out_prefix>")
}
events_file <- args[1]
out_prefix <- args[2]

df <- read.table(events_file, header=TRUE, sep="\t", quote="", check.names=FALSE)
df <- df[df$TE_class_plot %in% c("Non-TE", "DNA", "LTR"), ]
df$TE_group <- factor(df$TE_class_plot, levels=c("Non-TE", "DNA", "LTR"),
                      labels=c("Non-TE PAV", "DNA-TE PAV", "LTR-TE PAV"))
df$y <- as.numeric(df$log2_distance_plus1)
cols <- c("Non-TE PAV"="#2f9f8c", "DNA-TE PAV"="#df8f37", "LTR-TE PAV"="#43b7ce")

groups <- levels(df$TE_group)
vals <- lapply(groups, function(g) df$y[df$TE_group == g])
names(vals) <- groups

pairs <- list(c("Non-TE PAV", "DNA-TE PAV"),
              c("DNA-TE PAV", "LTR-TE PAV"),
              c("Non-TE PAV", "LTR-TE PAV"))
test_rows <- lapply(pairs, function(p) {
  wt <- wilcox.test(vals[[p[1]]], vals[[p[2]]], exact=FALSE)
  data.frame(group1=p[1], group2=p[2], p_value=wt$p.value, stringsAsFactors=FALSE)
})
tests <- do.call(rbind, test_rows)
write.table(tests, paste0(out_prefix, ".wilcoxon.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

format_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 2.2e-16) return("P < 2.2e-16")
  paste0("P = ", format(p, scientific=TRUE, digits=2))
}

draw_violin <- function(x, at, col, width=0.43, ymax=22.5) {
  if (length(x) < 2) return()
  d <- density(x, from=0, to=ymax, n=512, adjust=0.9)
  dx <- d$y / max(d$y) * width
  polygon(c(at - dx, rev(at + dx)), c(d$x, rev(d$x)),
          col=adjustcolor(col, alpha.f=0.18), border=col, lwd=1.4)
  boxplot(x, at=at, add=TRUE, outline=FALSE, boxwex=0.18, axes=FALSE,
          col=adjustcolor("white", alpha.f=0.85), border=col, whiskcol=col,
          staplecol=col, medcol=col, lwd=1.0)
  points(rep(at, min(length(x), 900)), sample(x, min(length(x), 900)),
         pch=16, col=adjustcolor(col, alpha.f=0.14), cex=0.18)
}

plot_one <- function(file, device=c("pdf", "png", "svg")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width=4.2, height=4.6, family="Times")
  if (device == "png") png(file, width=1260, height=1380, res=300, type="cairo")
  if (device == "svg") svg(file, width=4.2, height=4.6, family="Times")
  op <- par(no.readonly=TRUE)
  on.exit({par(op); dev.off()}, add=TRUE)

  ymax <- 25
  par(mar=c(5.8, 4.8, 1.0, 0.8), xaxs="i", yaxs="i", family="Times")
  plot(NA, xlim=c(0.45, 3.55), ylim=c(-0.3, ymax), axes=FALSE, xlab="", ylab="")
  for (i in seq_along(groups)) {
    draw_violin(vals[[groups[i]]], i, cols[groups[i]], ymax=ymax)
  }
  axis(2, at=seq(0, 20, by=5), las=1, cex.axis=0.95)
  axis(1, at=1:3, labels=FALSE, tick=FALSE)
  text(1:3, par("usr")[3] - 0.75, labels=groups, srt=50, adj=1, xpd=NA, cex=1.0)
  mtext(expression(log[2]*"(distance to gene (bp) + 1)"), side=2, line=3.0, cex=1.05)
  box(bty="l", lwd=1.1)

  # Pairwise brackets.
  bracket <- function(x1, x2, y, label) {
    segments(x1, y, x2, y, lwd=0.9)
    segments(x1, y, x1, y-0.22, lwd=0.9)
    segments(x2, y, x2, y-0.22, lwd=0.9)
    text((x1+x2)/2, y+0.28, label, cex=0.78)
  }
  bracket(1, 2, 20.7, format_p(tests$p_value[tests$group1=="Non-TE PAV" & tests$group2=="DNA-TE PAV"]))
  bracket(2, 3, 22.0, format_p(tests$p_value[tests$group1=="DNA-TE PAV" & tests$group2=="LTR-TE PAV"]))
  bracket(1, 3, 23.35, format_p(tests$p_value[tests$group1=="Non-TE PAV" & tests$group2=="LTR-TE PAV"]))
}

plot_one(paste0(out_prefix, ".pdf"), "pdf")
plot_one(paste0(out_prefix, ".png"), "png")
plot_one(paste0(out_prefix, ".svg"), "svg")
