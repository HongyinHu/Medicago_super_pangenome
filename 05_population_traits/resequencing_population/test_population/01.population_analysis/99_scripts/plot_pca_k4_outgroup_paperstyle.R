base <- "path/to/project/38.medicago_resequence/test_population/01.population_analysis"
pca_file <- file.path(base, "02_pca/sativa182.pca.tsv")
grp_file <- file.path(base, "06_summary/pca_group_assignment_K4.tsv")
eig_file <- file.path(base, "02_pca/sativa182.pca.eigenval")
out_prefix <- file.path(base, "05_plots/sativa182.PCA_K4_outgroup.paperstyle")
assign_out <- file.path(base, "06_summary/pca_group_assignment_K4_with_outgroup.tsv")

pca <- read.table(pca_file, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
grp <- read.table(grp_file, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
d <- merge(pca, grp[, c("sample", "pca_group")], by.x="IID", by.y="sample", all.x=TRUE, sort=FALSE)
d$plot_group <- sub("PCA_Group", "C", d$pca_group)
d$plot_group[grepl("^(SRR|DRR)", d$IID)] <- "out"
d$plot_group <- factor(d$plot_group, levels=c("out", "C1", "C2", "C3", "C4"))
write.table(d[, c("IID", "plot_group", "pca_group", "PC1", "PC2", "PC3")], assign_out, sep="\t", quote=FALSE, row.names=FALSE)

vars <- c(NA, NA, NA)
if (file.exists(eig_file)) {
  ev <- scan(eig_file, quiet=TRUE)
  vars <- round(ev[1:3] / sum(ev) * 100, 2)
}
cols <- c(out="black", C1="#ef6a6a", C2="#405a93", C3="#159c89", C4="#b786bd")
pchs <- c(out=16, C1=16, C2=15, C3=18, C4=17)

lims <- function(x) {
  r <- range(x, na.rm=TRUE); pad <- diff(r) * 0.10
  c(r[1]-pad, r[2]+pad)
}
xlim <- lims(d$PC1)
ylim12 <- lims(d$PC2)
ylim13 <- lims(d$PC3)

panel_plot <- function(x, y, xlab, ylab, panel, xlim, ylim) {
  plot(d[[x]], d[[y]], type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
       bty="l", cex.lab=1.05, cex.axis=0.9, las=1, xaxs="i", yaxs="i")
  mtext(panel, side=3, adj=-0.18, line=0.25, font=2, cex=1.25)
  for (g in levels(d$plot_group)) {
    idx <- which(d$plot_group == g)
    if (length(idx)) points(d[[x]][idx], d[[y]][idx], pch=pchs[g], col=cols[g], bg=cols[g], cex=0.78)
  }
  legend("right", inset=0.01, legend=levels(d$plot_group), col=cols, pch=pchs,
         pt.cex=0.95, bty="n", cex=0.9, y.intersp=0.82, x.intersp=0.7)
}

pdf(paste0(out_prefix, ".combined.pdf"), width=8.6, height=3.8)
par(mfrow=c(1,2), mar=c(3.7,4.1,1.0,0.5), mgp=c(2.35,0.65,0), tcl=-0.25)
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a", xlim, ylim12)
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b", xlim, ylim13)
dev.off()

png(paste0(out_prefix, ".combined.png"), width=2400, height=1050, res=300)
par(mfrow=c(1,2), mar=c(3.7,4.1,1.0,0.5), mgp=c(2.35,0.65,0), tcl=-0.25)
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a", xlim, ylim12)
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b", xlim, ylim13)
dev.off()

pdf(paste0(out_prefix, ".PC1_PC2.pdf"), width=4.3, height=3.8)
par(mar=c(3.7,4.1,1.0,0.5), mgp=c(2.35,0.65,0), tcl=-0.25)
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a", xlim, ylim12)
dev.off()

pdf(paste0(out_prefix, ".PC1_PC3.pdf"), width=4.3, height=3.8)
par(mar=c(3.7,4.1,1.0,0.5), mgp=c(2.35,0.65,0), tcl=-0.25)
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b", xlim, ylim13)
dev.off()

cat("group counts:\n")
print(table(d$plot_group, useNA="ifany"))
