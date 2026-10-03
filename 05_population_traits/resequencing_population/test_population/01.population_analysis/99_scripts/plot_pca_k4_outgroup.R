base <- "path/to/project/38.medicago_resequence/test_population/01.population_analysis"
pca_file <- file.path(base, "02_pca/sativa182.pca.tsv")
grp_file <- file.path(base, "06_summary/pca_group_assignment_K4.tsv")
eig_file <- file.path(base, "02_pca/sativa182.pca.eigenval")
out_prefix <- file.path(base, "05_plots/sativa182.PCA_K4_outgroup")
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

panel_plot <- function(x, y, xlab, ylab, panel="") {
  plot(d[[x]], d[[y]], type="n", xlab=xlab, ylab=ylab, bty="l", cex.lab=1.25, cex.axis=1.1)
  if (nzchar(panel)) mtext(panel, side=3, adj=-0.18, line=0.6, font=2, cex=1.5)
  for (g in levels(d$plot_group)) {
    idx <- which(d$plot_group == g)
    if (length(idx)) points(d[[x]][idx], d[[y]][idx], pch=pchs[g], col=cols[g], bg=cols[g], cex=0.95)
  }
  legend("right", legend=levels(d$plot_group), col=cols, pch=pchs, pt.cex=1.15, bty="n", cex=1.15, y.intersp=0.9)
}

pdf(paste0(out_prefix, ".PC1_PC2.pdf"), width=5.2, height=4.5)
par(mar=c(4.5,4.8,1.5,1))
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a")
dev.off()

pdf(paste0(out_prefix, ".PC1_PC3.pdf"), width=5.2, height=4.5)
par(mar=c(4.5,4.8,1.5,1))
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b")
dev.off()

pdf(paste0(out_prefix, ".combined.pdf"), width=10.4, height=4.5)
par(mfrow=c(1,2), mar=c(4.5,4.8,1.5,1), oma=c(0,0,0,0))
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a")
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b")
dev.off()

png(paste0(out_prefix, ".combined.png"), width=2600, height=1100, res=250)
par(mfrow=c(1,2), mar=c(4.5,4.8,1.5,1), oma=c(0,0,0,0))
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a")
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b")
dev.off()

png(paste0(out_prefix, ".PC1_PC2.png"), width=1300, height=1100, res=250)
par(mar=c(4.5,4.8,1.5,1))
panel_plot("PC1", "PC2", paste0("PC1 (", vars[1], "%)"), paste0("PC2 (", vars[2], "%)"), "a")
dev.off()

png(paste0(out_prefix, ".PC1_PC3.png"), width=1300, height=1100, res=250)
par(mar=c(4.5,4.8,1.5,1))
panel_plot("PC1", "PC3", paste0("PC1 (", vars[1], "%)"), paste0("PC3 (", vars[3], "%)"), "b")
dev.off()

cat("group counts:\n")
print(table(d$plot_group, useNA="ifany"))
cat("outputs:\n")
cat(paste0(out_prefix, c(".combined.png", ".combined.pdf", ".PC1_PC2.png", ".PC1_PC3.png"), collapse="\n"), "\n")
cat("assignment:\n", assign_out, "\n")
