args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript redraw_revised_C3_to_C2_plain.R <out_dir>")
}

out_dir <- normalizePath(args[1], mustWork = TRUE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required.")
}
library(ape)

pca <- read.delim(file.path(out_dir, "00_groups", "pca_group_assignment_K4_revised_C3_to_C2_with_outgroup.tsv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
pca$IID <- as.character(pca$IID)
for (cc in c("PC1", "PC2", "PC3")) pca[[cc]] <- as.numeric(pca[[cc]])

counts <- read.delim(file.path(out_dir, "00_groups", "revised_group_counts.tsv"), stringsAsFactors = FALSE)
palette <- c(out = "#111111", C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
shape <- c(out = 16, C1 = 16, C2 = 15, C3 = 18, C4 = 17)
group_order <- c("out", "C1", "C2", "C3", "C4")
legend_labels <- paste0(counts$group, " (n=", counts$n, ")")

draw_pca <- function(xpc, ypc, prefix) {
  png(paste0(prefix, ".png"), width = 1500, height = 1250, res = 220, type = "cairo")
  par(mar = c(4.2, 4.2, 1.1, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in group_order) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g],
           bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  legend("topleft", legend = legend_labels,
         pch = shape[counts$group], col = palette[counts$group],
         pt.bg = palette[counts$group], bty = "n", cex = 0.82, pt.cex = 1.0)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 6.0, height = 5.1, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.1, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in group_order) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g],
           bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  legend("topleft", legend = legend_labels,
         pch = shape[counts$group], col = palette[counts$group],
         pt.bg = palette[counts$group], bty = "n", cex = 0.82, pt.cex = 1.0)
  dev.off()
}

draw_pca("PC1", "PC2", file.path(out_dir, "01_pca", "sativa182_revised_C3_to_C2.PCA_PC1_PC2"))
draw_pca("PC1", "PC3", file.path(out_dir, "01_pca", "sativa182_revised_C3_to_C2.PCA_PC1_PC3"))

ml_tree <- read.tree(file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.ML.rooted_by_outgroup.treefile"))
nj_tree <- read.tree(file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.NJ.rooted_by_outgroup.treefile"))

plot_one_tree <- function(tree, title) {
  tip_group <- pca$plot_group[match(tree$tip.label, pca$IID)]
  tip_col <- palette[tip_group]
  tip_col[is.na(tip_col)] <- "#777777"
  plot.phylo(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE,
             edge.width = 0.55, edge.color = "#303030")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index], pch = 16, col = tip_col, cex = 0.50)
  mtext(title, side = 3, line = -1.0, cex = 1.15, font = 2)
}

plot_single_tree <- function(tree, title, prefix) {
  png(paste0(prefix, ".png"), width = 1900, height = 1900, res = 300, type = "cairo")
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[counts$group],
         bty = "n", cex = 0.72, pt.cex = 1.1)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 6.3, height = 6.3, useDingbats = FALSE)
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[counts$group],
         bty = "n", cex = 0.72, pt.cex = 1.1)
  dev.off()
}

plot_combined <- function(prefix) {
  png(paste0(prefix, ".png"), width = 3100, height = 1650, res = 300, type = "cairo")
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree revised")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[counts$group],
         bty = "n", cex = 0.62, pt.cex = 1.0)
  plot_one_tree(nj_tree, "NJ tree revised")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[counts$group],
         bty = "n", cex = 0.62, pt.cex = 1.0)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 10.4, height = 5.5, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree revised")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[counts$group],
         bty = "n", cex = 0.62, pt.cex = 1.0)
  plot_one_tree(nj_tree, "NJ tree revised")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[counts$group],
         bty = "n", cex = 0.62, pt.cex = 1.0)
  dev.off()
}

plot_single_tree(ml_tree, "ML tree revised", file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.ML_tree"))
plot_single_tree(nj_tree, "NJ tree revised", file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.NJ_tree"))
plot_combined(file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.ML_NJ_tree.combined"))
