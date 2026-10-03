args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript revised_groups_pca_trees.R <main_dir> <out_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- normalizePath(args[2], mustWork = TRUE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required.")
}
library(ape)

reassign_ids <- readLines(file.path(out_dir, "00_groups", "reassigned_C3_to_C2.ids"))
reassign_ids <- reassign_ids[nzchar(reassign_ids)]

pca_file <- file.path(main_dir, "06_summary", "pca_group_assignment_K4_with_outgroup.tsv")
pca <- read.delim(pca_file, stringsAsFactors = FALSE, check.names = FALSE)
pca$IID <- as.character(pca$IID)
for (cc in intersect(c("PC1", "PC2", "PC3"), names(pca))) pca[[cc]] <- as.numeric(pca[[cc]])
pca$old_plot_group <- pca$plot_group
pca$plot_group <- ifelse(pca$IID %in% reassign_ids, "C2", pca$plot_group)
pca$reassigned_from_C3_to_C2 <- pca$IID %in% reassign_ids

write.table(pca, file.path(out_dir, "00_groups", "pca_group_assignment_K4_revised_C3_to_C2_with_outgroup.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

counts <- as.data.frame(table(factor(pca$plot_group, levels = c("out", "C1", "C2", "C3", "C4"))))
colnames(counts) <- c("group", "n")
counts <- counts[counts$n > 0, ]
write.table(counts, file.path(out_dir, "00_groups", "revised_group_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (g in c("C1", "C2", "C3", "C4")) {
  ids <- pca$IID[pca$plot_group == g]
  writeLines(ids, file.path(out_dir, "00_groups", paste0(g, ".keep")))
}

palette <- c(out = "#111111", C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
shape <- c(out = 16, C1 = 16, C2 = 15, C3 = 18, C4 = 17)
group_order <- c("out", "C1", "C2", "C3", "C4")
legend_labels <- paste0(counts$group, " (n=", counts$n, ")")

plot_pca_panel <- function(xpc, ypc, prefix) {
  png(paste0(prefix, ".png"), width = 1500, height = 1250, res = 220, type = "cairo")
  par(mar = c(4.2, 4.2, 1.1, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in group_order) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g],
           bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  idx_re <- pca$reassigned_from_C3_to_C2
  points(pca[[xpc]][idx_re], pca[[ypc]][idx_re], pch = 21, col = "#f0b000", bg = NA, cex = 1.65, lwd = 1.7)
  legend("topright", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(shape[counts$group], 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, nrow(counts)), 1.35))
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 6.0, height = 5.1, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.1, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in group_order) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g],
           bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  idx_re <- pca$reassigned_from_C3_to_C2
  points(pca[[xpc]][idx_re], pca[[ypc]][idx_re], pch = 21, col = "#f0b000", bg = NA, cex = 1.65, lwd = 1.7)
  legend("topright", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(shape[counts$group], 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, nrow(counts)), 1.35))
  dev.off()
}

plot_pca_panel("PC1", "PC2", file.path(out_dir, "01_pca", "sativa182_revised_C3_to_C2.PCA_PC1_PC2"))
plot_pca_panel("PC1", "PC3", file.path(out_dir, "01_pca", "sativa182_revised_C3_to_C2.PCA_PC1_PC3"))

ml_tree <- read.tree(file.path(main_dir, "08_tree_K4_with_outgroup", "02_ml", "sativa182_pruned_varsites.K4_with_outgroup.ML.rooted_by_outgroup.treefile"))
nj_tree <- read.tree(file.path(main_dir, "08_tree_K4_with_outgroup", "01_nj", "sativa182_K4_with_outgroup.pruned.NJ.rooted_by_outgroup.treefile"))
write.tree(ml_tree, file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.ML.rooted_by_outgroup.treefile"))
write.tree(nj_tree, file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.NJ.rooted_by_outgroup.treefile"))

plot_one_tree <- function(tree, title) {
  tip_group <- pca$plot_group[match(tree$tip.label, pca$IID)]
  is_reassigned <- tree$tip.label %in% reassign_ids
  tip_col <- palette[tip_group]
  tip_col[is.na(tip_col)] <- "#777777"
  plot.phylo(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE,
             edge.width = 0.55, edge.color = "#303030")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index],
         pch = ifelse(is_reassigned, 21, 16),
         col = ifelse(is_reassigned, "#f0b000", tip_col),
         bg = ifelse(is_reassigned, NA, tip_col),
         cex = ifelse(is_reassigned, 0.86, 0.50),
         lwd = ifelse(is_reassigned, 1.1, 0.6))
  mtext(title, side = 3, line = -1.0, cex = 1.15, font = 2)
}

plot_single_tree <- function(tree, title, prefix) {
  png(paste0(prefix, ".png"), width = 1900, height = 1900, res = 300, type = "cairo")
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(rep(16, nrow(counts)), 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.72, pt.cex = c(rep(1.1, nrow(counts)), 1.2))
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 6.3, height = 6.3, useDingbats = FALSE)
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(rep(16, nrow(counts)), 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.72, pt.cex = c(rep(1.1, nrow(counts)), 1.2))
  dev.off()
}

plot_combined_trees <- function(prefix) {
  png(paste0(prefix, ".png"), width = 3100, height = 1650, res = 300, type = "cairo")
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree revised")
  legend("bottomleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(rep(16, nrow(counts)), 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.62, pt.cex = c(rep(1.0, nrow(counts)), 1.1))
  plot_one_tree(nj_tree, "NJ tree revised")
  legend("bottomleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(rep(16, nrow(counts)), 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.62, pt.cex = c(rep(1.0, nrow(counts)), 1.1))
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 10.4, height = 5.5, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree revised")
  legend("bottomleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(rep(16, nrow(counts)), 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.62, pt.cex = c(rep(1.0, nrow(counts)), 1.1))
  plot_one_tree(nj_tree, "NJ tree revised")
  legend("bottomleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(rep(16, nrow(counts)), 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.62, pt.cex = c(rep(1.0, nrow(counts)), 1.1))
  dev.off()
}

plot_single_tree(ml_tree, "ML tree revised", file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.ML_tree"))
plot_single_tree(nj_tree, "NJ tree revised", file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.NJ_tree"))
plot_combined_trees(file.path(out_dir, "02_trees", "sativa182_revised_C3_to_C2.ML_NJ_tree.combined"))

cat("Revised group counts:\n")
print(counts)
