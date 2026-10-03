args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_ml_nj_trees.R <main_dir> <out_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- normalizePath(args[2], mustWork = TRUE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required for tree reading, NJ construction, and plotting.")
}
library(ape)

group_file <- file.path(out_dir, "00_groups", "K4_noout.keep_with_group.tsv")
groups <- read.delim(group_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(groups) <- c("FID", "IID", "group")
groups$IID <- as.character(groups$IID)
groups$group <- factor(groups$group, levels = c("C1", "C2", "C3", "C4"))

palette <- c(C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
group_counts <- table(groups$group)
legend_labels <- paste0(names(group_counts), " (n=", as.integer(group_counts), ")")

get_groups_for_tree <- function(tree) {
  as.character(groups$group[match(tree$tip.label, groups$IID)])
}

drop_non_k4 <- function(tree) {
  drop <- setdiff(tree$tip.label, groups$IID)
  if (length(drop) > 0) tree <- drop.tip(tree, drop)
  tree
}

ml_file <- file.path(main_dir, "04_tree", "iqtree_repair", "sativa182_pruned_varsites.treefile")
ml_tree <- read.tree(ml_file)
ml_tree <- drop_non_k4(ml_tree)
write.tree(ml_tree, file = file.path(out_dir, "02_ml", "sativa182_pruned_varsites.K4_noout.ML.treefile"))

dist_prefix <- file.path(out_dir, "01_nj", "sativa182_K4_noout.pruned")
ids <- read.table(paste0(dist_prefix, ".mdist.id"), stringsAsFactors = FALSE)
if (ncol(ids) >= 2) {
  sample_ids <- as.character(ids[[2]])
} else {
  sample_ids <- as.character(ids[[1]])
}
dist_df <- read.table(paste0(dist_prefix, ".mdist"), check.names = FALSE)
dist_mat <- as.matrix(dist_df)
if (nrow(dist_mat) != length(sample_ids) || ncol(dist_mat) != length(sample_ids)) {
  stop("PLINK distance matrix dimensions do not match .dist.id.")
}
mode(dist_mat) <- "numeric"
rownames(dist_mat) <- sample_ids
colnames(dist_mat) <- sample_ids
nj_tree <- nj(as.dist(dist_mat))
write.tree(nj_tree, file = file.path(out_dir, "01_nj", "sativa182_K4_noout.pruned.NJ.treefile"))

plot_one_tree <- function(tree, main_title) {
  tip_group <- get_groups_for_tree(tree)
  tip_col <- palette[tip_group]
  tip_col[is.na(tip_col)] <- "#444444"
  plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    no.margin = TRUE,
    edge.width = 0.55,
    edge.color = "#303030"
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  points(pp$xx[seq_len(Ntip(tree))], pp$yy[seq_len(Ntip(tree))],
         pch = 16, cex = 0.55, col = tip_col)
  mtext(main_title, side = 3, line = -1.0, cex = 1.15, font = 2)
}

plot_single <- function(tree, title, prefix) {
  png(paste0(prefix, ".png"), width = 1800, height = 1800, res = 300, type = "cairo")
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.85, pt.cex = 1.3)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 6.2, height = 6.2, useDingbats = FALSE)
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.85, pt.cex = 1.3)
  dev.off()
}

plot_combined <- function(prefix) {
  png(paste0(prefix, ".png"), width = 3000, height = 1600, res = 300, type = "cairo")
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.70, pt.cex = 1.15)
  plot_one_tree(nj_tree, "NJ tree")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.70, pt.cex = 1.15)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 10.2, height = 5.4, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.70, pt.cex = 1.15)
  plot_one_tree(nj_tree, "NJ tree")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.70, pt.cex = 1.15)
  dev.off()
}

plot_dir <- file.path(out_dir, "03_plots")
plot_single(ml_tree, "ML tree", file.path(plot_dir, "sativa182_K4_noout.ML_tree"))
plot_single(nj_tree, "NJ tree", file.path(plot_dir, "sativa182_K4_noout.NJ_tree"))
plot_combined(file.path(plot_dir, "sativa182_K4_noout.ML_NJ_tree.combined"))

writeLines(c(
  file.path(out_dir, "02_ml", "sativa182_pruned_varsites.K4_noout.ML.treefile"),
  file.path(out_dir, "01_nj", "sativa182_K4_noout.pruned.NJ.treefile"),
  file.path(plot_dir, "sativa182_K4_noout.ML_tree.png"),
  file.path(plot_dir, "sativa182_K4_noout.NJ_tree.png"),
  file.path(plot_dir, "sativa182_K4_noout.ML_NJ_tree.combined.png")
))
