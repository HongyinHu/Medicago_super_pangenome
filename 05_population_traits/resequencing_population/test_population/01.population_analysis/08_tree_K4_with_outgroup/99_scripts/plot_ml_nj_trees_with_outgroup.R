args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_ml_nj_trees_with_outgroup.R <main_dir> <out_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- normalizePath(args[2], mustWork = TRUE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required for tree reading, NJ construction, and plotting.")
}
library(ape)

group_file <- file.path(out_dir, "00_groups", "K4_with_outgroup.keep_with_group.tsv")
groups <- read.delim(group_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(groups) <- c("FID", "IID", "group")
groups$IID <- as.character(groups$IID)
groups$group <- factor(groups$group, levels = c("out", "C1", "C2", "C3", "C4"))

palette <- c(out = "#111111", C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
group_counts <- table(groups$group)
group_counts <- group_counts[group_counts > 0]
legend_labels <- paste0(names(group_counts), " (n=", as.integer(group_counts), ")")

drop_missing_groups <- function(tree) {
  drop <- setdiff(tree$tip.label, groups$IID)
  if (length(drop) > 0) tree <- drop.tip(tree, drop)
  tree
}

root_with_outgroup <- function(tree) {
  out_tips <- intersect(tree$tip.label, groups$IID[as.character(groups$group) == "out"])
  if (length(out_tips) == 0) return(tree)
  tryCatch(
    root(tree, outgroup = out_tips, resolve.root = TRUE),
    error = function(e) tree
  )
}

get_groups_for_tree <- function(tree) {
  as.character(groups$group[match(tree$tip.label, groups$IID)])
}

ml_file <- file.path(main_dir, "04_tree", "iqtree_repair", "sativa182_pruned_varsites.treefile")
ml_tree <- read.tree(ml_file)
ml_tree <- drop_missing_groups(ml_tree)
ml_tree <- root_with_outgroup(ml_tree)
write.tree(ml_tree, file = file.path(out_dir, "02_ml", "sativa182_pruned_varsites.K4_with_outgroup.ML.rooted_by_outgroup.treefile"))

dist_prefix <- file.path(out_dir, "01_nj", "sativa182_K4_with_outgroup.pruned")
ids <- read.table(paste0(dist_prefix, ".mdist.id"), stringsAsFactors = FALSE)
if (ncol(ids) >= 2) {
  sample_ids <- as.character(ids[[2]])
} else {
  sample_ids <- as.character(ids[[1]])
}
dist_df <- read.table(paste0(dist_prefix, ".mdist"), check.names = FALSE)
dist_mat <- as.matrix(dist_df)
if (nrow(dist_mat) != length(sample_ids) || ncol(dist_mat) != length(sample_ids)) {
  stop("PLINK distance matrix dimensions do not match .mdist.id.")
}
mode(dist_mat) <- "numeric"
rownames(dist_mat) <- sample_ids
colnames(dist_mat) <- sample_ids
nj_tree <- nj(as.dist(dist_mat))
nj_tree <- root_with_outgroup(nj_tree)
write.tree(nj_tree, file = file.path(out_dir, "01_nj", "sativa182_K4_with_outgroup.pruned.NJ.rooted_by_outgroup.treefile"))

plot_one_tree <- function(tree, main_title) {
  tip_group <- get_groups_for_tree(tree)
  tip_col <- palette[tip_group]
  tip_col[is.na(tip_col)] <- "#777777"
  plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    no.margin = TRUE,
    edge.width = 0.55,
    edge.color = "#303030"
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index], pch = 16, cex = 0.54, col = tip_col)
  mtext(main_title, side = 3, line = -1.0, cex = 1.15, font = 2)
}

plot_single <- function(tree, title, prefix) {
  png(paste0(prefix, ".png"), width = 1900, height = 1900, res = 300, type = "cairo")
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.82, pt.cex = 1.25)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 6.3, height = 6.3, useDingbats = FALSE)
  par(mar = c(0, 0, 1, 0))
  plot_one_tree(tree, title)
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.82, pt.cex = 1.25)
  dev.off()
}

plot_combined <- function(prefix) {
  png(paste0(prefix, ".png"), width = 3100, height = 1650, res = 300, type = "cairo")
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree with outgroup")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.68, pt.cex = 1.12)
  plot_one_tree(nj_tree, "NJ tree with outgroup")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.68, pt.cex = 1.12)
  dev.off()

  pdf(paste0(prefix, ".pdf"), width = 10.4, height = 5.5, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), oma = c(0, 0, 0, 0))
  plot_one_tree(ml_tree, "ML tree with outgroup")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.68, pt.cex = 1.12)
  plot_one_tree(nj_tree, "NJ tree with outgroup")
  legend("bottomleft", legend = legend_labels, pch = 16, col = palette[names(group_counts)],
         bty = "n", cex = 0.68, pt.cex = 1.12)
  dev.off()
}

plot_dir <- file.path(out_dir, "03_plots")
plot_single(ml_tree, "ML tree with outgroup", file.path(plot_dir, "sativa182_K4_with_outgroup.ML_tree"))
plot_single(nj_tree, "NJ tree with outgroup", file.path(plot_dir, "sativa182_K4_with_outgroup.NJ_tree"))
plot_combined(file.path(plot_dir, "sativa182_K4_with_outgroup.ML_NJ_tree.combined"))

writeLines(c(
  file.path(out_dir, "02_ml", "sativa182_pruned_varsites.K4_with_outgroup.ML.rooted_by_outgroup.treefile"),
  file.path(out_dir, "01_nj", "sativa182_K4_with_outgroup.pruned.NJ.rooted_by_outgroup.treefile"),
  file.path(plot_dir, "sativa182_K4_with_outgroup.ML_tree.png"),
  file.path(plot_dir, "sativa182_K4_with_outgroup.NJ_tree.png"),
  file.path(plot_dir, "sativa182_K4_with_outgroup.ML_NJ_tree.combined.png")
))
