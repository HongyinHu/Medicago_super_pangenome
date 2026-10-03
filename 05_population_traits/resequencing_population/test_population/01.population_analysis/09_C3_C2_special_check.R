args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript analyze_C3_C2_special_check.R <main_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- file.path(main_dir, "09_C3_C2_special_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required.")
}
library(ape)

pca_file <- file.path(main_dir, "06_summary", "pca_group_assignment_K4_with_outgroup.tsv")
group_file <- file.path(main_dir, "08_tree_K4_with_outgroup", "00_groups", "K4_with_outgroup.keep_with_group.tsv")
ml_file <- file.path(main_dir, "08_tree_K4_with_outgroup", "02_ml", "sativa182_pruned_varsites.K4_with_outgroup.ML.rooted_by_outgroup.treefile")
nj_file <- file.path(main_dir, "08_tree_K4_with_outgroup", "01_nj", "sativa182_K4_with_outgroup.pruned.NJ.rooted_by_outgroup.treefile")

pca <- read.delim(pca_file, stringsAsFactors = FALSE, check.names = FALSE)
groups <- read.delim(group_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(groups) <- c("FID", "IID", "plot_group")
groups$IID <- as.character(groups$IID)
pca$IID <- as.character(pca$IID)
pca <- merge(pca, groups[, c("IID", "plot_group")], by = "IID", suffixes = c("", ".tree"))
pca$plot_group <- pca$plot_group.tree
pca$plot_group.tree <- NULL

numeric_cols <- intersect(c("PC1", "PC2", "PC3"), names(pca))
for (cc in numeric_cols) pca[[cc]] <- as.numeric(pca[[cc]])

candidate_stats <- function(tree_file, tree_name, k = 5) {
  tree <- read.tree(tree_file)
  tree_groups <- pca$plot_group[match(tree$tip.label, pca$IID)]
  names(tree_groups) <- tree$tip.label
  dist_mat <- cophenetic.phylo(tree)
  c3_ids <- tree$tip.label[tree_groups == "C3"]

  rows <- lapply(c3_ids, function(id) {
    d <- dist_mat[id, ]
    d <- d[names(d) != id]
    ordered <- names(sort(d))
    nn <- ordered[seq_len(min(k, length(ordered)))]
    nn_groups <- tree_groups[nn]
    data.frame(
      IID = id,
      tree = tree_name,
      nearest_group = nn_groups[1],
      nearest_id = nn[1],
      nearest_dist = unname(d[nn[1]]),
      n_C2_in_5NN = sum(nn_groups == "C2", na.rm = TRUE),
      n_C3_in_5NN = sum(nn_groups == "C3", na.rm = TRUE),
      n_out_in_5NN = sum(nn_groups == "out", na.rm = TRUE),
      fiveNN_groups = paste(nn_groups, collapse = ","),
      fiveNN_ids = paste(nn, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ml_stats <- candidate_stats(ml_file, "ML", k = 5)
nj_stats <- candidate_stats(nj_file, "NJ", k = 5)
stats <- rbind(ml_stats, nj_stats)

centroids <- aggregate(cbind(PC1, PC2, PC3) ~ plot_group, data = pca[pca$plot_group %in% c("C1", "C2", "C3", "C4"), ], mean)
get_centroid <- function(group) as.numeric(centroids[centroids$plot_group == group, c("PC1", "PC2", "PC3")])
c2_centroid <- get_centroid("C2")
c3_centroid <- get_centroid("C3")

pca_rows <- pca[pca$plot_group == "C3", c("IID", "plot_group", "PC1", "PC2", "PC3")]
pca_rows$dist_to_C2_centroid_PC123 <- apply(pca_rows[, c("PC1", "PC2", "PC3")], 1, function(x) sqrt(sum((x - c2_centroid)^2)))
pca_rows$dist_to_C3_centroid_PC123 <- apply(pca_rows[, c("PC1", "PC2", "PC3")], 1, function(x) sqrt(sum((x - c3_centroid)^2)))
pca_rows$closer_centroid_PC123 <- ifelse(
  pca_rows$dist_to_C2_centroid_PC123 < pca_rows$dist_to_C3_centroid_PC123,
  "C2",
  "C3"
)

wide <- reshape(
  stats[, c("IID", "tree", "nearest_group", "nearest_id", "n_C2_in_5NN", "n_C3_in_5NN", "fiveNN_groups", "fiveNN_ids")],
  idvar = "IID",
  timevar = "tree",
  direction = "wide"
)
candidate <- merge(wide, pca_rows, by = "IID", all.x = TRUE)
candidate$tree_support_C2 <- rowSums(candidate[, c("n_C2_in_5NN.ML", "n_C2_in_5NN.NJ")] >= 3, na.rm = TRUE)
candidate$nearest_support_C2 <- rowSums(candidate[, c("nearest_group.ML", "nearest_group.NJ")] == "C2", na.rm = TRUE)
candidate$recommendation <- ifelse(
  candidate$tree_support_C2 >= 1 & candidate$closer_centroid_PC123 == "C2",
  "strong_reassign_to_C2",
  ifelse(candidate$tree_support_C2 >= 1 | candidate$nearest_support_C2 >= 1,
         "tree_intermediate_keep_review",
         "keep_C3")
)

candidate <- candidate[order(
  -candidate$tree_support_C2,
  -candidate$nearest_support_C2,
  candidate$dist_to_C2_centroid_PC123 - candidate$dist_to_C3_centroid_PC123
), ]

write.table(stats, file.path(out_dir, "C3_tree_5NN_stats.long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(candidate, file.path(out_dir, "C3_C2_special_candidates.with_PCA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

special <- candidate[candidate$tree_support_C2 >= 1 | candidate$nearest_support_C2 >= 1, ]
write.table(special, file.path(out_dir, "C3_tree_near_C2_special_samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

palette <- c(out = "#111111", C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
shape <- c(out = 16, C1 = 16, C2 = 15, C3 = 18, C4 = 17)

plot_pca <- function(xpc, ypc, suffix) {
  x <- pca[[xpc]]
  y <- pca[[ypc]]
  png(file.path(out_dir, paste0("PCA_", suffix, "_highlight_C3_near_C2.png")), width = 1600, height = 1350, res = 220, type = "cairo")
  par(mar = c(4.2, 4.2, 1.2, 1.0))
  plot(x, y, type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in names(palette)) {
    idx <- pca$plot_group == g
    points(x[idx], y[idx], pch = shape[g], col = palette[g], bg = palette[g], cex = ifelse(g == "out", 0.9, 0.85))
  }
  if (nrow(special) > 0) {
    idx <- match(special$IID, pca$IID)
    points(x[idx], y[idx], pch = 21, col = "#f0b000", bg = NA, cex = 1.75, lwd = 1.7)
    text(x[idx], y[idx], labels = special$IID, pos = 3, cex = 0.62, col = "#333333")
  }
  legend("topright", legend = c(names(palette), "C3 near C2 in tree"),
         pch = c(shape, 21), col = c(palette, "#f0b000"),
         pt.bg = c(palette, NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, 5), 1.35))
  dev.off()

  pdf(file.path(out_dir, paste0("PCA_", suffix, "_highlight_C3_near_C2.pdf")), width = 6.0, height = 5.1, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.2, 1.0))
  plot(x, y, type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in names(palette)) {
    idx <- pca$plot_group == g
    points(x[idx], y[idx], pch = shape[g], col = palette[g], bg = palette[g], cex = ifelse(g == "out", 0.9, 0.85))
  }
  if (nrow(special) > 0) {
    idx <- match(special$IID, pca$IID)
    points(x[idx], y[idx], pch = 21, col = "#f0b000", bg = NA, cex = 1.75, lwd = 1.7)
    text(x[idx], y[idx], labels = special$IID, pos = 3, cex = 0.62, col = "#333333")
  }
  legend("topright", legend = c(names(palette), "C3 near C2 in tree"),
         pch = c(shape, 21), col = c(palette, "#f0b000"),
         pt.bg = c(palette, NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, 5), 1.35))
  dev.off()
}

plot_pca("PC1", "PC2", "PC1_PC2")
plot_pca("PC1", "PC3", "PC1_PC3")

cat("Wrote:\n")
cat(file.path(out_dir, "C3_C2_special_candidates.with_PCA.tsv"), "\n")
cat(file.path(out_dir, "C3_tree_near_C2_special_samples.tsv"), "\n")
cat(file.path(out_dir, "PCA_PC1_PC2_highlight_C3_near_C2.png"), "\n")
cat(file.path(out_dir, "PCA_PC1_PC3_highlight_C3_near_C2.png"), "\n")
