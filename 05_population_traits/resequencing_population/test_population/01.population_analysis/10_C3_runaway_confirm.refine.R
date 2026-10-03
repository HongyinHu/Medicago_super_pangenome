args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript refine_C3_runaway_candidates.R <main_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
base_dir <- file.path(main_dir, "10_C3_runaway_confirm")
out_dir <- file.path(main_dir, "10_C3_runaway_confirm", "refined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required.")
}
library(ape)

parse_counts <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  vals <- setNames(rep(0L, 5), c("out", "C1", "C2", "C3", "C4"))
  for (p in parts) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    vals[kv[1]] <- as.integer(kv[2])
  }
  vals
}

clades <- read.delim(file.path(base_dir, "C3_top_pure_clades_with_sister_context.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
metrics <- read.delim(file.path(base_dir, "C3_runaway_candidate_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

is_local_c2_c3 <- vapply(seq_len(nrow(clades)), function(i) {
  pc <- parse_counts(clades$parent_counts[i])
  sc <- parse_counts(clades$sister_counts[i])
  pc["out"] == 0 && pc["C1"] == 0 && pc["C4"] == 0 &&
    pc["C2"] > 0 && pc["C3"] > 0 &&
    sc["C2"] > 0 && sc["out"] == 0 && sc["C1"] == 0 && sc["C4"] == 0
}, logical(1))

local_clades <- clades[is_local_c2_c3 & clades$n_C3 <= 12, ]
ml_tip_sets <- strsplit(local_clades$tips[local_clades$tree == "ML"], ",", fixed = TRUE)
nj_tip_sets <- strsplit(local_clades$tips[local_clades$tree == "NJ"], ",", fixed = TRUE)

shared_clade_ids <- character(0)
if (length(ml_tip_sets) > 0 && length(nj_tip_sets) > 0) {
  for (m in ml_tip_sets) {
    for (n in nj_tip_sets) {
      if (setequal(m, n)) {
        shared_clade_ids <- union(shared_clade_ids, m)
      }
    }
  }
}

singleton_ids <- metrics$IID[
  !(metrics$IID %in% shared_clade_ids) &
    metrics$nearest_group.ML == "C2" &
    metrics$nearest_group.NJ == "C2" &
    metrics$nearest_C2_rank.ML <= 2 &
    metrics$nearest_C2_rank.NJ <= 2
]

final_ids <- unique(c(shared_clade_ids, singleton_ids))
evidence <- data.frame(
  IID = final_ids,
  evidence = ifelse(final_ids %in% shared_clade_ids,
                    "C3_clade_sister_to_C2_in_both_ML_NJ",
                    "singleton_nearest_C2_in_both_ML_NJ"),
  stringsAsFactors = FALSE
)

final <- merge(evidence, metrics, by = "IID", all.x = TRUE)
final <- final[order(final$evidence, final$IID), ]

write.table(local_clades, file.path(out_dir, "local_C3_clades_sister_to_C2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final, file.path(out_dir, "C3_runaway_refined_confirmed_samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(final_ids, file.path(out_dir, "C3_runaway_refined_confirmed_ids.txt"))

group_file <- file.path(main_dir, "08_tree_K4_with_outgroup", "00_groups", "K4_with_outgroup.keep_with_group.tsv")
pca_file <- file.path(main_dir, "06_summary", "pca_group_assignment_K4_with_outgroup.tsv")
ml_file <- file.path(main_dir, "08_tree_K4_with_outgroup", "02_ml", "sativa182_pruned_varsites.K4_with_outgroup.ML.rooted_by_outgroup.treefile")
nj_file <- file.path(main_dir, "08_tree_K4_with_outgroup", "01_nj", "sativa182_K4_with_outgroup.pruned.NJ.rooted_by_outgroup.treefile")

groups <- read.delim(group_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(groups) <- c("FID", "IID", "plot_group")
groups$IID <- as.character(groups$IID)
groups$plot_group <- as.character(groups$plot_group)

pca <- read.delim(pca_file, stringsAsFactors = FALSE, check.names = FALSE)
pca$IID <- as.character(pca$IID)
for (cc in c("PC1", "PC2", "PC3")) pca[[cc]] <- as.numeric(pca[[cc]])
pca$plot_group <- groups$plot_group[match(pca$IID, groups$IID)]

palette <- c(out = "#111111", C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
shape <- c(out = 16, C1 = 16, C2 = 15, C3 = 18, C4 = 17)

plot_pca_highlight <- function(xpc, ypc, suffix) {
  png(file.path(out_dir, paste0("PCA_", suffix, "_refined_C3_runaway.png")), width = 1650, height = 1350, res = 220, type = "cairo")
  par(mar = c(4.2, 4.2, 1.2, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in names(palette)) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g], bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  idx <- match(final_ids, pca$IID)
  points(pca[[xpc]][idx], pca[[ypc]][idx], pch = 21, col = "#f0b000", bg = NA, cex = 1.65, lwd = 1.7)
  text(pca[[xpc]][idx], pca[[ypc]][idx], labels = final_ids, pos = 3, cex = 0.58, col = "#333333")
  legend("topright", legend = c(names(palette), "confirmed C3 outliers"),
         pch = c(shape, 21), col = c(palette, "#f0b000"),
         pt.bg = c(palette, NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, 5), 1.35))
  dev.off()

  pdf(file.path(out_dir, paste0("PCA_", suffix, "_refined_C3_runaway.pdf")), width = 6.2, height = 5.1, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.2, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in names(palette)) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g], bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  idx <- match(final_ids, pca$IID)
  points(pca[[xpc]][idx], pca[[ypc]][idx], pch = 21, col = "#f0b000", bg = NA, cex = 1.65, lwd = 1.7)
  text(pca[[xpc]][idx], pca[[ypc]][idx], labels = final_ids, pos = 3, cex = 0.58, col = "#333333")
  legend("topright", legend = c(names(palette), "confirmed C3 outliers"),
         pch = c(shape, 21), col = c(palette, "#f0b000"),
         pt.bg = c(palette, NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, 5), 1.35))
  dev.off()
}

plot_tree_highlight <- function(tree_file, title, prefix) {
  tree <- read.tree(tree_file)
  tip_group <- groups$plot_group[match(tree$tip.label, groups$IID)]
  names(tip_group) <- tree$tip.label
  tip_col <- palette[tip_group]
  tip_col[is.na(tip_col)] <- "#777777"
  is_final <- tree$tip.label %in% final_ids

  png(file.path(out_dir, paste0(prefix, "_refined_C3_runaway.png")), width = 1900, height = 1900, res = 300, type = "cairo")
  par(mar = c(0, 0, 1, 0))
  plot.phylo(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE, edge.width = 0.55, edge.color = "#303030")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index], pch = ifelse(is_final, 21, 16),
         col = ifelse(is_final, "#f0b000", tip_col), bg = ifelse(is_final, NA, tip_col),
         cex = ifelse(is_final, 0.95, 0.50), lwd = ifelse(is_final, 1.2, 0.6))
  cand_idx <- which(is_final)
  text(pp$xx[cand_idx], pp$yy[cand_idx], labels = tree$tip.label[cand_idx], cex = 0.48, pos = 3, col = "#222222")
  mtext(title, side = 3, line = -1.0, cex = 1.12, font = 2)
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_refined_C3_runaway.pdf")), width = 6.3, height = 6.3, useDingbats = FALSE)
  par(mar = c(0, 0, 1, 0))
  plot.phylo(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE, edge.width = 0.55, edge.color = "#303030")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index], pch = ifelse(is_final, 21, 16),
         col = ifelse(is_final, "#f0b000", tip_col), bg = ifelse(is_final, NA, tip_col),
         cex = ifelse(is_final, 0.95, 0.50), lwd = ifelse(is_final, 1.2, 0.6))
  cand_idx <- which(is_final)
  text(pp$xx[cand_idx], pp$yy[cand_idx], labels = tree$tip.label[cand_idx], cex = 0.48, pos = 3, col = "#222222")
  mtext(title, side = 3, line = -1.0, cex = 1.12, font = 2)
  dev.off()
}

plot_pca_highlight("PC1", "PC2", "PC1_PC2")
plot_pca_highlight("PC1", "PC3", "PC1_PC3")
plot_tree_highlight(ml_file, "ML tree refined C3 outliers", "ML_tree")
plot_tree_highlight(nj_file, "NJ tree refined C3 outliers", "NJ_tree")

cat("Refined confirmed IDs:\n")
cat(paste(final_ids, collapse = "\n"), "\n")
cat("Files:\n")
cat(file.path(out_dir, "C3_runaway_refined_confirmed_samples.tsv"), "\n")
cat(file.path(out_dir, "PCA_PC1_PC2_refined_C3_runaway.png"), "\n")
