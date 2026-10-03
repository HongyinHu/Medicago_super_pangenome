args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript confirm_C3_runaway_in_C2_region.R <main_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- file.path(main_dir, "10_C3_runaway_confirm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required.")
}
library(ape)

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
for (cc in intersect(c("PC1", "PC2", "PC3"), names(pca))) pca[[cc]] <- as.numeric(pca[[cc]])
pca$plot_group <- groups$plot_group[match(pca$IID, groups$IID)]

palette <- c(out = "#111111", C1 = "#ef6a6a", C2 = "#405a93", C3 = "#159c89", C4 = "#b786bd")
shape <- c(out = 16, C1 = 16, C2 = 15, C3 = 18, C4 = 17)

make_tree_helpers <- function(tree) {
  n_tip <- length(tree$tip.label)
  children <- split(tree$edge[, 2], tree$edge[, 1])
  parent <- setNames(tree$edge[, 1], tree$edge[, 2])

  tips_under <- function(node) {
    if (node <= n_tip) return(tree$tip.label[node])
    ch <- children[[as.character(node)]]
    unlist(lapply(ch, tips_under), use.names = FALSE)
  }

  list(n_tip = n_tip, children = children, parent = parent, tips_under = tips_under)
}

group_counts_text <- function(tips) {
  gg <- groups$plot_group[match(tips, groups$IID)]
  tt <- table(factor(gg, levels = c("out", "C1", "C2", "C3", "C4")))
  paste(paste0(names(tt), "=", as.integer(tt)), collapse = ";")
}

tree_tip_metrics <- function(tree_file, tree_name) {
  tree <- read.tree(tree_file)
  tip_group <- groups$plot_group[match(tree$tip.label, groups$IID)]
  names(tip_group) <- tree$tip.label
  dist_mat <- cophenetic.phylo(tree)
  c3_ids <- tree$tip.label[tip_group == "C3"]

  rows <- lapply(c3_ids, function(id) {
    d <- dist_mat[id, ]
    d <- d[names(d) != id]
    ordered <- names(sort(d))
    ordered_group <- tip_group[ordered]
    c2_rank <- which(ordered_group == "C2")[1]
    c2_id <- ordered[c2_rank]
    data.frame(
      IID = id,
      tree = tree_name,
      nearest_group = ordered_group[1],
      nearest_id = ordered[1],
      nearest_C2_rank = c2_rank,
      nearest_C2_id = c2_id,
      n_C2_in_10NN = sum(ordered_group[1:10] == "C2", na.rm = TRUE),
      n_C2_in_20NN = sum(ordered_group[1:20] == "C2", na.rm = TRUE),
      n_C3_in_10NN = sum(ordered_group[1:10] == "C3", na.rm = TRUE),
      first10_groups = paste(ordered_group[1:10], collapse = ","),
      first10_ids = paste(ordered[1:10], collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

tree_c3_clades <- function(tree_file, tree_name) {
  tree <- read.tree(tree_file)
  helper <- make_tree_helpers(tree)
  internal_nodes <- (helper$n_tip + 1):(helper$n_tip + tree$Nnode)

  rows <- list()
  for (node in internal_nodes) {
    tips <- helper$tips_under(node)
    gg <- groups$plot_group[match(tips, groups$IID)]
    if (length(tips) == 0 || !all(gg == "C3", na.rm = TRUE)) next

    parent_node <- helper$parent[as.character(node)]
    parent_tips <- character(0)
    sister_tips <- character(0)
    is_top_pure_c3 <- TRUE

    if (!is.na(parent_node)) {
      parent_tips <- helper$tips_under(as.integer(parent_node))
      parent_groups <- groups$plot_group[match(parent_tips, groups$IID)]
      is_top_pure_c3 <- !all(parent_groups == "C3", na.rm = TRUE)

      sib_nodes <- setdiff(helper$children[[as.character(parent_node)]], node)
      if (length(sib_nodes) > 0) {
        sister_tips <- unlist(lapply(sib_nodes, helper$tips_under), use.names = FALSE)
      }
    }

    if (!is_top_pure_c3) next
    sister_groups <- groups$plot_group[match(sister_tips, groups$IID)]
    parent_groups <- groups$plot_group[match(parent_tips, groups$IID)]

    rows[[length(rows) + 1]] <- data.frame(
      tree = tree_name,
      node = node,
      n_C3 = length(tips),
      tips = paste(tips, collapse = ","),
      sister_counts = group_counts_text(sister_tips),
      parent_counts = group_counts_text(parent_tips),
      sister_has_C2 = any(sister_groups == "C2", na.rm = TRUE),
      parent_has_C2 = any(parent_groups == "C2", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

ml_metrics <- tree_tip_metrics(ml_file, "ML")
nj_metrics <- tree_tip_metrics(nj_file, "NJ")
metrics <- rbind(ml_metrics, nj_metrics)
write.table(metrics, file.path(out_dir, "C3_tip_neighbor_metrics.long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ml_clades <- tree_c3_clades(ml_file, "ML")
nj_clades <- tree_c3_clades(nj_file, "NJ")
clades <- rbind(ml_clades, nj_clades)
write.table(clades, file.path(out_dir, "C3_top_pure_clades_with_sister_context.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

wide <- reshape(metrics, idvar = "IID", timevar = "tree", direction = "wide")

candidate_from_clades <- character(0)
if (nrow(clades) > 0) {
  candidate_clades <- clades[
    clades$n_C3 <= 12 &
      (clades$sister_has_C2 | clades$parent_has_C2) &
      grepl("C2=[1-9]", clades$sister_counts),
  ]
  if (nrow(candidate_clades) > 0) {
    candidate_from_clades <- unique(unlist(strsplit(candidate_clades$tips, ",")))
  }
}

wide$ML_C2_near <- wide$nearest_C2_rank.ML <= 15 | wide$n_C2_in_20NN.ML >= 3
wide$NJ_C2_near <- wide$nearest_C2_rank.NJ <= 15 | wide$n_C2_in_20NN.NJ >= 3
wide$in_C2_adjacent_C3_clade <- wide$IID %in% candidate_from_clades
wide$candidate_level <- ifelse(
  wide$ML_C2_near & wide$NJ_C2_near,
  "both_trees_C2_near",
  ifelse(wide$in_C2_adjacent_C3_clade, "C3_clade_adjacent_to_C2", "not_C2_near")
)
candidate_ids <- unique(wide$IID[wide$candidate_level != "not_C2_near"])

centroids <- aggregate(cbind(PC1, PC2, PC3) ~ plot_group, data = pca[pca$plot_group %in% c("C1", "C2", "C3", "C4"), ], mean)
get_centroid <- function(g) as.numeric(centroids[centroids$plot_group == g, c("PC1", "PC2", "PC3")])
c2_centroid <- get_centroid("C2")
c3_centroid <- get_centroid("C3")

summary_rows <- merge(wide, pca[, c("IID", "plot_group", "PC1", "PC2", "PC3")], by = "IID", all.x = TRUE)
summary_rows$dist_to_C2_centroid_PC123 <- apply(summary_rows[, c("PC1", "PC2", "PC3")], 1, function(x) sqrt(sum((x - c2_centroid)^2)))
summary_rows$dist_to_C3_centroid_PC123 <- apply(summary_rows[, c("PC1", "PC2", "PC3")], 1, function(x) sqrt(sum((x - c3_centroid)^2)))
summary_rows$PCA_centroid_closer <- ifelse(summary_rows$dist_to_C2_centroid_PC123 < summary_rows$dist_to_C3_centroid_PC123, "C2", "C3")

fam_file <- file.path(main_dir, "03_admixture", "sativa182.admix.fam")
q_file <- file.path(main_dir, "03_admixture", "sativa182.admix.4.Q")
if (file.exists(fam_file) && file.exists(q_file)) {
  fam <- read.table(fam_file, stringsAsFactors = FALSE)
  q <- read.table(q_file, stringsAsFactors = FALSE)
  colnames(q) <- paste0("Q", seq_len(ncol(q)))
  admix <- cbind(IID = as.character(fam[[2]]), q)
  admix <- merge(admix, groups[, c("IID", "plot_group")], by = "IID")
  q_cols <- paste0("Q", seq_len(ncol(q)))
  means <- aggregate(admix[admix$plot_group %in% c("C1", "C2", "C3", "C4"), q_cols],
                     by = list(plot_group = admix$plot_group[admix$plot_group %in% c("C1", "C2", "C3", "C4")]), mean)
  component_to_group <- data.frame(component = q_cols, mapped_group = NA_character_, stringsAsFactors = FALSE)
  for (qc in q_cols) {
    component_to_group$mapped_group[component_to_group$component == qc] <- means$plot_group[which.max(means[[qc]])]
  }
  q_c2 <- component_to_group$component[component_to_group$mapped_group == "C2"]
  q_c3 <- component_to_group$component[component_to_group$mapped_group == "C3"]
  admix$ADMIX_C2_component <- if (length(q_c2)) rowSums(admix[, q_c2, drop = FALSE]) else NA_real_
  admix$ADMIX_C3_component <- if (length(q_c3)) rowSums(admix[, q_c3, drop = FALSE]) else NA_real_
  admix$ADMIX_max_component <- q_cols[max.col(admix[, q_cols], ties.method = "first")]
  admix$ADMIX_max_group <- component_to_group$mapped_group[match(admix$ADMIX_max_component, component_to_group$component)]
  summary_rows <- merge(summary_rows, admix[, c("IID", "ADMIX_C2_component", "ADMIX_C3_component", "ADMIX_max_group")], by = "IID", all.x = TRUE)
}

summary_rows <- summary_rows[order(
  summary_rows$candidate_level,
  summary_rows$nearest_C2_rank.ML,
  summary_rows$nearest_C2_rank.NJ
), ]
write.table(summary_rows, file.path(out_dir, "C3_runaway_candidate_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

confirmed <- summary_rows[summary_rows$IID %in% candidate_ids, ]
write.table(confirmed, file.path(out_dir, "C3_runaway_confirmed_samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_pca_highlight <- function(xpc, ypc, suffix) {
  png(file.path(out_dir, paste0("PCA_", suffix, "_highlight_C3_runaway.png")), width = 1650, height = 1350, res = 220, type = "cairo")
  par(mar = c(4.2, 4.2, 1.2, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in names(palette)) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g], bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  if (length(candidate_ids) > 0) {
    idx <- match(candidate_ids, pca$IID)
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = 21, col = "#f0b000", bg = NA, cex = 1.7, lwd = 1.7)
    text(pca[[xpc]][idx], pca[[ypc]][idx], labels = candidate_ids, pos = 3, cex = 0.62, col = "#333333")
  }
  legend("topright", legend = c(names(palette), "C3 candidate"),
         pch = c(shape, 21), col = c(palette, "#f0b000"),
         pt.bg = c(palette, NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, 5), 1.35))
  dev.off()

  pdf(file.path(out_dir, paste0("PCA_", suffix, "_highlight_C3_runaway.pdf")), width = 6.2, height = 5.1, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.2, 1.0))
  plot(pca[[xpc]], pca[[ypc]], type = "n", xlab = xpc, ylab = ypc, bty = "l", cex.lab = 1.2)
  for (g in names(palette)) {
    idx <- pca$plot_group == g
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = shape[g], col = palette[g], bg = palette[g], cex = ifelse(g == "out", 0.85, 0.82))
  }
  if (length(candidate_ids) > 0) {
    idx <- match(candidate_ids, pca$IID)
    points(pca[[xpc]][idx], pca[[ypc]][idx], pch = 21, col = "#f0b000", bg = NA, cex = 1.7, lwd = 1.7)
    text(pca[[xpc]][idx], pca[[ypc]][idx], labels = candidate_ids, pos = 3, cex = 0.62, col = "#333333")
  }
  legend("topright", legend = c(names(palette), "C3 candidate"),
         pch = c(shape, 21), col = c(palette, "#f0b000"),
         pt.bg = c(palette, NA), bty = "n", cex = 0.82, pt.cex = c(rep(1.0, 5), 1.35))
  dev.off()
}

plot_tree_highlight <- function(tree_file, tree_title, prefix) {
  tree <- read.tree(tree_file)
  tip_group <- groups$plot_group[match(tree$tip.label, groups$IID)]
  names(tip_group) <- tree$tip.label
  tip_col <- palette[tip_group]
  tip_col[is.na(tip_col)] <- "#777777"
  tip_pch <- ifelse(tree$tip.label %in% candidate_ids, 21, 16)
  tip_bg <- ifelse(tree$tip.label %in% candidate_ids, "#f0b000", tip_col)
  tip_cex <- ifelse(tree$tip.label %in% candidate_ids, 0.95, 0.50)

  png(file.path(out_dir, paste0(prefix, ".highlight_C3_runaway.png")), width = 1900, height = 1900, res = 300, type = "cairo")
  par(mar = c(0, 0, 1, 0))
  plot.phylo(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE, edge.width = 0.55, edge.color = "#303030")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index], pch = tip_pch, col = ifelse(tree$tip.label %in% candidate_ids, "#f0b000", tip_col),
         bg = tip_bg, cex = tip_cex, lwd = ifelse(tree$tip.label %in% candidate_ids, 1.2, 0.6))
  cand_idx <- which(tree$tip.label %in% candidate_ids)
  text(pp$xx[cand_idx], pp$yy[cand_idx], labels = tree$tip.label[cand_idx], cex = 0.48, pos = 3, col = "#222222")
  mtext(tree_title, side = 3, line = -1.0, cex = 1.15, font = 2)
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, ".highlight_C3_runaway.pdf")), width = 6.3, height = 6.3, useDingbats = FALSE)
  par(mar = c(0, 0, 1, 0))
  plot.phylo(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE, edge.width = 0.55, edge.color = "#303030")
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_index <- seq_len(Ntip(tree))
  points(pp$xx[tip_index], pp$yy[tip_index], pch = tip_pch, col = ifelse(tree$tip.label %in% candidate_ids, "#f0b000", tip_col),
         bg = tip_bg, cex = tip_cex, lwd = ifelse(tree$tip.label %in% candidate_ids, 1.2, 0.6))
  cand_idx <- which(tree$tip.label %in% candidate_ids)
  text(pp$xx[cand_idx], pp$yy[cand_idx], labels = tree$tip.label[cand_idx], cex = 0.48, pos = 3, col = "#222222")
  mtext(tree_title, side = 3, line = -1.0, cex = 1.15, font = 2)
  dev.off()
}

plot_pca_highlight("PC1", "PC2", "PC1_PC2")
plot_pca_highlight("PC1", "PC3", "PC1_PC3")
plot_tree_highlight(ml_file, "ML tree with C3 candidates", "ML_tree")
plot_tree_highlight(nj_file, "NJ tree with C3 candidates", "NJ_tree")

cat("Confirmed candidate IDs:\n")
cat(paste(candidate_ids, collapse = "\n"), "\n")
cat("Files:\n")
cat(file.path(out_dir, "C3_runaway_confirmed_samples.tsv"), "\n")
cat(file.path(out_dir, "C3_runaway_candidate_summary.tsv"), "\n")
cat(file.path(out_dir, "PCA_PC1_PC2_highlight_C3_runaway.png"), "\n")
