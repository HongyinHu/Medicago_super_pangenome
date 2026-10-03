#!/usr/bin/env bash
set -euo pipefail

MAIN="path/to/project/38.medicago_resequence/test_population/01.population_analysis"
OUT="$MAIN/11_revised_C3_to_C2"
RSCRIPT="path/to/home/anaconda3/envs/biosofeware/bin/Rscript"
VCFTOOLS="path/to/home/anaconda3/envs/biosofeware/bin/vcftools"
VCF="$MAIN/07_pi_fst_K4/01_vcf/sativa182.filtered.vcf"

mkdir -p "$OUT/00_groups" "$OUT/01_pca" "$OUT/02_trees" "$OUT/03_pi_fst/02_pi" "$OUT/03_pi_fst/03_fst" "$OUT/04_plots" "$OUT/99_scripts"

cat > "$OUT/00_groups/reassigned_C3_to_C2.ids" <<'IDS'
453
M9_1
492
M3_1
493
M20_1
M11_1
Ms489
474
IDS

cat > "$OUT/99_scripts/revised_groups_pca_trees.R" <<'RSCRIPT'
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
RSCRIPT

"$RSCRIPT" "$OUT/99_scripts/revised_groups_pca_trees.R" "$MAIN" "$OUT"

for g in C1 C2 C3 C4; do
  "$VCFTOOLS" --vcf "$VCF" --keep "$OUT/00_groups/$g.keep" --site-pi --out "$OUT/03_pi_fst/02_pi/$g" \
    > "$OUT/03_pi_fst/02_pi/$g.stdout.log" 2>&1
done

for pair in C1:C2 C1:C3 C1:C4 C2:C3 C2:C4 C3:C4; do
  g1="${pair%:*}"
  g2="${pair#*:}"
  "$VCFTOOLS" --vcf "$VCF" --weir-fst-pop "$OUT/00_groups/$g1.keep" --weir-fst-pop "$OUT/00_groups/$g2.keep" \
    --out "$OUT/03_pi_fst/03_fst/${g1}_vs_${g2}" \
    > "$OUT/03_pi_fst/03_fst/${g1}_vs_${g2}.stdout.log" 2>&1
done

cat > "$OUT/99_scripts/summarize_plot_pi_fst.R" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript summarize_plot_pi_fst.R <out_dir>")
}

out_dir <- normalizePath(args[1], mustWork = TRUE)
groups <- c("C1", "C2", "C3", "C4")
counts <- read.delim(file.path(out_dir, "00_groups", "revised_group_counts.tsv"), stringsAsFactors = FALSE)
counts <- counts[counts$group %in% groups, ]

pi_rows <- lapply(groups, function(g) {
  f <- file.path(out_dir, "03_pi_fst", "02_pi", paste0(g, ".sites.pi"))
  x <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  vals <- suppressWarnings(as.numeric(x$PI))
  data.frame(
    group = g,
    n_samples = counts$n[match(g, counts$group)],
    n_sites = length(vals),
    n_valid_sites = sum(is.finite(vals)),
    n_nan_sites = sum(!is.finite(vals)),
    pi_snp_mean_valid_sites = mean(vals[is.finite(vals)]),
    stringsAsFactors = FALSE
  )
})
pi_summary <- do.call(rbind, pi_rows)
write.table(pi_summary, file.path(out_dir, "03_pi_fst", "pi_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parse_fst_log <- function(log_file) {
  lines <- readLines(log_file, warn = FALSE)
  mean_line <- grep("Weir and Cockerham mean Fst estimate", lines, value = TRUE)
  weighted_line <- grep("Weir and Cockerham weighted Fst estimate", lines, value = TRUE)
  mean_fst <- as.numeric(sub(".*: *", "", mean_line[length(mean_line)]))
  weighted_fst <- as.numeric(sub(".*: *", "", weighted_line[length(weighted_line)]))
  c(mean_fst = mean_fst, weighted_fst = weighted_fst)
}

pairs <- list(c("C1", "C2"), c("C1", "C3"), c("C1", "C4"), c("C2", "C3"), c("C2", "C4"), c("C3", "C4"))
fst_rows <- lapply(pairs, function(p) {
  log_file <- file.path(out_dir, "03_pi_fst", "03_fst", paste0(p[1], "_vs_", p[2], ".stdout.log"))
  vals <- parse_fst_log(log_file)
  data.frame(group1 = p[1], group2 = p[2], mean_fst = vals["mean_fst"], weighted_fst = vals["weighted_fst"], stringsAsFactors = FALSE)
})
fst_summary <- do.call(rbind, fst_rows)
write.table(fst_summary, file.path(out_dir, "03_pi_fst", "fst_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_dir <- file.path(out_dir, "04_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

coords <- data.frame(
  group = groups,
  x = c(-1, 0, 1, 0),
  y = c(0, 1, 0, -1),
  color = c("#ef6a6a", "#405a93", "#159c89", "#b786bd"),
  stringsAsFactors = FALSE
)
coords$pi <- pi_summary$pi_snp_mean_valid_sites[match(coords$group, pi_summary$group)]
pi_range <- range(coords$pi, finite = TRUE)
coords$radius <- if (diff(pi_range) == 0) 0.31 else 0.25 + 0.08 * (coords$pi - pi_range[1]) / diff(pi_range)

edge_key <- function(a, b) paste(sort(c(a, b)), collapse = "-")
fst_summary$key <- mapply(edge_key, fst_summary$group1, fst_summary$group2)
fst_lookup <- setNames(fst_summary$weighted_fst, fst_summary$key)
edges <- data.frame(
  group1 = c("C1", "C2", "C1", "C1", "C2", "C3"),
  group2 = c("C2", "C3", "C3", "C4", "C4", "C4"),
  label_x = c(-0.72, 0.72, -0.18, -0.72, 0.40, 0.72),
  label_y = c(0.68, 0.68, 0.23, -0.68, -0.02, -0.68),
  stringsAsFactors = FALSE
)
edges$key <- mapply(edge_key, edges$group1, edges$group2)
edges$fst <- fst_lookup[edges$key]

draw_label <- function(x, y, label, cex = 1.35) {
  pad_x <- 0.045
  pad_y <- 0.035
  w <- strwidth(label, cex = cex)
  h <- strheight(label, cex = cex)
  rect(x - w / 2 - pad_x, y - h / 2 - pad_y, x + w / 2 + pad_x, y + h / 2 + pad_y, col = "white", border = NA)
  text(x, y, label, cex = cex, col = "#111111")
}

draw_panel <- function() {
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = c(-1.45, 1.45), ylim = c(-1.35, 1.35), asp = 1, axes = FALSE, xlab = "", ylab = "", bty = "n")
  for (i in seq_len(nrow(edges))) {
    p1 <- coords[coords$group == edges$group1[i], ]
    p2 <- coords[coords$group == edges$group2[i], ]
    segments(p1$x, p1$y, p2$x, p2$y, lty = 2, lwd = 1.4, col = "#222222")
  }
  for (i in seq_len(nrow(edges))) {
    draw_label(edges$label_x[i], edges$label_y[i], sprintf("%.3f", edges$fst[i]))
  }
  for (i in seq_len(nrow(coords))) {
    symbols(coords$x[i], coords$y[i], circles = coords$radius[i], inches = FALSE,
            add = TRUE, bg = coords$color[i], fg = coords$color[i], lwd = 1.5)
    text(coords$x[i], coords$y[i] + 0.055, coords$group[i], col = "white", cex = 1.25, font = 2)
    text(coords$x[i], coords$y[i] - 0.085, sprintf("%.3f", coords$pi[i]), col = "white", cex = 1.12, font = 2)
  }
}

png(file.path(plot_dir, "pi_fst_revised_C3_to_C2_network.png"), width = 1200, height = 1000, res = 220, type = "cairo")
draw_panel()
dev.off()
pdf(file.path(plot_dir, "pi_fst_revised_C3_to_C2_network.pdf"), width = 5.5, height = 4.8, useDingbats = FALSE)
draw_panel()
dev.off()

print(pi_summary)
print(fst_summary)
RSCRIPT

"$RSCRIPT" "$OUT/99_scripts/summarize_plot_pi_fst.R" "$OUT"

find "$OUT" -maxdepth 4 -type f | sort
