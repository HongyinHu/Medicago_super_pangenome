args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript fix_revised_pca_legend.R <out_dir>")
}

out_dir <- normalizePath(args[1], mustWork = TRUE)
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
  idx_re <- pca$reassigned_from_C3_to_C2
  points(pca[[xpc]][idx_re], pca[[ypc]][idx_re], pch = 21, col = "#f0b000", bg = NA, cex = 1.65, lwd = 1.7)
  legend("topleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(shape[counts$group], 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.82,
         pt.cex = c(rep(1.0, nrow(counts)), 1.35))
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
  legend("topleft", legend = c(legend_labels, "C3->C2 reassigned"),
         pch = c(shape[counts$group], 21), col = c(palette[counts$group], "#f0b000"),
         pt.bg = c(palette[counts$group], NA), bty = "n", cex = 0.82,
         pt.cex = c(rep(1.0, nrow(counts)), 1.35))
  dev.off()
}

draw_pca("PC1", "PC2", file.path(out_dir, "01_pca", "sativa182_revised_C3_to_C2.PCA_PC1_PC2"))
draw_pca("PC1", "PC3", file.path(out_dir, "01_pca", "sativa182_revised_C3_to_C2.PCA_PC1_PC3"))
