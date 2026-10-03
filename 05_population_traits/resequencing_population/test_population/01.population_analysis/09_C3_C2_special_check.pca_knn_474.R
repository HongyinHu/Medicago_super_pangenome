args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript pca_knn_474.R <main_dir>")
main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- file.path(main_dir, "09_C3_C2_special_check")
pca <- read.delim(file.path(main_dir, "06_summary", "pca_group_assignment_K4_with_outgroup.tsv"), stringsAsFactors = FALSE)
pca$IID <- as.character(pca$IID)
for (cc in c("PC1", "PC2", "PC3")) pca[[cc]] <- as.numeric(pca[[cc]])

knn <- function(cols, label) {
  mat <- as.matrix(pca[, cols])
  rownames(mat) <- pca$IID
  d <- sqrt(rowSums((t(t(mat) - mat["474", ]))^2))
  d <- d[names(d) != "474"]
  ids <- names(sort(d))[1:10]
  data.frame(
    space = label,
    rank = seq_along(ids),
    IID = ids,
    plot_group = pca$plot_group[match(ids, pca$IID)],
    distance = unname(d[ids]),
    stringsAsFactors = FALSE
  )
}
res <- rbind(knn(c("PC1", "PC2"), "PC1_PC2"), knn(c("PC1", "PC2", "PC3"), "PC1_PC2_PC3"))
write.table(res, file.path(out_dir, "sample_474_PCA_nearest_neighbors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
