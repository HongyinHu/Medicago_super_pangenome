args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript extract_admixture_K4_special.R <main_dir>")
}

main_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- file.path(main_dir, "09_C3_C2_special_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fam <- read.table(file.path(main_dir, "03_admixture", "sativa182.admix.fam"), stringsAsFactors = FALSE)
q <- read.table(file.path(main_dir, "03_admixture", "sativa182.admix.4.Q"), stringsAsFactors = FALSE)
colnames(q) <- paste0("Q", seq_len(ncol(q)))
admix <- cbind(IID = as.character(fam[[2]]), q)

groups <- read.delim(file.path(main_dir, "06_summary", "pca_group_assignment_K4_with_outgroup.tsv"), stringsAsFactors = FALSE)
groups$IID <- as.character(groups$IID)
admix <- merge(admix, groups[, c("IID", "plot_group")], by = "IID")

ingroup <- admix[admix$plot_group %in% c("C1", "C2", "C3", "C4"), ]
q_cols <- paste0("Q", seq_len(ncol(q)))
means <- aggregate(ingroup[, q_cols], by = list(plot_group = ingroup$plot_group), mean)

component_to_group <- data.frame(
  component = q_cols,
  mapped_group = NA_character_,
  mapped_group_mean = NA_real_,
  stringsAsFactors = FALSE
)
for (qc in q_cols) {
  ii <- which.max(means[[qc]])
  component_to_group[component_to_group$component == qc, "mapped_group"] <- means$plot_group[ii]
  component_to_group[component_to_group$component == qc, "mapped_group_mean"] <- means[[qc]][ii]
}

admix$C2_component_sum <- rowSums(admix[, component_to_group$component[component_to_group$mapped_group == "C2"], drop = FALSE])
admix$C3_component_sum <- rowSums(admix[, component_to_group$component[component_to_group$mapped_group == "C3"], drop = FALSE])
admix$max_component <- q_cols[max.col(admix[, q_cols], ties.method = "first")]
admix$max_component_group <- component_to_group$mapped_group[match(admix$max_component, component_to_group$component)]

special_file <- file.path(out_dir, "C3_tree_near_C2_special_samples.tsv")
if (file.exists(special_file)) {
  special <- read.delim(special_file, stringsAsFactors = FALSE, check.names = FALSE)
  special_q <- merge(special[, "IID", drop = FALSE], admix, by = "IID", all.x = TRUE)
} else {
  special_q <- admix[admix$IID == "474", ]
}

write.table(means, file.path(out_dir, "admixture_K4_component_means_by_PCA_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(component_to_group, file.path(out_dir, "admixture_K4_component_to_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(special_q, file.path(out_dir, "C3_tree_near_C2_special_samples.ADMIXTURE_K4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

print(component_to_group)
print(special_q)
