suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: plot_mixed_illumina_expression_sv_distance.R WORK MANIFEST TX2GENE SV_CONTEXT PAV_MATRIX")
}

work <- args[[1]]
manifest_file <- args[[2]]
tx2gene_file <- args[[3]]
sv_context_file <- args[[4]]
pav_file <- args[[5]]

ref_group <- "genome_Msa"
tables_dir <- file.path(work, "tables")
figures_dir <- file.path(work, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
tx2gene <- read.delim(tx2gene_file, stringsAsFactors = FALSE, check.names = FALSE)

message("Reading Salmon quant.sf files")
gene_tables <- list()
for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[i]
  qfile <- file.path(work, "quants", sample_id, "quant.sf")
  if (!file.exists(qfile)) {
    stop("Missing quant.sf: ", qfile)
  }
  q <- read.delim(qfile, stringsAsFactors = FALSE, check.names = FALSE)
  q <- q[, c("Name", "TPM")]
  q <- merge(q, tx2gene, by.x = "Name", by.y = "tx_id", all.x = FALSE, all.y = FALSE)
  g <- aggregate(TPM ~ gene_id, data = q, FUN = sum)
  colnames(g)[2] <- sample_id
  gene_tables[[sample_id]] <- g
}

gene_tpm <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), gene_tables)
gene_tpm[is.na(gene_tpm)] <- 0
write.table(gene_tpm, file.path(tables_dir, "gene_expression_TPM_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

groups <- unique(manifest$group)
if (!(ref_group %in% groups)) {
  stop("Reference group not found in manifest: ", ref_group)
}
target_groups <- setdiff(groups, ref_group)
if (length(target_groups) == 0) {
  stop("No target groups available")
}

gene_summary <- gene_tpm
group_mean_cols <- c()
for (grp in groups) {
  samples <- manifest$sample_id[manifest$group == grp]
  colname <- paste0("mean_TPM_", grp)
  gene_summary[[colname]] <- rowMeans(gene_summary[, samples, drop = FALSE])
  group_mean_cols <- c(group_mean_cols, colname)
}

write.table(gene_summary[, c("gene_id", group_mean_cols)],
            file.path(tables_dir, "gene_expression_group_mean_TPM.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Reading SV context and PAV matrix")
sv_context <- read.delim(sv_context_file, stringsAsFactors = FALSE, check.names = FALSE)
pav <- read.delim(pav_file, stringsAsFactors = FALSE, check.names = FALSE)
sv_context$distance_to_gene_bp <- suppressWarnings(as.numeric(sv_context$distance_to_gene_bp))

needed_pav_cols <- c("SV_ID", target_groups)
missing_pav_cols <- setdiff(needed_pav_cols, colnames(pav))
if (length(missing_pav_cols) > 0) {
  stop("PAV matrix missing columns: ", paste(missing_pav_cols, collapse = ", "))
}

near_sv <- sv_context %>%
  filter(!is.na(nearest_gene_id),
         nearest_gene_id != "",
         !is.na(distance_to_gene_bp),
         distance_to_gene_bp <= 2000) %>%
  select(SV_ID, nearest_gene_id, distance_to_gene_bp)

near_pav <- near_sv %>%
  inner_join(pav[, needed_pav_cols], by = "SV_ID")

present_as_logical <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(!is.na(x) & x > 0)
  }
  xx <- tolower(as.character(x))
  !is.na(xx) & xx %in% c("1", "true", "yes", "present", "presence")
}

ref_mean_col <- paste0("mean_TPM_", ref_group)
all_compare <- list()
gene_sv_status_tables <- list()

for (grp in target_groups) {
  target_mean_col <- paste0("mean_TPM_", grp)
  with_sv_genes <- unique(near_pav$nearest_gene_id[present_as_logical(near_pav[[grp]])])
  tmp <- gene_summary %>%
    transmute(
      gene_id = gene_id,
      comparison = paste(ref_group, "vs", grp),
      target_group = grp,
      mean_TPM_ref = .data[[ref_mean_col]],
      mean_TPM_target = .data[[target_mean_col]],
      log2FC_target_vs_ref = log2(.data[[target_mean_col]] + 1) - log2(.data[[ref_mean_col]] + 1),
      expression_divergence = abs(log2FC_target_vs_ref),
      expressed_any = .data[[ref_mean_col]] > 1 | .data[[target_mean_col]] > 1,
      SV_within_2kb = gene_id %in% with_sv_genes,
      sv_status = ifelse(SV_within_2kb, "With SV", "Without SV")
    )
  all_compare[[grp]] <- tmp
  gene_sv_status_tables[[grp]] <- tmp[, c("gene_id", "comparison", "target_group", "SV_within_2kb", "sv_status")]
}

all_compare_df <- bind_rows(all_compare)
write.table(all_compare_df, file.path(tables_dir, "genome_Msa_vs_targets_expression_divergence.all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(gene_sv_status_tables), file.path(tables_dir, "gene_species_SV2kb_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_df <- all_compare_df %>%
  filter(expressed_any) %>%
  mutate(
    sv_status = factor(sv_status, levels = c("Without SV", "With SV")),
    comparison = factor(comparison, levels = paste(ref_group, "vs", target_groups))
  )
write.table(plot_df, file.path(tables_dir, "genome_Msa_vs_targets_expression_divergence.expressed_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(20260704)
boot_n <- as.integer(Sys.getenv("BOOT_N", "10000"))
bootstrap_ci <- function(x, n = boot_n) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(boot_mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  b <- replicate(n, mean(sample(x, size = length(x), replace = TRUE)))
  c(boot_mean = mean(x), ci_low = unname(quantile(b, 0.025)), ci_high = unname(quantile(b, 0.975)))
}

summary_tbl <- plot_df %>%
  group_by(comparison, target_group, sv_status) %>%
  summarise(
    n_genes = n(),
    mean_divergence = mean(expression_divergence),
    median_divergence = median(expression_divergence),
    q25 = quantile(expression_divergence, 0.25),
    q75 = quantile(expression_divergence, 0.75),
    .groups = "drop"
  )

boot_tbl <- do.call(rbind, lapply(split(plot_df, list(plot_df$comparison, plot_df$sv_status), drop = TRUE), function(d) {
  ci <- bootstrap_ci(d$expression_divergence)
  data.frame(comparison = as.character(d$comparison[1]),
             sv_status = as.character(d$sv_status[1]),
             t(ci),
             row.names = NULL)
}))

p_tbl <- do.call(rbind, lapply(split(plot_df, plot_df$comparison), function(d) {
  pval <- tryCatch(
    wilcox.test(expression_divergence ~ sv_status, data = d, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  data.frame(comparison = as.character(d$comparison[1]), wilcox_p = pval, row.names = NULL)
}))

summary_tbl <- summary_tbl %>%
  left_join(boot_tbl, by = c("comparison", "sv_status")) %>%
  left_join(p_tbl, by = "comparison") %>%
  mutate(
    bootstrap_n = boot_n,
    n_ref_RNA = sum(manifest$group == ref_group),
    n_target_RNA = vapply(target_group, function(g) sum(manifest$group == g), integer(1))
  )

write.table(summary_tbl, file.path(tables_dir, "genome_Msa_vs_targets_expression_divergence.summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

label_tbl <- p_tbl %>%
  mutate(
    comparison = factor(comparison, levels = paste(ref_group, "vs", target_groups)),
    label = ifelse(is.na(wilcox_p), "P = NA",
                   ifelse(wilcox_p < 2.2e-16, "P < 2.2e-16",
                          paste0("P = ", format(wilcox_p, digits = 2, scientific = TRUE))))
  )
ymax_by_comp <- plot_df %>%
  group_by(comparison) %>%
  summarise(y = quantile(expression_divergence, 0.995, na.rm = TRUE), .groups = "drop")
label_tbl <- left_join(label_tbl, ymax_by_comp, by = "comparison")

mean_df <- summary_tbl %>%
  mutate(sv_status = factor(sv_status, levels = c("Without SV", "With SV")),
         comparison = factor(comparison, levels = paste(ref_group, "vs", target_groups)))

cols <- c("Without SV" = "#1F78B4", "With SV" = "#D8A700")

p <- ggplot(plot_df, aes(x = comparison, y = expression_divergence, color = sv_status)) +
  geom_boxplot(aes(group = interaction(comparison, sv_status)),
               position = position_dodge(width = 0.72),
               width = 0.52, fill = "white", linewidth = 0.5,
               outlier.size = 0.15, outlier.alpha = 0.15) +
  geom_crossbar(data = mean_df,
                aes(x = comparison, y = mean_divergence, ymin = mean_divergence,
                    ymax = mean_divergence, group = sv_status),
                position = position_dodge(width = 0.72),
                inherit.aes = FALSE, width = 0.48, color = "#D62728", linewidth = 0.35) +
  geom_text(data = label_tbl, aes(x = comparison, y = y, label = label),
            inherit.aes = FALSE, angle = 90, size = 2.4, hjust = 0, vjust = 0.5) +
  scale_color_manual(values = cols, name = NULL) +
  labs(
    x = NULL,
    y = expression("Expression divergence  |log"[2]*"(mean TPM + 1)|"),
    title = "Mixed-tissue Illumina RNA expression divergence",
    subtitle = "Each comparison uses species-specific panSVs within 2 kb of genes"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8, hjust = 0),
    axis.text.x = element_text(size = 7.2, angle = 38, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 18, 8, 10)
  )

pdf_file <- file.path(figures_dir, "genome_Msa_vs_8species_SV2kb_expression_divergence.pdf")
png_file <- file.path(figures_dir, "genome_Msa_vs_8species_SV2kb_expression_divergence.png")
tiff_file <- file.path(figures_dir, "genome_Msa_vs_8species_SV2kb_expression_divergence.tiff")

cairo_pdf(pdf_file, width = 7.2, height = 4.5)
print(p)
dev.off()

ragg::agg_png(png_file, width = 7.2, height = 4.5, units = "in", res = 600, background = "white")
print(p)
dev.off()

ragg::agg_tiff(tiff_file, width = 7.2, height = 4.5, units = "in", res = 600,
               compression = "lzw", background = "white")
print(p)
dev.off()

message("Wrote figure: ", pdf_file)

