suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: plot_pav_gene_structure_expression_divergence.R WORK SV_CONTEXT PAV_MATRIX")
}

work <- args[[1]]
sv_context_file <- args[[2]]
pav_file <- args[[3]]

tables_dir <- file.path(work, "tables")
figures_dir <- file.path(work, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(tables_dir, "genome_Msa_vs_targets_expression_divergence.expressed_genes.tsv")
if (!file.exists(expr_file)) {
  stop("Missing expression divergence table: ", expr_file)
}

message("Reading expression divergence table")
expr <- read.delim(expr_file, stringsAsFactors = FALSE, check.names = FALSE)
target_groups <- unique(expr$target_group)

message("Reading SV context")
ctx <- read.delim(sv_context_file, stringsAsFactors = FALSE, check.names = FALSE)
ctx$distance_to_gene_bp <- suppressWarnings(as.numeric(ctx$distance_to_gene_bp))

ctx_small <- ctx %>%
  filter(!is.na(nearest_gene_id), nearest_gene_id != "") %>%
  mutate(
    feature_category = case_when(
      category4 == "Gene body" & gene_body_detail == "exon_overlap" ~ "exon",
      category4 == "Gene body" & gene_body_detail == "intron_only" ~ "intron",
      category4 == "Upstream" & !is.na(distance_to_gene_bp) & distance_to_gene_bp <= 2000 ~ "2-kb upstream",
      category4 == "Downstream" & !is.na(distance_to_gene_bp) & distance_to_gene_bp <= 2000 ~ "2-kb downstream",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(feature_category)) %>%
  transmute(
    SV_ID = SV_ID,
    gene_id = nearest_gene_id,
    feature_category = feature_category,
    distance_to_gene_bp = distance_to_gene_bp
  ) %>%
  distinct()

message("Reading PAV matrix")
pav <- read.delim(pav_file, stringsAsFactors = FALSE, check.names = FALSE)
needed_cols <- c("SV_ID", target_groups)
missing_cols <- setdiff(needed_cols, colnames(pav))
if (length(missing_cols) > 0) {
  stop("PAV matrix missing columns: ", paste(missing_cols, collapse = ", "))
}

present_as_logical <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(!is.na(x) & x > 0)
  }
  xx <- tolower(as.character(x))
  !is.na(xx) & xx %in% c("1", "true", "yes", "present", "presence")
}

pav_long <- pav[, needed_cols] %>%
  pivot_longer(cols = all_of(target_groups), names_to = "target_group", values_to = "present") %>%
  filter(present_as_logical(present)) %>%
  select(SV_ID, target_group)

membership_multi <- pav_long %>%
  inner_join(ctx_small, by = "SV_ID") %>%
  distinct(target_group, gene_id, feature_category, SV_ID, distance_to_gene_bp)

expr_key <- expr %>%
  select(gene_id, comparison, target_group, expression_divergence, mean_TPM_ref, mean_TPM_target,
         log2FC_target_vs_ref)

membership_expr_multi <- membership_multi %>%
  inner_join(expr_key, by = c("target_group", "gene_id")) %>%
  mutate(feature_category = factor(feature_category,
                                   levels = c("exon", "intron", "2-kb downstream", "2-kb upstream")))

write.table(membership_expr_multi,
            file.path(tables_dir, "PAV_gene_structure_expression_divergence.multi_category.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

priority_tbl <- tibble(
  feature_category = c("exon", "intron", "2-kb upstream", "2-kb downstream"),
  feature_priority = c(1L, 2L, 3L, 4L)
)

membership_priority <- membership_multi %>%
  inner_join(priority_tbl, by = "feature_category") %>%
  group_by(target_group, gene_id) %>%
  arrange(feature_priority, distance_to_gene_bp, .by_group = TRUE) %>%
  summarise(
    feature_category = first(feature_category),
    feature_priority = first(feature_priority),
    n_feature_categories = n_distinct(feature_category),
    n_nearby_pav_events = n_distinct(SV_ID),
    min_distance_to_gene_bp = min(distance_to_gene_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(expr_key, by = c("target_group", "gene_id")) %>%
  mutate(feature_category = factor(feature_category,
                                   levels = c("exon", "intron", "2-kb downstream", "2-kb upstream")))

write.table(membership_priority,
            file.path(tables_dir, "PAV_gene_structure_expression_divergence.priority_category.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

multi_category_audit <- membership_expr_multi %>%
  group_by(target_group, gene_id) %>%
  summarise(n_feature_categories = n_distinct(feature_category), .groups = "drop") %>%
  count(target_group, n_feature_categories, name = "n_gene_species")
write.table(multi_category_audit,
            file.path(tables_dir, "PAV_gene_structure_expression_divergence.multicategory_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(20260705)
boot_n <- as.integer(Sys.getenv("BOOT_N", "10000"))
bootstrap_ci <- function(x, n = boot_n) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(boot_mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  b <- replicate(n, mean(sample(x, size = length(x), replace = TRUE)))
  c(boot_mean = mean(x), ci_low = unname(quantile(b, 0.025)), ci_high = unname(quantile(b, 0.975)))
}

summary_pooled <- membership_priority %>%
  group_by(feature_category) %>%
  summarise(
    n_gene_species = n(),
    n_genes = n_distinct(gene_id),
    mean_divergence = mean(expression_divergence),
    median_divergence = median(expression_divergence),
    q25 = quantile(expression_divergence, 0.25),
    q75 = quantile(expression_divergence, 0.75),
    .groups = "drop"
  )
boot_pooled <- do.call(rbind, lapply(levels(membership_priority$feature_category), function(cat) {
  vals <- membership_priority$expression_divergence[membership_priority$feature_category == cat]
  ci <- bootstrap_ci(vals)
  data.frame(feature_category = cat, t(ci), row.names = NULL)
}))
summary_pooled <- left_join(summary_pooled, boot_pooled, by = "feature_category") %>%
  mutate(bootstrap_n = boot_n)

write.table(summary_pooled,
            file.path(tables_dir, "PAV_gene_structure_expression_divergence.summary_pooled.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_by_species <- membership_priority %>%
  group_by(comparison, target_group, feature_category) %>%
  summarise(
    n_gene_species = n(),
    n_genes = n_distinct(gene_id),
    mean_divergence = mean(expression_divergence),
    median_divergence = median(expression_divergence),
    q25 = quantile(expression_divergence, 0.25),
    q75 = quantile(expression_divergence, 0.75),
    .groups = "drop"
  )
write.table(summary_by_species,
            file.path(tables_dir, "PAV_gene_structure_expression_divergence.summary_by_species.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pairwise_pooled <- as.data.frame(pairwise.wilcox.test(
  membership_priority$expression_divergence,
  membership_priority$feature_category,
  p.adjust.method = "BH",
  exact = FALSE
)$p.value)
pairwise_pooled$feature_category_2 <- rownames(pairwise_pooled)
pairwise_pooled <- pairwise_pooled %>%
  pivot_longer(cols = -feature_category_2, names_to = "feature_category_1", values_to = "p_adj_BH") %>%
  filter(!is.na(p_adj_BH)) %>%
  select(feature_category_1, feature_category_2, p_adj_BH)
write.table(pairwise_pooled,
            file.path(tables_dir, "PAV_gene_structure_expression_divergence.pairwise_wilcox_pooled.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

format_p <- function(p) {
  ifelse(is.na(p), "P = NA",
         ifelse(p < 2.2e-16, "P < 2.2e-16",
                paste0("P = ", format(p, digits = 2, scientific = TRUE))))
}

get_p <- function(cat1, cat2) {
  x <- pairwise_pooled %>%
    filter((feature_category_1 == cat1 & feature_category_2 == cat2) |
             (feature_category_1 == cat2 & feature_category_2 == cat1))
  if (nrow(x) == 0) return(NA_real_)
  x$p_adj_BH[1]
}

ymax <- quantile(membership_priority$expression_divergence, 0.995, na.rm = TRUE)
brackets <- data.frame(
  x = c(1, 1, 1),
  xend = c(2, 3, 4),
  y = ymax * c(1.05, 1.16, 1.27),
  label = c(format_p(get_p("exon", "intron")),
            format_p(get_p("exon", "2-kb downstream")),
            format_p(get_p("exon", "2-kb upstream")))
)

palette_fill <- c(
  "exon" = "#F8766D",
  "intron" = "#B79F00",
  "2-kb downstream" = "#00BA38",
  "2-kb upstream" = "#00BFC4"
)

base_theme <- theme_classic(base_size = 8.5) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text.x = element_text(angle = 36, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8),
    strip.text = element_text(size = 7.5, face = "bold"),
    panel.grid = element_blank()
  )

p_pooled <- ggplot(membership_priority, aes(x = feature_category, y = expression_divergence)) +
  geom_violin(aes(fill = feature_category), color = "grey25", linewidth = 0.35, width = 0.88,
              scale = "width", trim = FALSE, alpha = 0.18) +
  geom_boxplot(aes(color = feature_category), width = 0.20, fill = "white", linewidth = 0.45,
               outlier.size = 0.25, outlier.alpha = 0.35) +
  geom_segment(data = brackets, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = brackets, aes(x = x, xend = x, y = y, yend = y - ymax * 0.025),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = brackets, aes(x = xend, xend = xend, y = y, yend = y - ymax * 0.025),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = brackets, aes(x = (x + xend) / 2, y = y + ymax * 0.025, label = label),
            inherit.aes = FALSE, size = 2.7) +
  scale_fill_manual(values = palette_fill) +
  scale_color_manual(values = palette_fill) +
  coord_cartesian(ylim = c(0, max(brackets$y) * 1.10), clip = "off") +
  labs(
    x = NULL,
    y = expression("Expression divergence  |log"[2]*"(mean TPM + 1)|"),
    title = "Expression divergence by PAV position relative to genes",
    subtitle = "Priority category per gene-species pair: exon > intron > 2-kb upstream > 2-kb downstream"
  ) +
  base_theme

p_by_species <- ggplot(membership_priority, aes(x = feature_category, y = expression_divergence)) +
  geom_violin(aes(fill = feature_category), color = "grey25", linewidth = 0.25, width = 0.88,
              scale = "width", trim = FALSE, alpha = 0.18) +
  geom_boxplot(aes(color = feature_category), width = 0.18, fill = "white", linewidth = 0.35,
               outlier.size = 0.12, outlier.alpha = 0.20) +
  facet_wrap(~ comparison, ncol = 4) +
  scale_fill_manual(values = palette_fill) +
  scale_color_manual(values = palette_fill) +
  labs(
    x = NULL,
    y = expression("Expression divergence  |log"[2]*"(mean TPM + 1)|"),
    title = "Expression divergence by PAV genic position across comparisons"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.2))

save_plot <- function(plot, stem, width_in, height_in, dpi = 600) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in",
                res = dpi, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in",
                 res = dpi, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

save_plot(p_pooled,
          file.path(figures_dir, "PAV_gene_structure_expression_divergence.pooled"),
          width_in = 4.4, height_in = 4.2)
save_plot(p_by_species,
          file.path(figures_dir, "PAV_gene_structure_expression_divergence.by_species"),
          width_in = 7.6, height_in = 5.8)

message("Wrote PAV gene-structure expression divergence figures and tables")

