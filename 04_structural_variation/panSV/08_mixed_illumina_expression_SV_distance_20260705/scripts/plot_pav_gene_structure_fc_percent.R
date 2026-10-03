suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: plot_pav_gene_structure_fc_percent.R WORK")
}

work <- args[[1]]
tables_dir <- file.path(work, "tables")
figures_dir <- file.path(work, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

infile <- file.path(tables_dir, "PAV_gene_structure_expression_divergence.priority_category.tsv")
if (!file.exists(infile)) {
  stop("Missing input table: ", infile)
}

df <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE)

plot_df <- df %>%
  mutate(
    abs_log2FC = abs(log2FC_target_vs_ref),
    fold_change = 2 ^ abs_log2FC,
    fc_bin = case_when(
      fold_change < 1.5 ~ "1-1.5",
      fold_change < 2 ~ "1.5-2",
      fold_change < 3 ~ "2-3",
      TRUE ~ ">3"
    ),
    feature_label = case_when(
      feature_category == "exon" ~ "Coding region",
      feature_category == "intron" ~ "Intron",
      feature_category == "2-kb upstream" ~ "2 kb upstream",
      feature_category == "2-kb downstream" ~ "2 kb downstream",
      TRUE ~ feature_category
    )
  ) %>%
  mutate(
    feature_label = factor(feature_label,
                           levels = c("Coding region", "2 kb upstream", "2 kb downstream", "Intron")),
    fc_bin = factor(fc_bin, levels = c("1-1.5", "1.5-2", "2-3", ">3"))
  )

summary_count <- plot_df %>%
  count(feature_label, fc_bin, name = "n_gene_species") %>%
  group_by(feature_label) %>%
  mutate(
    total_gene_species = sum(n_gene_species),
    percent = 100 * n_gene_species / total_gene_species
  ) %>%
  ungroup()

write.table(summary_count,
            file.path(tables_dir, "PAV_gene_structure_foldchange_percent.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_by_species <- plot_df %>%
  count(comparison, target_group, feature_label, fc_bin, name = "n_gene_species") %>%
  group_by(comparison, target_group, feature_label) %>%
  mutate(
    total_gene_species = sum(n_gene_species),
    percent = 100 * n_gene_species / total_gene_species
  ) %>%
  ungroup()

write.table(summary_by_species,
            file.path(tables_dir, "PAV_gene_structure_foldchange_percent.by_species.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

palette_fc <- c(
  "1-1.5" = "#23A891",
  "1.5-2" = "#63BDD0",
  "2-3" = "#FFD447",
  ">3" = "#9B92C7"
)

base_theme <- theme_classic(base_size = 8.5) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(3.2, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 8, 5, 5)
  )

p <- ggplot(summary_count, aes(x = feature_label, y = percent, fill = fc_bin)) +
  geom_col(width = 0.62, color = NA) +
  scale_fill_manual(values = palette_fc, name = "Fold change") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Percent (%)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  base_theme

p_by_species <- ggplot(summary_by_species, aes(x = feature_label, y = percent, fill = fc_bin)) +
  geom_col(width = 0.62, color = NA) +
  facet_wrap(~ comparison, ncol = 4) +
  scale_fill_manual(values = palette_fc, name = "Fold change") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Percent (%)") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 7, face = "bold"))

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

save_plot(p, file.path(figures_dir, "PAV_gene_structure_foldchange_percent"),
          width_in = 2.4, height_in = 4.2)
save_plot(p_by_species, file.path(figures_dir, "PAV_gene_structure_foldchange_percent.by_species"),
          width_in = 7.6, height_in = 5.8)

message("Wrote fold-change percent figures and source tables")

