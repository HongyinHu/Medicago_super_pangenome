suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

out_dir <- "path/to/project/39.TE_type_soloLTR/output/publication_te_broad_classes"
fig_dir <- file.path(out_dir, "figures")
df <- read_tsv(file.path(out_dir, "tables", "edta_broad_class_percent.tsv"), show_col_types = FALSE) %>%
  mutate(
    percent_of_genome = as.numeric(percent_of_genome),
    species_order = as.integer(species_order),
    label = if_else(is.na(label) | label == "", species, label)
  )

totals <- df %>%
  group_by(species, label, species_order) %>%
  summarise(total_percent = sum(percent_of_genome), .groups = "drop")

# Factor levels are right-to-left for ggplot's default stack, so the visual bar reads:
# Other -> Helitron -> DNA -> LTR/unknown -> LTR/Gypsy -> LTR/Copia.
stack_levels <- c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA transposons", "Helitron", "Other / unclassified")
legend_breaks <- rev(stack_levels)
pal <- c(
  "Other / unclassified" = "#C7636B",
  "Helitron" = "#E79757",
  "DNA transposons" = "#D6B84A",
  "LTR/unknown" = "#86A77B",
  "LTR/Gypsy" = "#2F7F8F",
  "LTR/Copia" = "#8F79A5"
)

df <- df %>% mutate(category = factor(category, levels = stack_levels))

plot_one <- function(levels_y, prefix, width_mm = 170, height_mm = 112) {
  plot_df <- df %>% mutate(label_plot = factor(label, levels = rev(levels_y)))
  total_df <- totals %>% mutate(label_plot = factor(label, levels = rev(levels_y)))
  xmax <- max(90, ceiling(max(total_df$total_percent, na.rm = TRUE) / 10) * 10)
  p <- ggplot(plot_df, aes(x = percent_of_genome, y = label_plot, fill = category)) +
    geom_col(width = 0.70, colour = "white", linewidth = 0.10) +
    geom_text(
      data = total_df,
      aes(x = total_percent + 1.0, y = label_plot, label = sprintf("%.1f", total_percent)),
      inherit.aes = FALSE,
      hjust = 0,
      size = 1.75,
      family = "Arial",
      colour = "#4B4B4B"
    ) +
    scale_fill_manual(values = pal, breaks = legend_breaks) +
    scale_x_continuous(
      limits = c(0, xmax + 7),
      breaks = seq(0, xmax, by = 10),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = "Percent of assembly (%)", y = NULL) +
    guides(fill = guide_legend(title = NULL, override.aes = list(linewidth = 0), keyheight = unit(3.2, "mm"), keywidth = unit(3.2, "mm"))) +
    theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = "#333333"),
      axis.ticks.x = element_line(linewidth = 0.25, colour = "#333333"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6.3, colour = "#2B2B2B"),
      axis.text.y = element_text(size = 5.9, colour = "#202020"),
      axis.title.x = element_text(size = 6.9, margin = margin(t = 4)),
      legend.position = "right",
      legend.text = element_text(size = 5.8, colour = "#202020"),
      legend.margin = margin(l = 2, r = 0, unit = "pt"),
      panel.grid.major.x = element_line(linewidth = 0.16, colour = "#E9E9E9"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 4, 3, 3)
    )
  w <- width_mm / 25.4; h <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = w, height = h, family = "Arial")
  print(p); grDevices::dev.off()
  ragg::agg_png(paste0(prefix, ".png"), width = w, height = h, units = "in", res = 700, background = "white")
  print(p); grDevices::dev.off()
  ragg::agg_tiff(paste0(prefix, ".tiff"), width = w, height = h, units = "in", res = 700, background = "white", compression = "lzw")
  print(p); grDevices::dev.off()
  grDevices::svg(paste0(prefix, ".svg"), width = w, height = h, family = "Arial")
  print(p); grDevices::dev.off()
}

levels_project <- totals %>% arrange(species_order, species) %>% pull(label)
plot_one(levels_project, file.path(fig_dir, "edta_broad_te_composition.manuscript_project_order"))
levels_total <- totals %>% arrange(total_percent, species_order, species) %>% pull(label)
plot_one(levels_total, file.path(fig_dir, "edta_broad_te_composition.manuscript_sorted_by_total"))
