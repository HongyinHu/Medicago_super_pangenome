suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(forcats)
})

out_dir <- "path/to/project/39.TE_type_soloLTR/output/publication_te_broad_classes"
table_path <- file.path(out_dir, "tables", "edta_broad_class_percent.tsv")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(table_path, show_col_types = FALSE) %>%
  mutate(
    category = factor(category, levels = c(
      "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA transposons", "Helitron", "Other / unclassified"
    )),
    label = if_else(is.na(label) | label == "", species, label),
    percent_of_genome = as.numeric(percent_of_genome),
    species_order = as.integer(species_order)
  )

totals <- df %>%
  group_by(species, label, species_order) %>%
  summarise(total_percent = sum(percent_of_genome), .groups = "drop")

# A restrained, colorblind-aware palette: LTRs stay in cool families, DNA/Helitron in warm families, residual class muted.
pal <- c(
  "LTR/Copia" = "#8F79A5",
  "LTR/Gypsy" = "#2F7F8F",
  "LTR/unknown" = "#86A77B",
  "DNA transposons" = "#D6B84A",
  "Helitron" = "#E79757",
  "Other / unclassified" = "#C7636B"
)

base_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#303030"),
    axis.ticks = element_line(linewidth = 0.30, colour = "#303030"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 6.8, colour = "#2B2B2B"),
    axis.text.y = element_text(size = 6.2, colour = "#1F1F1F", face = "italic"),
    axis.title.x = element_text(size = 7.2, colour = "#1F1F1F", margin = margin(t = 5)),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.2, colour = "#1F1F1F"),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width = unit(3.5, "mm"),
    legend.margin = margin(l = 4, r = 0, unit = "pt"),
    panel.grid.major.x = element_line(linewidth = 0.18, colour = "#E6E6E6"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8.5, face = "bold", colour = "#1F1F1F", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.4, colour = "#5C5C5C", margin = margin(b = 5)),
    plot.margin = margin(5, 5, 5, 5)
  )

make_plot <- function(dat, totals_dat, label_levels, title, subtitle) {
  dat <- dat %>% mutate(label_plot = factor(label, levels = rev(label_levels)))
  totals_dat <- totals_dat %>% mutate(label_plot = factor(label, levels = rev(label_levels)))
  xmax <- ceiling(max(totals_dat$total_percent, na.rm = TRUE) / 10) * 10
  xmax <- max(70, xmax)
  ggplot(dat, aes(x = percent_of_genome, y = label_plot, fill = category)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.12) +
    geom_text(
      data = totals_dat,
      aes(x = total_percent + xmax * 0.012, y = label_plot, label = sprintf("%.1f", total_percent)),
      inherit.aes = FALSE,
      hjust = 0,
      size = 1.9,
      family = "Arial",
      colour = "#3A3A3A"
    ) +
    scale_fill_manual(values = pal, breaks = names(pal)) +
    scale_x_continuous(
      limits = c(0, xmax * 1.09),
      breaks = seq(0, xmax, by = 10),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Percent of assembly (%)"
    ) +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(reverse = FALSE, byrow = TRUE)) +
    base_theme
}

save_plot <- function(p, prefix, width_mm = 176, height_mm = 118, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  pdf_file <- paste0(prefix, ".pdf")
  png_file <- paste0(prefix, ".png")
  tiff_file <- paste0(prefix, ".tiff")
  svg_file <- paste0(prefix, ".svg")
  grDevices::cairo_pdf(pdf_file, width = w, height = h, family = "Arial")
  print(p)
  grDevices::dev.off()
  ragg::agg_png(png_file, width = w, height = h, units = "in", res = dpi, background = "white")
  print(p)
  grDevices::dev.off()
  ragg::agg_tiff(tiff_file, width = w, height = h, units = "in", res = dpi, background = "white", compression = "lzw")
  print(p)
  grDevices::dev.off()
  grDevices::svg(svg_file, width = w, height = h, family = "Arial")
  print(p)
  grDevices::dev.off()
}

levels_project <- totals %>% arrange(species_order, species) %>% pull(label)
plot_project <- make_plot(
  df, totals, levels_project,
  "Transposable element composition across Medicago genomes",
  "EDTA annotations grouped into broad TE classes; values denote percent of assembled genome."
)
save_plot(plot_project, file.path(fig_dir, "edta_broad_te_composition.project_order"))

levels_total <- totals %>% arrange(total_percent, species_order, species) %>% pull(label)
plot_total <- make_plot(
  df, totals, levels_total,
  "Transposable element composition across Medicago genomes",
  "Samples sorted by total annotated TE fraction."
)
save_plot(plot_total, file.path(fig_dir, "edta_broad_te_composition.sorted_by_total"))

write_tsv(totals %>% arrange(desc(total_percent)), file.path(out_dir, "tables", "edta_broad_class_totals.tsv"))
cat("Wrote figures to ", fig_dir, "\n", sep = "")
