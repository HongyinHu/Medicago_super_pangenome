suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

out_dir <- "path/to/project/39.TE_type_soloLTR/output/publication_ltr_insertion_time"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
records <- read_tsv(file.path(out_dir, "tables", "ltr_insertion_time_records.tsv"), show_col_types = FALSE) %>%
  mutate(
    insertion_time_mya = as.numeric(insertion_time_mya),
    species_order = as.integer(species_order),
    label = if_else(is.na(label) | label == "", species, label),
    ltr_group = recode(ltr_class, Gypsy = "LTRs-Gypsy", Copia = "LTRs-Copia", Unknown = "LTRs-unknown")
  ) %>%
  filter(!is.na(insertion_time_mya), insertion_time_mya >= 0, insertion_time_mya <= 8)

bin_width <- 0.25
panel_levels <- c("LTRs-Gypsy", "LTRs-Copia", "LTRs-unknown", "LTRs")
panel_labels <- c("a  LTRs-Gypsy", "b  LTRs-Copia", "c  LTRs-unknown", "d  LTRs")

make_binned <- function(dat) {
  by_class <- dat %>%
    mutate(time_bin = floor(insertion_time_mya / bin_width) * bin_width + bin_width / 2) %>%
    count(species, label, species_order, ltr_group, time_bin, name = "count") %>%
    rename(panel = ltr_group)
  all_ltr <- dat %>%
    mutate(time_bin = floor(insertion_time_mya / bin_width) * bin_width + bin_width / 2) %>%
    count(species, label, species_order, time_bin, name = "count") %>%
    mutate(panel = "LTRs")
  bind_rows(by_class, all_ltr) %>%
    mutate(panel = factor(panel, levels = panel_levels, labels = panel_labels))
}

palette_main <- c(
  "genome_395" = "#1f77b4", "genome_436" = "#ff7f0e", "genome_457" = "#2ca02c",
  "genome_468" = "#d62728", "genome_474" = "#9467bd", "genome_482" = "#8c564b",
  "genome_M22" = "#e377c2", "Marc" = "#7f7f7f", "Mrut" = "#bcbd22",
  "Msat_Cae" = "#17becf", "Msat_ZM4" = "#111111", "genome_410" = "#6baed6",
  "genome_454" = "#fdae6b", "genome_461" = "#74c476", "genome_472" = "#fb6a4a",
  "genome_M46" = "#9e9ac8", "Mpol" = "#c49c94", "Mtru_R108" = "#f7b6d2",
  "genome_474_T2T" = "#636363", "Mtru_A17" = "#bdbdbd", "Msat_T2T" = "#969696"
)

plot_dataset <- function(dat, prefix, exclude_note = NULL) {
  binned <- make_binned(dat)
  levels_label <- dat %>% distinct(label, species_order) %>% arrange(species_order, label) %>% pull(label)
  binned <- binned %>% mutate(label = factor(label, levels = levels_label))
  pal <- palette_main[levels_label]
  missing_cols <- is.na(pal)
  if (any(missing_cols)) {
    fallback <- grDevices::hcl.colors(sum(missing_cols), palette = "Dark 3")
    pal[missing_cols] <- fallback
  }
  names(pal) <- levels_label
  size_map <- rep(0.55, length(levels_label)); names(size_map) <- levels_label
  if ("Msat_ZM4" %in% names(size_map)) size_map["Msat_ZM4"] <- 0.95
  alpha_map <- rep(0.82, length(levels_label)); names(alpha_map) <- levels_label
  if ("Msat_ZM4" %in% names(alpha_map)) alpha_map["Msat_ZM4"] <- 1.0

  p <- ggplot(binned, aes(x = time_bin, y = count, colour = label, group = label)) +
    geom_point(aes(alpha = label), size = 0.55, stroke = 0, show.legend = FALSE) +
    geom_smooth(aes(linewidth = label, alpha = label), method = "loess", formula = y ~ x,
                span = 0.28, se = FALSE, show.legend = TRUE) +
    facet_wrap(~panel, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = pal, name = NULL) +
    scale_linewidth_manual(values = size_map, guide = "none") +
    scale_alpha_manual(values = alpha_map, guide = "none") +
    scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 1), expand = expansion(mult = c(0, 0.01))) +
    labs(x = "Insertion time (Mya)", y = "Number of intact LTR-RTs") +
    guides(colour = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(linewidth = 1.0, alpha = 1))) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#262626"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#262626"),
      axis.text = element_text(size = 6.2, colour = "#222222"),
      axis.title = element_text(size = 7.0, colour = "#111111"),
      strip.background = element_blank(),
      strip.text = element_text(size = 7.2, face = "plain", colour = "#111111", margin = margin(b = 3)),
      panel.spacing.x = unit(11, "mm"),
      panel.spacing.y = unit(8, "mm"),
      legend.position = "bottom",
      legend.text = element_text(size = 5.4, colour = "#111111"),
      legend.key.width = unit(6, "mm"),
      legend.key.height = unit(2.8, "mm"),
      legend.margin = margin(t = 0, b = 0),
      plot.margin = margin(4, 5, 3, 5),
      panel.grid.major.y = element_line(linewidth = 0.14, colour = "#ECECEC"),
      panel.grid.minor = element_blank()
    )
  w <- 178 / 25.4; h <- 145 / 25.4
  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = w, height = h, family = "Arial"); print(p); grDevices::dev.off()
  ragg::agg_png(paste0(prefix, ".png"), width = w, height = h, units = "in", res = 700, background = "white"); print(p); grDevices::dev.off()
  ragg::agg_tiff(paste0(prefix, ".tiff"), width = w, height = h, units = "in", res = 700, background = "white", compression = "lzw"); print(p); grDevices::dev.off()
  grDevices::svg(paste0(prefix, ".svg"), width = w, height = h, family = "Arial"); print(p); grDevices::dev.off()

  write_tsv(binned %>% arrange(panel, species_order, label, time_bin), paste0(prefix, ".binned_counts.tsv"))
}

main_exclude <- c("genome_474_T2T", "genome_A17", "genome_Msa_T2T")
plot_dataset(records %>% filter(!species %in% main_exclude), file.path(fig_dir, "ltr_insertion_time.main18_with_Msat_ZM4"))
plot_dataset(records, file.path(fig_dir, "ltr_insertion_time.all21_with_Msat_ZM4"))

summary <- records %>%
  count(species, label, species_order, ltr_class, name = "intact_ltr_count") %>%
  arrange(species_order, species, ltr_class)
write_tsv(summary, file.path(out_dir, "tables", "ltr_insertion_time_class_counts.tsv"))
