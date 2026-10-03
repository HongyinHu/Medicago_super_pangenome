#!/usr/bin/env Rscript

base_dir <- "path/to/project/N_4.pod_spiny/15_Chr23997_Msa_reference_similarity_RBH_20260708"
data_dir <- file.path(base_dir, "data")
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(file.path(data_dir, "samples.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
profile <- read.delim(file.path(data_dir, "similarity_profile.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
regions <- read.delim(file.path(data_dir, "target_regions.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
gene <- read.delim(file.path(data_dir, "gene_model.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
events <- read.delim(file.path(data_dir, "event_regions.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

profile$x <- as.numeric(profile$x)
profile$similarity <- as.numeric(profile$similarity)
samples$plot_order <- as.numeric(samples$plot_order)

plot_samples <- samples[samples$sample %in% unique(profile$sample), ]
plot_samples <- plot_samples[order(plot_samples$plot_order), ]

ref_len <- 88982331 - 88979818 + 1
spiny_col <- "#F27678"
spineless_col <- "#6CBBC1"
candidate_col <- "#9B9B9B"
grey_line <- "#3A3A3A"
light_grey <- "#EFEFEF"
mid_grey <- "#D4D4D4"
out_prefix <- "Chr23997_Msa_reference_sequence_similarity_RBH"

draw_tracks <- function() {
  n <- nrow(plot_samples)
  track_gap <- 0.70
  track_h <- 0.44
  bottom_y <- 0.80
  y_pos <- bottom_y + rev(seq_len(n)) * track_gap
  names(y_pos) <- plot_samples$sample
  top_track <- max(y_pos) + track_h + 0.15
  gene_y <- top_track + 0.85
  y_max <- gene_y + 0.80

  par(mar = c(3.8, 4.5, 1.1, 4.3), xaxs = "i", yaxs = "i", family = "sans")
  plot.new()
  plot.window(xlim = c(-380, ref_len + 1100), ylim = c(0, y_max))

  intron <- events[events$event == "second_intron", ]
  core <- events[events$event == "core_INS_DEL", ]
  rect(intron$start, bottom_y - 0.20, intron$end, gene_y + 0.22, col = light_grey, border = NA)
  rect(core$start, bottom_y - 0.20, core$end, gene_y + 0.22, col = adjustcolor(mid_grey, alpha.f = 0.60), border = NA)

  segments(1, gene_y, ref_len, gene_y, col = grey_line, lwd = 1.2)
  arrows(ref_len, gene_y, ref_len - 170, gene_y, length = 0.08, lwd = 1.2, col = grey_line)
  for (i in seq_len(nrow(gene))) {
    rect(gene$start[i], gene_y - 0.11, gene$end[i], gene_y + 0.11, col = "#4A4A4A", border = NA)
  }
  text(ref_len / 2, gene_y + 0.56, "Chr23997 (genome_Msa reference)", cex = 0.90, font = 3)
  text(1, gene_y - 0.32, "Chr4 88,979,818", adj = c(0, 0.5), cex = 0.70)
  text(ref_len, gene_y - 0.32, "88,982,331", adj = c(1, 0.5), cex = 0.70)
  text((core$start + core$end) / 2, gene_y + 0.18, "222-bp INS/DEL", cex = 0.56, col = "#555555")

  group_names <- c("spiny", "spineless", "candidate")
  for (g_i in seq_len(length(group_names) - 1)) {
    idx_a <- which(plot_samples$group == group_names[g_i])
    idx_b <- which(plot_samples$group == group_names[g_i + 1])
    if (length(idx_a) > 0 && length(idx_b) > 0) {
      sep_y <- mean(c(y_pos[plot_samples$sample[max(idx_a)]], y_pos[plot_samples$sample[min(idx_b)]])) - 0.02
      segments(1, sep_y, ref_len, sep_y, col = "#BEBEBE", lwd = 0.8)
    }
  }

  for (i in seq_len(nrow(plot_samples))) {
    sample <- plot_samples$sample[i]
    dat <- profile[profile$sample == sample, ]
    dat <- dat[order(dat$x), ]
    y0 <- y_pos[sample]
    col <- if (plot_samples$group[i] == "spiny") spiny_col else if (plot_samples$group[i] == "spineless") spineless_col else candidate_col
    segments(1, y0, ref_len, y0, col = "#222222", lwd = 0.55)
    if (nrow(dat) > 1) {
      yy <- y0 + track_h * pmax(0, pmin(100, dat$similarity)) / 100
      polygon(c(dat$x[1], dat$x, dat$x[nrow(dat)]), c(y0, yy, y0), col = adjustcolor(col, alpha.f = 0.88), border = NA)
      lines(dat$x, yy, col = adjustcolor(col, alpha.f = 0.95), lwd = 0.55)
    }
    label <- sample
    if (plot_samples$confidence[i] == "best_hit_non_RBH") {
      label <- paste0(sample, "*")
    }
    text(ref_len + 35, y0 + track_h * 0.52, label, adj = c(0, 0.5), cex = 0.73)
  }

  spiny_idx <- which(plot_samples$group == "spiny")
  spineless_idx <- which(plot_samples$group == "spineless")
  candidate_idx <- which(plot_samples$group == "candidate")
  if (length(spiny_idx) > 0) {
    text(-185, mean(y_pos[plot_samples$sample[spiny_idx]]) + track_h * 0.45, "Spiny\nstrict RBH", srt = 90, cex = 0.68, col = spiny_col, font = 2)
  }
  if (length(spineless_idx) > 0) {
    text(-185, mean(y_pos[plot_samples$sample[spineless_idx]]) + track_h * 0.45, "Spineless\nstrict RBH", srt = 90, cex = 0.68, col = spineless_col, font = 2)
  }
  if (length(candidate_idx) > 0) {
    text(-185, mean(y_pos[plot_samples$sample[candidate_idx]]) + track_h * 0.45, "Best-hit\ncandidate", srt = 90, cex = 0.66, col = candidate_col, font = 2)
  }

  text(-315, mean(range(y_pos)) + track_h * 0.5, "Nucleotide sequence similarity (%)", srt = 90, cex = 0.92)
  axis(1, at = seq(0, 2400, 300), labels = seq(0, 2400, 300), pos = bottom_y - 0.33, lwd = 1, lwd.ticks = 1, cex.axis = 0.76)
  mtext("Position on genome_Msa Chr23997 region (bp)", side = 1, line = 2.3, cex = 0.83)

  legend_x <- ref_len + 250
  legend_y <- gene_y + 0.30
  rect(legend_x, legend_y - 0.12, legend_x + 55, legend_y + 0.08, col = spiny_col, border = NA)
  text(legend_x + 70, legend_y - 0.02, "Spiny strict RBH", adj = c(0, 0.5), cex = 0.62)
  rect(legend_x, legend_y - 0.42, legend_x + 55, legend_y - 0.22, col = spineless_col, border = NA)
  text(legend_x + 70, legend_y - 0.32, "Spineless strict RBH", adj = c(0, 0.5), cex = 0.62)
  rect(legend_x, legend_y - 0.72, legend_x + 55, legend_y - 0.52, col = candidate_col, border = NA)
  text(legend_x + 70, legend_y - 0.62, "Best-hit candidate", adj = c(0, 0.5), cex = 0.62)
  text(legend_x, legend_y - 1.02, "* non-RBH candidate", adj = c(0, 0.5), cex = 0.52, col = "#555555")
}

height_in <- max(7.2, 2.2 + nrow(plot_samples) * 0.38)
width_in <- 8.2

svg(file.path(fig_dir, paste0(out_prefix, ".svg")), width = width_in, height = height_in)
draw_tracks()
dev.off()

pdf(file.path(fig_dir, paste0(out_prefix, ".pdf")), width = width_in, height = height_in, useDingbats = FALSE)
draw_tracks()
dev.off()

png(file.path(fig_dir, paste0(out_prefix, ".png")), width = width_in, height = height_in, units = "in", res = 350)
draw_tracks()
dev.off()

try({
  tiff(file.path(fig_dir, paste0(out_prefix, ".tiff")), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
  draw_tracks()
  dev.off()
}, silent = TRUE)

skipped <- regions[!(regions$sample %in% plot_samples$sample), c("sample", "group", "status", "note")]
write.table(skipped, file.path(fig_dir, paste0(out_prefix, ".skipped.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Plotted samples:", nrow(plot_samples), "\n")
cat("Output dir:", fig_dir, "\n")
