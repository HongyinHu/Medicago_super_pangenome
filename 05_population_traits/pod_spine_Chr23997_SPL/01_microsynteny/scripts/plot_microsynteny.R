args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: plot_microsynteny.R <out_dir> <plot_prefix>")
}
out_dir <- args[[1]]
prefix <- args[[2]]

genes <- read.delim(file.path(out_dir, "neighboring_gene_table.tsv"), stringsAsFactors = FALSE)
links <- read.delim(file.path(out_dir, "synteny_links.tsv"), stringsAsFactors = FALSE)

row_y <- c(Msa = 3, Ath = 2, R108 = 1)
row_lab <- c(Msa = "Medicago Msa", Ath = "Arabidopsis SPL10", R108 = "Medicago R108")
genes$x <- as.numeric(genes$Order_index)
genes$y <- row_y[genes$Species]

palette <- c(
  "#4E79A7", "#59A14F", "#F28E2B", "#76B7B2", "#EDC948",
  "#B07AA1", "#9C755F", "#BAB0AC", "#86BCB6", "#E15759"
)
orthogroups <- unique(links$Orthogroup[links$Orthogroup != "SPL10_candidate"])
col_map <- setNames(palette[seq_along(orthogroups) %% length(palette) + 1], orthogroups)
col_map["SPL10_candidate"] <- "#C43B3B"
default_col <- "#D9D9D9"

gene_key <- paste(links$Species, links$Gene_ID, sep = "|")
og_by_gene <- setNames(links$Orthogroup, gene_key)
genes$key <- paste(genes$Species, genes$Gene_ID, sep = "|")
genes$og <- og_by_gene[genes$key]
genes$fill <- ifelse(is.na(genes$og), default_col, col_map[genes$og])

draw_panel <- function() {
  par(mar = c(2.2, 7.2, 1.0, 1.0), family = "sans")
  plot(NA, xlim = c(-1.55, 11.75), ylim = c(0.45, 3.55), axes = FALSE, xlab = "", ylab = "")
  for (sp in names(row_y)) {
    y <- row_y[sp]
    segments(0.55, y, 11.3, y, col = "#BDBDBD", lwd = 0.8)
    text(-1.48, y, row_lab[sp], adj = 0, cex = 0.72, font = ifelse(sp == "Ath", 2, 1))
  }

  link_pairs <- subset(links, Orthogroup != "SPL10_candidate")
  for (og in unique(link_pairs$Orthogroup)) {
    rows <- subset(link_pairs, Orthogroup == og)
    ath <- subset(rows, Species == "Ath")
    if (nrow(ath) == 0) next
    ath_gene <- ath$Gene_ID[1]
    g_ath <- subset(genes, Species == "Ath" & Gene_ID == ath_gene)
    if (nrow(g_ath) == 0) next
    for (sp in c("Msa", "R108")) {
      q <- subset(rows, Species == sp)
      if (nrow(q) == 0) next
      g_q <- subset(genes, Species == sp & Gene_ID == q$Gene_ID[1])
      if (nrow(g_q) == 0) next
      segments(g_q$x[1], g_q$y[1], g_ath$x[1], g_ath$y[1],
               col = adjustcolor(col_map[og], alpha.f = 0.35), lwd = 2)
    }
  }

  draw_gene <- function(x, y, strand, fill, border, lwd) {
    w <- 0.78
    h <- 0.17
    head <- 0.18
    if (strand == "-") {
      pts_x <- c(x + w / 2, x - w / 2 + head, x - w / 2, x - w / 2 + head, x + w / 2)
      pts_y <- c(y - h / 2, y - h / 2, y, y + h / 2, y + h / 2)
    } else {
      pts_x <- c(x - w / 2, x + w / 2 - head, x + w / 2, x + w / 2 - head, x - w / 2)
      pts_y <- c(y - h / 2, y - h / 2, y, y + h / 2, y + h / 2)
    }
    polygon(pts_x, pts_y, col = fill, border = border, lwd = lwd)
  }
  for (i in seq_len(nrow(genes))) {
    g <- genes[i, ]
    draw_gene(g$x, g$y, g$Strand, g$fill, "#333333", ifelse(g$Is_target == "1", 1.4, 0.55))
    label <- ifelse(g$Is_target == "1", g$Display_label, g$Gene_ID)
    text(g$x, g$y + 0.28, label, srt = 35, adj = 0, cex = ifelse(g$Is_target == "1", 0.66, 0.50))
  }
  axis(1, at = 1:11, labels = genes$Gene_order[genes$Species == "Ath"], tick = FALSE, line = -0.8, cex.axis = 0.72)
  box(col = NA)
}

pdf(paste0(prefix, ".pdf"), width = 9.2, height = 4.2, useDingbats = FALSE)
draw_panel()
dev.off()

png(paste0(prefix, ".png"), width = 3200, height = 1500, res = 350, type = "cairo")
draw_panel()
dev.off()
