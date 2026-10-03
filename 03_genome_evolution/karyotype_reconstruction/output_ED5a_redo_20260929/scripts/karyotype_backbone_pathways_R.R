#!/usr/bin/env Rscript

## R-only redraw of the AMK karyotype backbone and the evidence-ranked
## rearrangement pathways.  The pathways are schematic minimum models; they
## do not claim a unique historical order where the WGDI blocks cannot resolve it.

library(grid)

out_dir <- "path/to/project/37.karyotype_reconstruction/output_ED5a_redo_20260929/07_karyotype_evolution"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

AMK_COL <- c("#4A69AD", "#E51B23", "#8DBB1F", "#2FA6D3",
             "#328F65", "#F5C21B", "#A04AA4", "#A5601C")
CL_COL <- c(III = "#B23A55", II = "#E08A2C", I = "#3F8F7A", out = "#555555")
INK <- "#202020"
MUTED <- "#707070"
GREEN <- "#2E7D32"
LIGHT <- "#F7F7F7"

push <- function(xscale = c(0, 16), yscale = c(0, 12)) {
  pushViewport(viewport(xscale = xscale, yscale = yscale,
                        name = "figure", gp = gpar(fontfamily = "Arial")))
}

txt <- function(x, y, label, size = 9, col = INK, face = "plain",
                just = c("centre", "centre"), rot = 0) {
  grid.text(label, x = unit(x, "native"), y = unit(y, "native"),
            just = just, rot = rot,
            gp = gpar(fontfamily = "Arial", fontsize = size, col = col,
                      fontface = face))
}

line <- function(x, y, col = INK, lwd = 1.1, lty = 1, arrow_end = FALSE) {
  grid.lines(unit(x, "native"), unit(y, "native"),
             gp = gpar(col = col, lwd = lwd, lty = lty),
             arrow = if (arrow_end) arrow(type = "closed", length = unit(2, "mm")) else NULL)
}

box <- function(x, y, w, h, col = "#666666", lwd = 1, fill = "white", lty = 1) {
  grid.rect(x = unit(x + w / 2, "native"), y = unit(y + h / 2, "native"),
            width = unit(w, "native"), height = unit(h, "native"),
            gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
}

## x,y are the lower-left corner; segs is a data.frame with colour and frac.
chrom <- function(x, y, w = 0.10, h = 0.55, segs, border = "#777777", lwd = 0.55) {
  if (nrow(segs) == 0) segs <- data.frame(col = "white", frac = 1)
  grid.rect(unit(x + w / 2, "native"), unit(y + h / 2, "native"),
            width = unit(w, "native"), height = unit(h, "native"),
            gp = gpar(fill = LIGHT, col = NA))
  yy <- y
  for (i in seq_len(nrow(segs))) {
    hh <- h * segs$frac[i]
    grid.rect(unit(x + w / 2, "native"), unit(yy + hh / 2, "native"),
              width = unit(w, "native"), height = unit(hh, "native"),
              gp = gpar(fill = segs$col[i], col = "white", lwd = 0.25))
    yy <- yy + hh
  }
  grid.rect(unit(x + w / 2, "native"), unit(y + h / 2, "native"),
            width = unit(w, "native"), height = unit(h, "native"),
            gp = gpar(fill = NA, col = border, lwd = lwd))
}

amk_seg <- function(a) data.frame(col = AMK_COL[a], frac = 1)
mix <- function(a, f) data.frame(col = AMK_COL[a], frac = f)

amk_row <- function(x0, y0, scale = 1, labels = TRUE) {
  for (a in 1:8) {
    chrom(x0 + (a - 1) * 0.19 * scale, y0, w = 0.095 * scale,
          h = 0.48 * scale, amk_seg(a), border = NA)
    if (labels) txt(x0 + (a - 1) * 0.19 * scale + 0.047 * scale,
                    y0 - 0.13 * scale, as.character(a), 6.5 * scale, MUTED)
  }
}

plus <- function(x, y, size = 9) txt(x, y, "+", size = size, col = INK, face = "bold")

pair <- function(x, y, a, b, out = NULL, op = "RCT", scale = 1,
                 labels = TRUE, candidate = FALSE) {
  chrom(x, y, w = 0.11 * scale, h = 0.55 * scale, amk_seg(a), border = NA)
  plus(x + 0.23 * scale, y + 0.27 * scale, 8 * scale)
  chrom(x + 0.32 * scale, y, w = 0.11 * scale, h = 0.55 * scale, amk_seg(b), border = NA)
  line(c(x + 0.48 * scale, x + 0.90 * scale),
       c(y + 0.28 * scale, y + 0.28 * scale), col = INK, lwd = 0.8, arrow_end = TRUE)
  if (op == "RCT") {
    grid.lines(unit(c(x + 0.63, x + 0.80) * scale + x * (1 - scale), "native"),
               unit(c(y + 0.10, y + 0.46) * scale + y * (1 - scale), "native"),
               gp = gpar(col = AMK_COL[a], lwd = 2.0 * scale))
    grid.lines(unit(c(x + 0.63, x + 0.80) * scale + x * (1 - scale), "native"),
               unit(c(y + 0.46, y + 0.10) * scale + y * (1 - scale), "native"),
               gp = gpar(col = AMK_COL[b], lwd = 2.0 * scale))
  } else {
    txt(x + 0.72 * scale, y + 0.28 * scale, op, 6.5 * scale, GREEN,
        face = "bold")
  }
  if (is.null(out)) out <- c(a, b)
  ## outputs are deliberately schematic derivatives, not coordinate-scaled chromosomes
  chrom(x + 1.02 * scale, y, w = 0.11 * scale, h = 0.55 * scale,
        mix(c(a, b), c(0.55, 0.45)), border = NA)
  plus(x + 1.25 * scale, y + 0.27 * scale, 8 * scale)
  chrom(x + 1.34 * scale, y, w = 0.11 * scale, h = 0.55 * scale,
        mix(c(b, a), c(0.55, 0.45)), border = NA)
  if (labels) {
    txt(x + 1.075 * scale, y - 0.13 * scale, paste0("Chr", out[1]), 5.7 * scale, MUTED)
    txt(x + 1.395 * scale, y - 0.13 * scale, paste0("Chr", out[2]), 5.7 * scale, MUTED)
  }
}

draw_header <- function(label, title, x, y, col = INK) {
  txt(x, y, label, 11, col, face = "bold", just = c("left", "top"))
  txt(x + 0.24, y, title, 9.2, INK, face = "bold", just = c("left", "top"))
}

draw_tip <- function(x, name, lab, clade, tip_segments, changed = FALSE, x7 = FALSE) {
  y <- 6.18
  g <- tip_segments[tip_segments$species == lab, , drop = FALSE]
  chs <- sort(unique(as.integer(g$chr)))
  nbar <- if (length(chs)) length(chs) else if (x7) 7 else 8
  slot <- 0.72 / nbar
  for (i in seq_len(nbar)) {
    ch <- if (i <= length(chs)) chs[i] else i
    gg <- g[as.integer(g$chr) == ch, , drop = FALSE]
    total <- if (nrow(gg)) max(gg$chr_genes) else 1
    segs <- if (nrow(gg)) {
      data.frame(col = AMK_COL[as.integer(gg$AMK)], frac = gg$genes / total)
    } else data.frame(col = LIGHT, frac = 1)
    chrom(x + (i - 1) * slot, y, w = slot * 0.55, h = 0.38 * total / max(1, if (nrow(g)) max(g$chr_genes) else 1),
          segs, border = NA)
  }
  if (x7) {
    ## use a compact red outline to flag the seven-chromosome assembly
    box(x - 0.07, y - 0.06, 0.82, 0.56, col = "#C0392B", lwd = 1.0, fill = NA, lty = 2)
  }
  nm <- paste(strwrap(name, width = 17), collapse = "\n")
  txt(x + 0.33, 5.78, nm, 5.0, INK, face = "italic")
  txt(x + 0.33, 5.33, if (x7) "x = 7" else "x = 8", 5.0,
      if (x7) "#C0392B" else MUTED, face = "bold")
}

draw_backbone <- function() {
  tip_segments <- read.delim(file.path(out_dir, "segments_min30.tsv"),
                             sep = "\t", stringsAsFactors = FALSE)
  draw_header("a", "Phylogenetic backbone and independent karyotype changes", 0.18, 11.85)

  ## AMK box
  box(6.75, 10.65, 2.5, 0.9, col = INK, lwd = 1.0)
  amk_row(7.05, 11.08, scale = 0.82)
  txt(8.0, 10.77, "AMK (x = 8)", 9.0, INK, face = "bold")
  txt(6.65, 11.10, "Medicago", 7.4, INK, face = "italic", just = c("right", "centre"))

  ## Root and collapsed but explicitly bifurcating clade backbone:
  ## Clade-III | (Clade-II + Clade-I).
  line(c(8.0, 8.0), c(10.65, 10.35), col = INK, lwd = 1.2)
  line(c(0.75, 8.0), c(10.35, 10.35), col = MUTED, lwd = 0.9, lty = 2)
  line(c(2.1, 8.0), c(9.55, 9.55), col = CL_COL["III"], lwd = 1.2)
  line(c(8.0, 9.2), c(9.55, 9.55), col = "#777777", lwd = 1.0)
  line(c(9.2, 9.2), c(9.55, 8.95), col = "#777777", lwd = 1.0)
  line(c(6.15, 10.8), c(8.95, 8.95), col = "#777777", lwd = 1.0)
  line(c(2.1, 2.1), c(9.55, 8.95), col = CL_COL["III"], lwd = 1.2)
  line(c(6.15, 6.15), c(8.95, 8.95), col = CL_COL["II"], lwd = 1.2)
  line(c(10.8, 10.8), c(8.95, 8.95), col = CL_COL["I"], lwd = 1.2)
  txt(2.1, 9.15, "Clade-III", 8.2, CL_COL["III"], face = "bold")
  txt(6.15, 9.15, "Clade-II", 8.2, CL_COL["II"], face = "bold")
  txt(10.8, 9.15, "Clade-I", 8.2, CL_COL["I"], face = "bold")

  ## Outgroup is retained for rooting but has no polarized event call.
  line(c(0.75, 0.75), c(10.35, 6.95), col = MUTED, lwd = 0.9, lty = 2)
  txt(0.75, 6.77, "Melilotus albus\n(outgroup; direction unresolved)", 5.2, MUTED, face = "italic")
  draw_tip(0.43, "Melilotus albus", "Malbus", "out", tip_segments, FALSE, FALSE)

  ## Species positions and labels.  Dashed stems are unchanged AMK-like states.
  tips <- data.frame(
    x = c(1.35, 1.85, 2.35, 4.85, 5.35, 5.85, 6.35,
          8.05, 8.65, 9.25, 9.85, 10.45, 11.05, 11.65, 12.25, 12.85, 13.45, 14.05),
    lab = c("Mlan", "Mrut", "Marc", "Mrad", "Medg", "Mfis", "Mlup", "Msuf", "Mcar", "Msec",
            "Morb", "Mpol", "Mpra", "Mtru_R108", "Mmar", "Mcre", "Msat_zm4", "Msat_cae"),
    name = c("Medicago lanigera", "Medicago ruthenica", "Medicago archiducis-nicolai",
             "Medicago radiata", "Medicago edgeworthii", "Medicago fisheriana", "Medicago lupulina",
             "Medicago suffruticosa", "Medicago carstiensis", "Medicago secundiflora", "Medicago orbicularis",
             "Medicago polymorpha", "Medicago praecox", "Medicago truncatula R108", "Medicago marina",
             "Medicago cretacea", "Medicago sativa cv. Zhongmu-4", "Medicago sativa subsp. caerulea"),
    clade = c(rep("III", 3), rep("II", 3), rep("I", 12)),
    x7 = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  for (i in seq_len(nrow(tips))) {
    changed <- tips$name %in% c("Medicago lanigera", "Medicago radiata", "Medicago edgeworthii",
                                "Medicago fisheriana", "Medicago suffruticosa", "Medicago secundiflora",
                                "Medicago polymorpha", "Medicago praecox")
    lty <- if (changed[i]) 1 else 2
    line(c(tips$x[i] + 0.30, tips$x[i] + 0.30), c(8.95, 6.95),
         col = CL_COL[tips$clade[i]], lwd = if (changed[i]) 1.0 else 0.75, lty = lty)
    draw_tip(tips$x[i], tips$name[i], tips$lab[i], tips$clade[i], tip_segments, changed[i], tips$x7[i])
  }

  ## Shared and independent event annotations.
  txt(1.95, 8.25, "shared Tr\nAMK7 -> AMK5", 6.1, GREEN, face = "bold")
  txt(5.85, 8.15, "independent\nRCT", 5.9, GREEN, face = "bold")
  txt(6.35, 7.72, ">=6 RCT-like\n+ fusion", 5.9, GREEN, face = "bold")
  txt(10.95, 8.15, "2xRCT + EEJ?", 5.9, GREEN, face = "bold")
  txt(11.65, 7.72, "RCT + NCF?", 5.9, GREEN, face = "bold")

  ## Legend
  box(12.05, 10.45, 3.55, 0.55, col = NA, lwd = 0, fill = NA)
  line(c(12.2, 12.55), c(10.80, 10.80), col = INK, lwd = 1.0)
  txt(12.68, 10.80, "event-bearing branch", 5.6, INK, just = c("left", "centre"))
  line(c(12.2, 12.55), c(10.59, 10.59), col = MUTED, lwd = 0.9, lty = 2)
  txt(12.68, 10.59, "same AMK-like state", 5.6, MUTED, just = c("left", "centre"))

}

draw_module <- function(x, y, w, h, title, col, type = "rct") {
  box(x, y, w, h, col = col, lwd = 1.1, fill = "white")
  txt(x + 0.12, y + h - 0.05, title, 7.0, INK, face = "italic", just = c("left", "top"))
  if (type == "tr") {
    pair(x + 0.25, y + 0.53, 5, 7, out = c(5, 7), op = "Tr", scale = 0.68, labels = FALSE)
    txt(x + 0.25, y + 0.28, "Chr5 contains a small AMK7 segment; x remains 8", 5.3, MUTED, just = c("left", "centre"))
    txt(x + w - 0.12, y + h - 0.05, "shared; supported", 5.1, GREEN, face = "bold", just = c("right", "top"))
  } else if (type == "rct") {
    pair(x + 0.25, y + 0.46, 1, 6, out = c(1, 6), op = "RCT", scale = 0.68, labels = TRUE)
    txt(x + 0.25, y + 0.24, "minimum model; x = 8", 5.3, MUTED, just = c("left", "centre"))
  } else if (type == "dual") {
    a1 <- if (grepl("suffruticosa", title)) 1 else 2
    b1 <- 6
    a2 <- if (grepl("suffruticosa", title)) 1 else 7
    b2 <- if (grepl("suffruticosa", title)) 8 else 8
    pair(x + 0.20, y + 0.78, a1, b1, out = c(a1, b1), op = "RCT", scale = 0.58, labels = FALSE)
    pair(x + 0.20, y + 0.26, a2, b2, out = c(a2, b2), op = "RCT", scale = 0.58, labels = FALSE)
    txt(x + w - 0.12, y + 0.22, "two independent RCTs", 5.0, MUTED, just = c("right", "centre"))
  } else if (type == "mosaic") {
    txt(x + 0.14, y + h - 0.16, ">=6 RCT-like + fusion; order unresolved", 5.0, GREEN, face = "bold", just = c("left", "top"))
    ## Representative AMK pieces on seven chromosomes; no false exact order.
    for (i in 1:7) {
      seg <- if (i %% 2 == 0) mix(c((i %% 8) + 1, ((i + 2) %% 8) + 1), c(.55, .45)) else mix(c((i %% 8) + 1, ((i + 3) %% 8) + 1), c(.65, .35))
      chrom(x + 0.25 + (i - 1) * 0.17, y + 0.28, w = 0.085, h = 0.55, seg, border = NA)
      txt(x + 0.292 + (i - 1) * 0.17, y + 0.15, paste0("Chr", i), 4.7, MUTED)
    }
    txt(x + w - 0.12, y + 0.26, "Mfis x = 7", 5.2, "#C0392B", face = "bold", just = c("right", "centre"))
  } else if (type == "x7") {
    txt(x + 0.14, y + h - 0.16, "minimum candidate pathway; order unresolved", 5.1, GREEN, face = "bold", just = c("left", "top"))
    if (type == "x7") {
      txt(x + 0.35, y + 0.62, "AMK", 5.7, MUTED)
      inputs <- if (grepl("polymorpha", title)) c(3, 5, 6) else c(5, 6, 8)
      ix <- x + 0.65
      for (k in seq_along(inputs)) {
        chrom(ix + (k - 1) * 0.17, y + 0.48, w = 0.095, h = 0.55,
              amk_seg(inputs[k]), border = NA)
        txt(ix + 0.047 + (k - 1) * 0.17, y + 0.32, as.character(inputs[k]), 4.7, MUTED)
        if (k < length(inputs)) plus(ix + 0.122 + (k - 1) * 0.17, y + 0.75, 6.5)
      }
      end_in <- ix + (length(inputs) - 1) * 0.17 + 0.095
      line(c(end_in + 0.14, end_in + 0.64), c(y + 0.76, y + 0.76), col = INK, lwd = 0.8, arrow_end = TRUE)
      txt(end_in + 0.39, y + 0.96, "RCT?", 5.2, GREEN, face = "bold")
      line(c(end_in + 0.89, end_in + 1.39), c(y + 0.76, y + 0.76), col = INK, lwd = 0.8, lty = 2, arrow_end = TRUE)
      txt(end_in + 1.14, y + 0.96, if (title == "Medicago praecox") "NCF?" else "EEJ?", 5.2, GREEN, face = "bold")
      out1 <- if (grepl("polymorpha", title)) c(5, 6, 3) else c(5, 8, 6)
      out2 <- if (grepl("polymorpha", title)) c(6, 3, 5) else c(6, 8, 5)
      for (i in 1:2) {
        seg <- if (i == 1) mix(out1, c(.35, .30, .35)) else mix(out2, c(.35, .32, .33))
        chrom(end_in + 1.66 + (i - 1) * 0.20, y + 0.48, w = 0.095, h = 0.55, seg, border = NA)
      }
      txt(end_in + 2.04, y + 0.26, "x = 7", 5.2, "#C0392B", face = "bold")
    }
  }
}

draw_pathways <- function() {
  draw_header("b", "Detailed minimum pathways for informative rearrangements", 0.18, 4.82)
  txt(15.82, 4.82, "solid: segment-supported; dashed: candidate/order unresolved", 5.4, MUTED, just = c("right", "top"))

  draw_module(0.25, 3.15, 3.72, 1.35, "M. ruthenica + M. archiducis-nicolai", CL_COL["III"], "tr")
  draw_module(4.12, 3.15, 3.72, 1.35, "M. lanigera", CL_COL["III"], "rct")
  draw_module(7.99, 3.15, 3.72, 1.35, "M. secundiflora", CL_COL["I"], "dual")
  draw_module(11.86, 3.15, 3.89, 1.35, "M. suffruticosa", CL_COL["I"], "dual")

  draw_module(0.25, 1.40, 4.90, 1.45, "M. polymorpha (x = 7)", CL_COL["I"], "x7")
  draw_module(5.28, 1.40, 4.90, 1.45, "M. praecox (x = 7)", CL_COL["I"], "x7")
  draw_module(10.31, 1.40, 5.44, 1.45, "M. fisheriana (genome_M46_2; x = 7)", CL_COL["II"], "mosaic")

  txt(0.25, 0.88,
      "Numbers label AMK chromosomes; Chr# labels the resulting species chromosome. Intermediate states are schematic and most parsimonious.\n",
      5.6, MUTED, just = c("left", "top"))
  txt(0.25, 0.57,
      "Midpoint-filled intervals are display-only; exact breakpoint positions and event order require breakpoint-level validation.",
      5.6, MUTED, just = c("left", "top"))
}

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 12),
                        gp = gpar(fontfamily = "Arial")))
  draw_backbone()
  draw_pathways()
  upViewport()
}

save_one <- function(kind, file) {
  if (kind == "pdf") cairo_pdf(file, width = 16, height = 12, family = "Arial")
  if (kind == "svg") svg(file, width = 16, height = 12, pointsize = 10, family = "Arial")
  if (kind == "png") png(file, width = 16, height = 12, units = "in", res = 600, type = "cairo", bg = "white")
  draw_all()
  dev.off()
}

save_one("pdf", file.path(out_dir, "karyotype_backbone_pathways_R.pdf"))
save_one("svg", file.path(out_dir, "karyotype_backbone_pathways_R.svg"))
save_one("png", file.path(out_dir, "karyotype_backbone_pathways_R.png"))

