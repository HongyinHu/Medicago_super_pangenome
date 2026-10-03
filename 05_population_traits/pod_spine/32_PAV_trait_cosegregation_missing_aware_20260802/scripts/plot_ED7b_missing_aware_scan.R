#!/usr/bin/env Rscript
# Extended Data Fig. 7b (revised): genome-wide phenotype-segregation scan of high-confidence PAVs,
# drawn from the missing-aware reanalysis (N_4.pod_spiny/32_PAV_trait_cosegregation_missing_aware_20260802):
# NA never converted to absence; >=12/14 called and >=5 called per group; two-sided Fisher exact test.
# Chr23997 is shown with its IGV-corrected call (filled diamond) and its raw call (open diamond).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_ED7b_missing_aware_scan.R <scan.tsv.gz> <chr_len.tsv> <out_prefix>")
}
d <- read.delim(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
f <- read.delim(args[2], sep = "\t", header = FALSE, col.names = c("chrom", "len"), stringsAsFactors = FALSE)
out_prefix <- args[3]

# Raw (uncorrected) Chr23997 call: 7/7 DEL vs 3/7 DEL, Fisher P = 0.06993
# (32_.../results/Chr23997_association_raw_vs_IGV_corrected.tsv).
raw_target_p <- 0.06993006993006994

chroms <- paste0("Chr", 1:8)
f <- f[match(chroms, f$chrom), ]
f$offset <- c(0, cumsum(head(as.numeric(f$len), -1)))
off <- setNames(f$offset, f$chrom)
chr_total <- sum(as.numeric(f$len)) / 1e6
chr_mids <- (f$offset + f$len / 2) / 1e6

is_true <- function(x) x %in% c(TRUE, "True", "TRUE", "true", 1, "1")
d <- d[d$CHROM %in% chroms, ]
d$x <- (off[d$CHROM] + as.numeric(d$POS)) / 1e6
d$is_target <- is_true(d$target_Chr23997)
d$is_complete <- is_true(d$complete_separation) & !d$is_target
d$is_group80 <- is_true(d$group_specific80) & !d$is_complete & !d$is_target
target <- d[d$is_target, ]
stopifnot(nrow(target) == 1)

# Small vertical jitter for background and 80% points only (P values are discrete);
# complete-separation points and Chr23997 are drawn at their exact -log10(P).
set.seed(23997)
d$y <- d$neglog10p
jit <- d$is_group80 | !(d$is_complete | d$is_target)
d$y[jit] <- pmax(0, d$y[jit] + runif(sum(jit), -0.065, 0.065))
ymax <- 4.5

# Rasterize background points (per-chromosome alternating colours).
bg <- d[!(d$is_group80 | d$is_complete | d$is_target), ]
px_w <- 4200L; px_h <- 1450L
img <- matrix("#FFFFFF", nrow = px_h, ncol = px_w)
px <- pmin(px_w, pmax(1L, as.integer(bg$x / chr_total * (px_w - 1L)) + 1L))
py <- px_h - pmin(px_h, pmax(1L, as.integer(bg$y / ymax * (px_h - 1L)) + 1L)) + 1L
bg_cols <- ifelse(match(bg$CHROM, chroms) %% 2 == 1, "#88CCF2", "#F05F9D")
for (k in seq_along(px)) {
  img[py[k], px[k]] <- bg_cols[k]
  if (px[k] < px_w) img[py[k], px[k] + 1L] <- bg_cols[k]
}
ras <- as.raster(img)

col_g80 <- "#2B2B2B"; col_comp <- "#E69F00"; col_tgt <- "#D55E00"
plot_one <- function(device, file, width = 8.8, height = 3.3, res = 600) {
  if (device == "pdf") {
    cairo_pdf(file, width = width, height = height, family = "Times New Roman")
  } else {
    png(file, width = width, height = height, units = "in", res = res, type = "cairo", family = "Times New Roman")
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(3.4, 4.2, 0.8, 0.8), mgp = c(2.1, 0.5, 0), tcl = -0.3, xaxs = "i", yaxs = "i")
  plot(NA, NA, xlim = c(0, chr_total), ylim = c(0, ymax), axes = FALSE, xlab = "", ylab = "")
  rasterImage(ras, 0, 0, chr_total, ymax, interpolate = FALSE)
  abline(h = -log10(0.05), lty = 2, col = "#8C8C8C", lwd = 0.8)
  abline(h = -log10(0.01), lty = 3, col = "#555555", lwd = 0.8)
  g80 <- d[d$is_group80, ]
  points(g80$x, g80$y, pch = 16, cex = 0.38, col = adjustcolor(col_g80, 0.9))
  comp <- d[d$is_complete, ]
  points(comp$x, comp$y, pch = 16, cex = 0.62, col = col_comp)
  # Chr23997: raw call (open) and IGV-corrected call (filled), joined by a thin dashed segment.
  segments(target$x, -log10(raw_target_p), target$x, target$neglog10p - 0.07, lty = 2, lwd = 0.6, col = "#555555")
  points(target$x, -log10(raw_target_p), pch = 23, bg = "white", col = "black", cex = 0.95, lwd = 0.8)
  points(target$x, target$neglog10p, pch = 23, bg = col_tgt, col = "black", cex = 1.25, lwd = 0.9)
  arrows(target$x, target$neglog10p + 0.42, target$x, target$neglog10p + 0.1, length = 0.05, lwd = 0.8)
  text(target$x, target$neglog10p + 0.42, "Chr23997", pos = 3, cex = 0.85, offset = 0.25)
  # Threshold labels on a white box so they stay legible over points.
  for (p in c(0.05, 0.01)) {
    lab <- sprintf("P = %s", p)
    w <- strwidth(lab, cex = 0.72, font = 3); h <- strheight(lab, cex = 0.72, font = 3)
    xr <- chr_total * 0.997; yt <- -log10(p) - 0.04
    rect(xr - w - 0.4, yt - h * 1.5, xr + 0.2, yt, col = "white", border = NA)
    text(xr, yt - h * 0.75, lab, adj = c(1, 0.5), col = "#D62728", font = 3, cex = 0.72)
  }
  axis(2, at = 0:4, las = 1, lwd = 0.8, cex.axis = 0.9)
  axis(1, at = chr_mids, labels = 1:8, lwd = 0.8, cex.axis = 0.9)
  box(bty = "l", lwd = 0.8)
  title(xlab = "Chromosome", line = 2.0, cex.lab = 1.05)
  title(ylab = expression(-log[10](italic(P))), line = 2.2, cex.lab = 1.05)
  legend(x = chr_total, y = ymax, xjust = 1, yjust = 1, bty = "n", cex = 0.7, y.intersp = 1.1, ncol = 2,
         x.intersp = 0.8, text.width = c(chr_total * 0.235, chr_total * 0.165),
         legend = c("80% group-specific", "Complete phenotype-group separation",
                    "Chr23997 (IGV-corrected)", "Chr23997 (raw call)"),
         pch = c(16, 16, 23, 23), col = c(col_g80, col_comp, "black", "black"),
         pt.bg = c(NA, NA, col_tgt, "white"), pt.cex = c(0.7, 0.9, 1.1, 0.9))
}
plot_one("pdf", paste0(out_prefix, ".pdf"))
plot_one("png", paste0(out_prefix, ".png"))

summ <- data.frame(
  metric = c("tested_sites", "group_specific80_sites_incl_complete", "complete_separation_sites_incl_Chr23997",
             "sites_at_minimum_P", "minimum_P", "Chr23997_corrected_P", "Chr23997_raw_P"),
  value = c(nrow(d), sum(is_true(d$group_specific80)), sum(is_true(d$complete_separation)),
            sum(d$fisher_p <= min(d$fisher_p) * (1 + 1e-9)), min(d$fisher_p), target$fisher_p, raw_target_p))
write.table(summ, paste0(out_prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
