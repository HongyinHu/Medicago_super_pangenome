#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: plot_joint_gwas.R combined.tsv.gz reference.fa.fai outdir")

input <- args[1]
fai_path <- args[2]
outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(gzfile(input), stringsAsFactors = FALSE, check.names = FALSE)
fai <- read.delim(fai_path, header = FALSE, stringsAsFactors = FALSE)
colnames(fai)[1:2] <- c("chr", "length")
fai <- fai[fai$chr %in% paste0("Chr", 1:8), c("chr", "length")]
fai$chr <- factor(fai$chr, levels = paste0("Chr", 1:8))
fai <- fai[order(fai$chr), ]

dat <- dat[dat$chr %in% as.character(fai$chr) & is.finite(dat$minus_log10_p), ]
offset <- c(0, cumsum(fai$length[-nrow(fai)]))
names(offset) <- as.character(fai$chr)
dat$cum_pos <- dat$pos + offset[dat$chr]
centers <- offset + fai$length / 2
total_length <- sum(fai$length)
bonf <- -log10(0.05 / nrow(dat))
suggestive <- -log10(1 / nrow(dat))
target_x <- 88980942 + offset["Chr4"]

plot_manhattan <- function(filename, device = c("png", "tiff")) {
  device <- match.arg(device)
  if (device == "png") {
    png(filename, width = 3900, height = 1650, res = 300, type = "cairo")
  } else {
    tiff(filename, width = 3900, height = 1650, res = 300, compression = "lzw")
  }
  par(mar = c(4.8, 5.0, 1.2, 1.2), mgp = c(2.8, 0.8, 0), las = 1, family = "sans")
  ymax <- max(dat$minus_log10_p, bonf, na.rm = TRUE) * 1.10
  plot(NA, xlim = c(0, total_length), ylim = c(0, ymax), xaxt = "n",
       xlab = "Chromosome", ylab = expression(-log[10](italic(P))), bty = "l")
  abline(v = cumsum(fai$length)[-nrow(fai)], col = "#E5E5E5", lwd = 0.7)
  snp <- dat[dat$class == "SNP", ]
  indel <- dat[dat$class == "INDEL", ]
  sv <- dat[dat$class == "SV", ]
  points(snp$cum_pos, snp$minus_log10_p, pch = 20, cex = 0.18, col = "#AEB6BF80")
  points(indel$cum_pos, indel$minus_log10_p, pch = 20, cex = 0.45, col = "#E76F5190")
  points(sv$cum_pos, sv$minus_log10_p, pch = 21, cex = 0.65, bg = "#2A9D8F", col = "#1B6F66")
  abline(h = suggestive, col = "#777777", lty = 3, lwd = 0.9)
  abline(h = bonf, col = "#B2182B", lty = 2, lwd = 1.0)
  axis(1, at = centers, labels = 1:8, las = 1)
  arrows(target_x, ymax * 0.94, target_x, ymax * 0.82, length = 0.08, lwd = 0.9)
  text(target_x, ymax * 0.965, "Chr23997", cex = 0.78)
  legend("topright", legend = c("SNP", "INDEL", "SV", "Bonferroni", "1/N"),
         pch = c(20, 20, 21, NA, NA), pt.bg = c(NA, NA, "#2A9D8F", NA, NA),
         col = c("#7F8C8D", "#E76F51", "#1B6F66", "#B2182B", "#777777"),
         lty = c(NA, NA, NA, 2, 3), bty = "n", cex = 0.78)
  dev.off()
}

plot_manhattan(file.path(outdir, "SNP_INDEL_SV.manhattan.png"), "png")
plot_manhattan(file.path(outdir, "SNP_INDEL_SV.manhattan.tiff"), "tiff")

png(file.path(outdir, "SNP_INDEL_SV.QQ.png"), width = 3000, height = 1000, res = 300, type = "cairo")
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.2, 1.0), mgp = c(2.6, 0.8, 0), family = "sans")
cols <- c(SNP = "#6C5B7B", INDEL = "#E76F51", SV = "#2A9D8F")
for (klass in c("SNP", "INDEL", "SV")) {
  p <- sort(dat$p_wald[dat$class == klass & dat$p_wald > 0 & dat$p_wald <= 1])
  n <- length(p)
  keep <- unique(round(seq(1, n, length.out = min(n, 100000))))
  obs <- -log10(p[keep])
  expv <- -log10((keep - 0.5) / n)
  plot(expv, obs, pch = 20, cex = 0.25, col = paste0(cols[klass], "80"),
       xlab = expression(Expected~~-log[10](italic(P))),
       ylab = expression(Observed~~-log[10](italic(P))), main = klass, bty = "l")
  abline(0, 1, lty = 2, col = "#555555")
}
dev.off()

local <- dat[dat$chr == "Chr4" & dat$pos >= 87980831 & dat$pos <= 89981052, ]
png(file.path(outdir, "Chr23997.plusminus1Mb.png"), width = 2600, height = 1400, res = 300, type = "cairo")
par(mar = c(4.8, 5.0, 1.2, 1.0), mgp = c(2.8, 0.8, 0), family = "sans")
plot(NA, xlim = c(87.980831, 89.981052),
     ylim = c(0, max(local$minus_log10_p, na.rm = TRUE) * 1.12),
     xlab = "Chr4 position (Mb)", ylab = expression(-log[10](italic(P))), bty = "l")
rect(88.980831, par("usr")[3], 88.981052, par("usr")[4], col = "#F4A26135", border = NA)
for (klass in c("SNP", "INDEL", "SV")) {
  d <- local[local$class == klass, ]
  pch <- if (klass == "SV") 21 else 20
  points(d$pos / 1e6, d$minus_log10_p, pch = pch, cex = if (klass == "SNP") 0.35 else 0.75,
         col = cols[klass], bg = if (klass == "SV") cols[klass] else NA)
}
abline(v = c(88.980831, 88.981052), col = "#D55E00", lty = 2, lwd = 0.8)
legend("topright", legend = c("SNP", "INDEL", "SV", "validated 221-bp DEL"),
       pch = c(20, 20, 21, NA), lty = c(NA, NA, NA, 2),
       col = c(cols["SNP"], cols["INDEL"], cols["SV"], "#D55E00"),
       pt.bg = c(NA, NA, cols["SV"], NA), bty = "n", cex = 0.8)
dev.off()

cat("Plots written to", outdir, "\n")
