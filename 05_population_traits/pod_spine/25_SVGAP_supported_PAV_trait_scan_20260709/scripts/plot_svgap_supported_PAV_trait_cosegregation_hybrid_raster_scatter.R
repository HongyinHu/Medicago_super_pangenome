args <- commandArgs(trailingOnly=TRUE)
out <- ifelse(length(args) >= 1, args[1], "path/to/project/N_4.pod_spiny/25_SVGAP_supported_PAV_trait_scan_20260709")
res <- file.path(out, "results/svgap_supported_PAV_trait_cosegregation.tsv")
tgt <- file.path(out, "results/Chr23997_manual_IGV_corrected_target_point.tsv")
fai <- "path/to/project/N_3.call_SV/00_data/3.two_ref/genome_Msa.fa.fai"
figdir <- file.path(out, "figures")
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

chroms <- paste0("Chr", 1:8)
d <- read.delim(res, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
manual <- read.delim(tgt, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
f <- read.delim(fai, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(f)[1:2] <- c("chrom", "len")
f <- f[f$chrom %in% chroms, c("chrom", "len")]
f$chrom <- factor(f$chrom, levels=chroms)
f <- f[order(f$chrom),]
f$offset <- c(0, cumsum(head(f$len, -1)))
f$mid <- f$offset + f$len/2
off <- setNames(f$offset, as.character(f$chrom))

d <- d[d$CHROM %in% chroms,]
d$genome_x <- (off[d$CHROM] + as.numeric(d$POS)) / 1e6
d$is_group80 <- d$group_specific80 %in% c(TRUE, "True", "TRUE", "true", 1, "1")
d$is_complete <- d$complete_segregation %in% c(TRUE, "True", "TRUE", "true", 1, "1")
d$is_target_like <- d$target_like_Chr23997 %in% c(TRUE, "True", "TRUE", "true", 1, "1")
manual_x <- (off[manual$chrom[1]] + manual$pos[1]) / 1e6
manual_y <- manual$neglog10p_core[1]
set.seed(23997)
j_amp <- ifelse(d$is_group80 | d$is_complete, 0.065, 0.090)
d$y_plot <- pmax(0, d$neglog10p + runif(nrow(d), -j_amp, j_amp))
chr_total <- sum(f$len)/1e6
chr_breaks <- cumsum(f$len)/1e6
chr_mids <- f$mid/1e6
ymax <- max(4.0, max(d$neglog10p, manual_y, na.rm=TRUE) + 0.55)

# Build a raster layer for ordinary background points only. Highlighted points remain vector.
bg <- d[!(d$is_group80 | d$is_complete | d$is_target_like),]
px_w <- 4200L
px_h <- 1450L
img <- matrix("#FFFFFFFF", nrow=px_h, ncol=px_w)
cols_chr <- rep(c("#82C7F0", "#EF5C9A"), length.out=length(chroms))
# Convert to panel pixels. y=0 is bottom, matrix row 1 is top.
px <- pmin(px_w, pmax(1L, as.integer(bg$genome_x / chr_total * (px_w - 1L)) + 1L))
py_from_bottom <- pmin(px_h, pmax(1L, as.integer(bg$y_plot / ymax * (px_h - 1L)) + 1L))
py <- px_h - py_from_bottom + 1L
chr_idx <- match(bg$CHROM, chroms)
# Light alpha-like colors pre-blended on white, with larger effective point by 2 pixels.
point_cols <- ifelse(chr_idx %% 2 == 1, "#88CCF2", "#F05F9D")
for (k in seq_along(px)) {
  img[py[k], px[k]] <- point_cols[k]
  if (px[k] < px_w) img[py[k], px[k]+1L] <- point_cols[k]
}
ras <- as.raster(img)
rm(img); gc()

figbase <- file.path(figdir, "svgap_supported_PAV_Chr1_Chr8_trait_cosegregation.neglog10P.hybrid_raster_scatter")
plot_hybrid <- function(file, type=c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") {
    grDevices::pdf(file, width=8.8, height=3.4, useDingbats=FALSE)
  } else {
    grDevices::png(file, width=8.8, height=3.4, units="in", res=600)
  }
  par(mar=c(4.4,4.6,1.1,1.2), xaxs="i", yaxs="i")
  plot(NA, NA, xlim=c(0, chr_total), ylim=c(0, ymax), axes=FALSE,
       xlab="Chromosome", ylab=expression(-log[10](italic(P))))
  rasterImage(ras, xleft=0, ybottom=0, xright=chr_total, ytop=ymax, interpolate=FALSE)
  abline(h=-log10(0.05), lty=2, col="#8C8C8C", lwd=0.85)
  abline(h=-log10(0.01), lty=3, col="#555555", lwd=0.85)
  abline(v=chr_breaks[-length(chr_breaks)], col="#D5D5D5", lwd=0.5)
  g80 <- d[d$is_group80 & !d$is_complete,]
  if (nrow(g80) > 0) points(g80$genome_x, g80$y_plot, pch=20, cex=0.42, col=adjustcolor("#2B2B2B", alpha.f=0.90))
  comp <- d[d$is_complete,]
  if (nrow(comp) > 0) points(comp$genome_x, comp$y_plot, pch=20, cex=0.58, col=adjustcolor("#E69F00", alpha.f=0.96))
  target_like <- d[d$is_target_like,]
  if (nrow(target_like) > 0) points(target_like$genome_x, target_like$y_plot, pch=23, bg="white", col="#2B2B2B", cex=0.95, lwd=0.85)
  points(manual_x, manual_y, pch=23, bg="#D55E00", col="black", cex=1.12, lwd=0.9)
  text(manual_x, min(ymax - 0.12, manual_y + 0.34), labels="Chr23997", font=3, cex=0.78, pos=3)
  arrows(manual_x, min(ymax - 0.24, manual_y + 0.28), manual_x, manual_y + 0.05, length=0.06, lwd=0.8)
  axis(2, las=1, lwd=0.8, lwd.ticks=0.8)
  axis(1, at=chr_mids, labels=gsub("Chr", "", chroms), lwd=0.8, lwd.ticks=0.8)
  legend("topright", bty="n", cex=0.68,
         legend=c("All high-confidence PAVs (rasterized)", "80% group-specific", "Complete segregation", "Chr23997 manual IGV-corrected", "Target-like high-confidence PAV", "P = 0.05", "P = 0.01"),
         col=c("#82C7F0", "#2B2B2B", "#E69F00", "black", "#2B2B2B", "#8C8C8C", "#555555"),
         pt.bg=c("#82C7F0", NA, NA, "#D55E00", "white", NA, NA),
         pch=c(20,20,20,23,23,NA,NA), lty=c(NA,NA,NA,NA,NA,2,3),
         pt.cex=c(0.65,0.78,0.82,1.0,0.90,NA,NA))
  box(bty="l", lwd=0.8)
  dev.off()
}
plot_hybrid(paste0(figbase, ".pdf"), "pdf")
plot_hybrid(paste0(figbase, ".png"), "png")
cat("output_pdf\t", paste0(figbase, ".pdf"), "\n", sep="")
cat("output_png\t", paste0(figbase, ".png"), "\n", sep="")
cat("rasterized_background_points\t", nrow(bg), "\n", sep="")
cat("vector_group80_points\t", sum(d$is_group80 & !d$is_complete), "\n", sep="")
cat("vector_complete_points\t", sum(d$is_complete), "\n", sep="")
cat("vector_target_like_points\t", sum(d$is_target_like), "\n", sep="")
cat("manual_chr23997_neglog10p\t", manual_y, "\n", sep="")
