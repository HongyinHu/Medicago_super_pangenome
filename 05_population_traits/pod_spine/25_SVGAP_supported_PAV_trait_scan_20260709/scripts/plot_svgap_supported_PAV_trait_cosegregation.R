args <- commandArgs(trailingOnly=TRUE)
out <- ifelse(length(args) >= 1, args[1], "path/to/project/N_4.pod_spiny/25_SVGAP_supported_PAV_trait_scan_20260709")
res <- file.path(out, "results/svgap_supported_PAV_trait_cosegregation.tsv")
tgt <- file.path(out, "results/Chr23997_manual_IGV_corrected_target_point.tsv")
fai <- "path/to/project/N_3.call_SV/01.read_based_dualref_hifi/02_ref_prepare/Msa/ref.fa.fai"
d <- read.delim(res, sep="\t", header=TRUE, stringsAsFactors=FALSE)
manual <- read.delim(tgt, sep="\t", header=TRUE, stringsAsFactors=FALSE)
chroms <- paste0("Chr", 1:8)
f <- read.delim(fai, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(f)[1:2] <- c("chrom", "len")
f <- f[f$chrom %in% chroms, c("chrom", "len")]
f$chrom <- factor(f$chrom, levels=chroms)
f <- f[order(f$chrom),]
f$offset <- c(0, cumsum(head(f$len, -1)))
f$mid <- f$offset + f$len / 2
off <- setNames(f$offset, f$chrom)
d <- d[d$CHROM %in% chroms,]
d$genome_x <- (off[d$CHROM] + as.numeric(d$POS)) / 1e6
d$category <- "SVgap-supported PAVs"
d$category[d$group_specific80 == "True" | d$group_specific80 == TRUE] <- "80% group-specific"
d$category[d$complete_segregation == "True" | d$complete_segregation == TRUE] <- "Complete segregation"
set.seed(23997)
# Visual-only vertical jitter: keep exact P values in the source table unchanged.
# Larger jitter makes the discrete Fisher-test P-value rows easier to read.
j_amp <- ifelse(d$category == "SVgap-supported PAVs", 0.090, 0.065)
d$y_plot <- pmax(0, d$neglog10p + runif(nrow(d), -j_amp, j_amp))
manual_x <- (off[manual$chrom[1]] + manual$pos[1]) / 1e6
manual_y <- manual$neglog10p_core[1]
cols_chr <- rep(c("#72BDEB", "#EE4B8B"), length.out=length(chroms))
cols <- c("SVgap-supported PAVs"="#72BDEB", "80% group-specific"="#2B2B2B", "Complete segregation"="#E69F00")
figbase <- file.path(out, "figures/svgap_supported_PAV_Chr1_Chr8_trait_cosegregation.neglog10P.jittered")
plot_one <- function(device_fun, file, width, height) {
  device_fun(file, width=width, height=height)
  par(mar=c(4.4,4.6,1.1,1.2), xaxs="i", yaxs="i")
  ymax <- max(4.0, max(d$neglog10p, manual_y, na.rm=TRUE) + 0.50)
  plot(d$genome_x, d$y_plot, type="n", xlab="Chromosome", ylab=expression(-log[10](italic(P))),
       xlim=c(0, sum(f$len)/1e6), ylim=c(0, ymax), axes=FALSE, cex.lab=1.05)
  axis(2, las=1, lwd=0.8, lwd.ticks=0.8)
  axis(1, at=f$mid/1e6, labels=gsub("Chr", "", as.character(f$chrom)), lwd=0.8, lwd.ticks=0.8)
  abline(h=-log10(0.05), lty=2, col="#8C8C8C", lwd=0.9)
  abline(h=-log10(0.01), lty=3, col="#555555", lwd=0.9)
  brks <- cumsum(f$len)/1e6
  abline(v=brks[-length(brks)], col="#D0D0D0", lwd=0.5)
  for (i in seq_along(chroms)) {
    idx <- which(d$CHROM == chroms[i] & d$category == "SVgap-supported PAVs")
    if (length(idx)) points(d$genome_x[idx], d$y_plot[idx], pch=20, cex=0.34, col=adjustcolor(cols_chr[i], alpha.f=0.58))
  }
  for (cat in c("80% group-specific", "Complete segregation")) {
    idx <- which(d$category == cat)
    if (length(idx)) points(d$genome_x[idx], d$y_plot[idx], pch=20, cex=0.68, col=adjustcolor(cols[cat], alpha.f=0.92))
  }
  raw_tgt <- which(d$target_like_Chr23997 == "True" | d$target_like_Chr23997 == TRUE)
  if (length(raw_tgt)) points(d$genome_x[raw_tgt], d$y_plot[raw_tgt], pch=23, bg="white", col="#2B2B2B", cex=0.92, lwd=0.9)
  points(manual_x, manual_y, pch=23, bg="#D55E00", col="black", cex=1.10, lwd=0.9)
  text(manual_x, min(ymax - 0.12, manual_y + 0.34), labels="Chr23997", font=3, cex=0.78, pos=3)
  arrows(manual_x, min(ymax - 0.24, manual_y + 0.28), manual_x, manual_y + 0.05, length=0.06, lwd=0.8)
  legend("topright", bty="n", cex=0.70,
         legend=c("SVgap-supported PAVs", "80% group-specific", "Complete segregation", "Chr23997 manual IGV-corrected", "Target-like high-confidence PAV", "P = 0.05", "P = 0.01", "Y jitter for display only"),
         col=c("#72BDEB", cols[2], cols[3], "black", "#2B2B2B", "#8C8C8C", "#555555", "#666666"),
         pt.bg=c("#72BDEB", cols[2], cols[3], "#D55E00", "white", NA, NA, NA),
         pch=c(20,20,20,23,23,NA,NA,NA), lty=c(NA,NA,NA,NA,NA,2,3,1),
         pt.cex=c(0.75,0.85,0.85,1.0,0.9,NA,NA,NA))
  box(bty="l", lwd=0.8)
  dev.off()
}
plot_one(function(file,width,height) grDevices::pdf(file, width=width, height=height, useDingbats=FALSE), paste0(figbase, ".pdf"), 8.8, 3.4)
plot_one(function(file,width,height) grDevices::png(file, width=width, height=height, units="in", res=600), paste0(figbase, ".png"), 8.8, 3.4)
