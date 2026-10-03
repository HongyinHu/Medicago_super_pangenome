args <- commandArgs(trailingOnly=TRUE)
out <- ifelse(length(args) >= 1, args[1], "path/to/project/N_4.pod_spiny/22_readmapping_rawSV_trait_cosegregation_scan_20260708")
res <- file.path(out, "results/readmapping_rawSV_genomewide_fast.breakpoint_bin_clusters_with_fisher.tsv")
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
d <- d[d$chrom %in% chroms,]
d$genome_x <- (off[d$chrom] + d$pos) / 1e6
d$category <- "All raw read-mapping SV clusters"
d$category[d$group_specific80 == "True" | d$group_specific80 == TRUE] <- "80% group-specific"
d$category[d$complete_segregation == "True" | d$complete_segregation == TRUE] <- "Complete segregation"
set.seed(23997)
# Visual-only vertical jitter to avoid artificial-looking rows from discrete Fisher exact-test p-values.
# The source table keeps the exact, unjittered fisher_p and neglog10p values.
j_amp <- ifelse(d$category == "All raw read-mapping SV clusters", 0.035, 0.025)
d$y_plot <- pmax(0, d$neglog10p + runif(nrow(d), -j_amp, j_amp))
manual_x <- (off[manual$chrom[1]] + manual$pos[1]) / 1e6
manual_y <- manual$neglog10p_core[1]
cols_chr <- rep(c("#8DD3F7", "#EF4A8A"), length.out=length(chroms))
cols <- c("All raw read-mapping SV clusters"="#8DD3F7", "80% group-specific"="#2B2B2B", "Complete segregation"="#E69F00")
figbase <- file.path(out, "figures/readmapping_rawSV_genomewide_fast_Chr1_Chr8_cosegregation.neglog10P.jittered")
plot_one <- function(device_fun, file, width, height) {
  device_fun(file, width=width, height=height)
  par(mar=c(4.4,4.6,1.1,1.2), xaxs="i", yaxs="i")
  ymax <- max(4.2, max(d$neglog10p, manual_y, na.rm=TRUE) + 0.45)
  plot(d$genome_x, d$y_plot, type="n", xlab="Chromosome", ylab=expression(-log[10](italic(P))),
       xlim=c(0, sum(f$len)/1e6), ylim=c(0, ymax), axes=FALSE, cex.lab=1.05)
  axis(2, las=1, lwd=0.8, lwd.ticks=0.8)
  axis(1, at=f$mid/1e6, labels=gsub("Chr", "", as.character(f$chrom)), lwd=0.8, lwd.ticks=0.8)
  abline(h=-log10(0.05), lty=2, col="#9E9E9E", lwd=0.9)
  brks <- cumsum(f$len)/1e6
  abline(v=brks[-length(brks)], col="#D0D0D0", lwd=0.5)
  for (i in seq_along(chroms)) {
    idx <- which(d$chrom == chroms[i] & d$category == "All raw read-mapping SV clusters")
    if (length(idx)) points(d$genome_x[idx], d$y_plot[idx], pch=20, cex=0.24, col=adjustcolor(cols_chr[i], alpha.f=0.45))
  }
  for (cat in c("80% group-specific", "Complete segregation")) {
    idx <- which(d$category == cat)
    if (length(idx)) points(d$genome_x[idx], d$y_plot[idx], pch=20, cex=0.55, col=adjustcolor(cols[cat], alpha.f=0.9))
  }
  raw_tgt <- which(d$target_like_raw_cluster == "True" | d$target_like_raw_cluster == TRUE)
  if (length(raw_tgt)) points(d$genome_x[raw_tgt], d$y_plot[raw_tgt], pch=23, bg="white", col="#2B2B2B", cex=0.85, lwd=0.9)
  points(manual_x, manual_y, pch=23, bg="#D55E00", col="black", cex=1.05, lwd=0.9)
  text(manual_x, min(ymax - 0.12, manual_y + 0.34), labels="Chr23997", font=3, cex=0.78, pos=3)
  arrows(manual_x, min(ymax - 0.24, manual_y + 0.28), manual_x, manual_y + 0.05, length=0.06, lwd=0.8)
  legend("topright", bty="n", cex=0.72,
         legend=c("All raw read-mapping SV clusters", "80% group-specific", "Complete segregation", "Chr23997 manual IGV-corrected", "Target-like raw cluster", "Y jitter for display only"),
         col=c("#8DD3F7", cols[2], cols[3], "black", "#2B2B2B", "#666666"),
         pt.bg=c("#8DD3F7", cols[2], cols[3], "#D55E00", "white", NA),
         pch=c(20,20,20,23,23,NA), lty=c(NA,NA,NA,NA,NA,1), pt.cex=c(0.65,0.75,0.75,1.0,0.9,NA))
  box(bty="l", lwd=0.8)
  dev.off()
}
plot_one(function(file,width,height) grDevices::pdf(file, width=width, height=height, useDingbats=FALSE), paste0(figbase, ".pdf"), 8.8, 3.4)
plot_one(function(file,width,height) grDevices::png(file, width=width, height=height, units="in", res=600), paste0(figbase, ".png"), 8.8, 3.4)
