args <- commandArgs(trailingOnly=TRUE)
out <- ifelse(length(args) >= 1, args[1], "path/to/project/N_4.pod_spiny/22_readmapping_rawSV_trait_cosegregation_scan_20260708")
res <- file.path(out, "results/readmapping_rawSV_Chr4_86_92Mb.strict_clusters_with_fisher.tsv")
tgt <- file.path(out, "results/Chr23997_manual_IGV_corrected_target_point.tsv")
d <- read.delim(res, sep="\t", header=TRUE, stringsAsFactors=FALSE)
manual <- read.delim(tgt, sep="\t", header=TRUE, stringsAsFactors=FALSE)
if (nrow(d) == 0) stop("No strict clusters available")
d$pos_mb <- d$pos / 1e6
d$category <- "All raw read-mapping SV clusters"
d$category[d$group_specific80 == "True" | d$group_specific80 == TRUE] <- "80% group-specific"
d$category[d$complete_segregation == "True" | d$complete_segregation == TRUE] <- "Complete segregation"
# target-like raw cluster, if any, gets its own symbol layer later
manual_x <- manual$pos[1] / 1e6
manual_y <- manual$neglog10p_core[1]
figbase <- file.path(out, "figures/readmapping_rawSV_Chr4_86_92Mb_cosegregation.neglog10P")
cols <- c("All raw read-mapping SV clusters"="#8DD3F7", "80% group-specific"="#2B2B2B", "Complete segregation"="#E69F00")
plot_one <- function(device_fun, file, width, height) {
  device_fun(file, width=width, height=height)
  par(mar=c(4.2,4.4,1.0,1.0), xaxs="i", yaxs="i")
  y <- pmin(d$neglog10p, max(6, max(d$neglog10p, manual_y, na.rm=TRUE)))
  plot(d$pos_mb, y, type="n", xlab="Chr4 position (Mb)", ylab=expression(-log[10](italic(P))),
       xlim=c(86,92), ylim=c(0, max(2.6, max(y, manual_y, na.rm=TRUE)+0.3)), cex.lab=1.05)
  abline(h=-log10(0.05), lty=2, col="#9E9E9E", lwd=1)
  # plot from low-priority to high-priority
  for (cat in names(cols)) {
    idx <- which(d$category == cat)
    if (length(idx)) {
      pch <- ifelse(cat == "All raw read-mapping SV clusters", 20, 20)
      cex <- ifelse(cat == "All raw read-mapping SV clusters", 0.45, 0.65)
      points(d$pos_mb[idx], d$neglog10p[idx], pch=pch, cex=cex, col=cols[cat])
    }
  }
  raw_tgt <- which(d$target_like_raw_cluster == "True" | d$target_like_raw_cluster == TRUE)
  if (length(raw_tgt)) {
    points(d$pos_mb[raw_tgt], d$neglog10p[raw_tgt], pch=23, bg="white", col="#2B2B2B", cex=1.0, lwd=1.0)
  }
  points(manual_x, manual_y, pch=23, bg="#D55E00", col="black", cex=1.15, lwd=0.9)
  text(manual_x, manual_y+0.28, labels="Chr23997", font=3, cex=0.8)
  arrows(manual_x, manual_y+0.22, manual_x, manual_y+0.04, length=0.07, lwd=0.8)
  legend("topright", bty="n", cex=0.72,
         legend=c("All raw read-mapping SV clusters", "80% group-specific", "Complete segregation", "Chr23997 manual IGV-corrected", "Target-like raw cluster"),
         col=c(cols[1], cols[2], cols[3], "black", "#2B2B2B"),
         pt.bg=c(cols[1], cols[2], cols[3], "#D55E00", "white"),
         pch=c(20,20,20,23,23), pt.cex=c(0.7,0.8,0.8,1.1,1.0))
  dev.off()
}
# PDF/SVG/PNG. Use base R devices only.
plot_one(function(file,width,height) grDevices::pdf(file, width=width, height=height, useDingbats=FALSE), paste0(figbase, ".pdf"), 7.2, 3.0)
plot_one(function(file,width,height) grDevices::svg(file, width=width, height=height), paste0(figbase, ".svg"), 7.2, 3.0)
plot_one(function(file,width,height) grDevices::png(file, width=width, height=height, units="in", res=600), paste0(figbase, ".png"), 7.2, 3.0)
# Raw P-value version, where lower is stronger; useful because the user previously asked for p-value axis.
d$rawp_plot <- pmax(d$fisher_p, 1e-6)
figbase2 <- file.path(out, "figures/readmapping_rawSV_Chr4_86_92Mb_cosegregation.rawP")
plot_rawp <- function(device_fun, file, width, height) {
  device_fun(file, width=width, height=height)
  par(mar=c(4.2,4.4,1.0,1.0), xaxs="i", yaxs="i")
  plot(d$pos_mb, d$fisher_p, type="n", xlab="Chr4 position (Mb)", ylab=expression(italic(P)), xlim=c(86,92), ylim=c(1,0), cex.lab=1.05)
  abline(h=0.05, lty=2, col="#9E9E9E", lwd=1)
  for (cat in names(cols)) {
    idx <- which(d$category == cat)
    if (length(idx)) points(d$pos_mb[idx], d$fisher_p[idx], pch=20, cex=ifelse(cat==names(cols)[1],0.45,0.65), col=cols[cat])
  }
  raw_tgt <- which(d$target_like_raw_cluster == "True" | d$target_like_raw_cluster == TRUE)
  if (length(raw_tgt)) points(d$pos_mb[raw_tgt], d$fisher_p[raw_tgt], pch=23, bg="white", col="#2B2B2B", cex=1.0, lwd=1.0)
  points(manual_x, manual$fisher_p_core[1], pch=23, bg="#D55E00", col="black", cex=1.15, lwd=0.9)
  text(manual_x, manual$fisher_p_core[1]+0.08, labels="Chr23997", font=3, cex=0.8)
  legend("bottomright", bty="n", cex=0.72,
         legend=c("All raw read-mapping SV clusters", "80% group-specific", "Complete segregation", "Chr23997 manual IGV-corrected", "Target-like raw cluster"),
         col=c(cols[1], cols[2], cols[3], "black", "#2B2B2B"), pt.bg=c(cols[1], cols[2], cols[3], "#D55E00", "white"), pch=c(20,20,20,23,23), pt.cex=c(0.7,0.8,0.8,1.1,1.0))
  dev.off()
}
plot_rawp(function(file,width,height) grDevices::pdf(file, width=width, height=height, useDingbats=FALSE), paste0(figbase2, ".pdf"), 7.2, 3.0)
plot_rawp(function(file,width,height) grDevices::svg(file, width=width, height=height), paste0(figbase2, ".svg"), 7.2, 3.0)
plot_rawp(function(file,width,height) grDevices::png(file, width=width, height=height, units="in", res=600), paste0(figbase2, ".png"), 7.2, 3.0)
