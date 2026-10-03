args <- commandArgs(trailingOnly=TRUE)
out <- ifelse(length(args) >= 1, args[1], "path/to/project/N_4.pod_spiny/25_SVGAP_supported_PAV_trait_scan_20260709")
res <- file.path(out, "results/svgap_supported_PAV_trait_cosegregation.tsv")
tgt <- file.path(out, "results/Chr23997_manual_IGV_corrected_target_point.tsv")
fai <- "path/to/project/N_3.call_SV/00_data/3.two_ref/genome_Msa.fa.fai"
figdir <- file.path(out, "figures")
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

bin_size <- 100000L
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

# Keep only target chromosomes and define visual classes.
d <- d[d$CHROM %in% chroms,]
d$bin <- floor((as.numeric(d$POS) - 1) / bin_size)
d$bin_start <- d$bin * bin_size + 1
d$bin_end <- pmin((d$bin + 1) * bin_size, f$len[match(d$CHROM, f$chrom)])
d$x_mid <- (off[d$CHROM] + (d$bin_start + d$bin_end)/2) / 1e6
d$chr_idx <- match(d$CHROM, chroms)
d$base_class <- ifelse(d$chr_idx %% 2 == 1, "odd_chr", "even_chr")
d$is_group80 <- d$group_specific80 %in% c(TRUE, "True", "TRUE", "true", 1, "1")
d$is_complete <- d$complete_segregation %in% c(TRUE, "True", "TRUE", "true", 1, "1")
d$is_target_like <- d$target_like_Chr23997 %in% c(TRUE, "True", "TRUE", "true", 1, "1")

# One bar per 100-kb window: use the strongest Fisher signal in that window.
split_key <- paste(d$CHROM, d$bin, sep="\t")
idx_by_bin <- split(seq_len(nrow(d)), split_key)
bin_rows <- lapply(idx_by_bin, function(ii) {
  yy <- d$neglog10p[ii]
  jj <- ii[which.max(yy)]
  data.frame(
    CHROM=d$CHROM[jj],
    bin=d$bin[jj],
    bin_start=d$bin_start[jj],
    bin_end=d$bin_end[jj],
    x_mid=d$x_mid[jj],
    y=max(yy, na.rm=TRUE),
    n_pav=length(ii),
    has_group80=any(d$is_group80[ii]),
    has_complete=any(d$is_complete[ii]),
    has_target_like=any(d$is_target_like[ii]),
    top_FinalSV_ID=d$FinalSV_ID[jj],
    top_SVTYPE=d$SVTYPE[jj],
    top_fisher_p=d$fisher_p[jj],
    stringsAsFactors=FALSE
  )
})
b <- do.call(rbind, bin_rows)
b$chr_idx <- match(b$CHROM, chroms)
b$bar_col <- ifelse(b$chr_idx %% 2 == 1, "#79BFEA", "#EE4B8B")
# Keep all window bars in chromosome colors; overlay group-specific and complete windows as black/orange ticks.

write.table(b, file=file.path(out, "results/svgap_supported_PAV_trait_cosegregation.binned100kb.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

manual_x <- (off[manual$chrom[1]] + manual$pos[1]) / 1e6
manual_y <- manual$neglog10p_core[1]
chr_total <- sum(f$len)/1e6
chr_breaks <- cumsum(f$len)/1e6
chr_mids <- f$mid/1e6
bar_width <- (bin_size/1e6) * 0.92

plot_bars <- function(file, type=c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") {
    grDevices::pdf(file, width=8.8, height=3.2, useDingbats=FALSE)
  } else {
    grDevices::png(file, width=8.8, height=3.2, units="in", res=600)
  }
  par(mar=c(4.2,4.5,0.8,1.0), xaxs="i", yaxs="i")
  ymax <- max(4.0, b$y, manual_y, na.rm=TRUE) + 0.55
  plot(NA, NA, xlim=c(0, chr_total), ylim=c(0, ymax), axes=FALSE,
       xlab="Chromosome", ylab=expression(-log[10](italic(P))))
  # Vertical bar-like Manhattan view; far fewer graphical objects than plotting all PAVs.
  segments(b$x_mid, 0, b$x_mid, b$y, col=adjustcolor(b$bar_col, alpha.f=0.78), lwd=1.05, lend="butt")
  # Highlight bins containing 80% group-specific or complete segregation PAVs without drawing all candidates.
  g80 <- b[b$has_group80 & !b$has_complete,]
  if (nrow(g80) > 0) segments(g80$x_mid, 0, g80$x_mid, g80$y, col="#2B2B2B", lwd=1.05, lend="butt")
  comp <- b[b$has_complete,]
  if (nrow(comp) > 0) segments(comp$x_mid, 0, comp$x_mid, comp$y, col="#E69F00", lwd=1.20, lend="butt")
  abline(h=-log10(0.05), lty=2, col="#8C8C8C", lwd=0.8)
  abline(h=-log10(0.01), lty=3, col="#555555", lwd=0.8)
  abline(v=chr_breaks[-length(chr_breaks)], col="#D9D9D9", lwd=0.5)
  points(manual_x, manual_y, pch=23, bg="#D55E00", col="black", cex=1.05, lwd=0.8)
  text(manual_x, min(ymax - 0.12, manual_y + 0.35), labels="Chr23997", font=3, cex=0.78, pos=3)
  arrows(manual_x, min(ymax - 0.25, manual_y + 0.29), manual_x, manual_y + 0.05, length=0.055, lwd=0.8)
  axis(2, las=1, lwd=0.8, lwd.ticks=0.8)
  axis(1, at=chr_mids, labels=gsub("Chr", "", chroms), lwd=0.8, lwd.ticks=0.8)
  mtext("P = 0.05", side=4, at=-log10(0.05), las=1, line=0.10, cex=0.50, col="#666666")
  mtext("P = 0.01", side=4, at=-log10(0.01), las=1, line=0.10, cex=0.50, col="#444444")
  legend("topright", bty="n", cex=0.68,
         legend=c("100-kb windows, max PAV signal", "80% group-specific window", "Complete segregation window", "Chr23997 manual IGV-corrected"),
         col=c("#79BFEA", "#2B2B2B", "#E69F00", "black"),
         pt.bg=c(NA, NA, NA, "#D55E00"),
         pch=c(NA, NA, NA, 23), lty=c(1,1,1,NA), lwd=c(2.2,2.2,2.4,NA), pt.cex=c(NA,NA,NA,1.0))
  box(bty="l", lwd=0.8)
  dev.off()
}
figbase <- file.path(figdir, "svgap_supported_PAV_Chr1_Chr8_trait_cosegregation.neglog10P.binned100kb_bars")
plot_bars(paste0(figbase, ".pdf"), "pdf")
plot_bars(paste0(figbase, ".png"), "png")

cat("output_pdf\t", paste0(figbase, ".pdf"), "\n", sep="")
cat("output_png\t", paste0(figbase, ".png"), "\n", sep="")
cat("source_binned_tsv\t", file.path(out, "results/svgap_supported_PAV_trait_cosegregation.binned100kb.tsv"), "\n", sep="")
cat("input_PAV\t", nrow(d), "\n", sep="")
cat("binned_windows\t", nrow(b), "\n", sep="")
cat("bin_size_bp\t", bin_size, "\n", sep="")
cat("manual_chr23997_neglog10p\t", manual_y, "\n", sep="")
