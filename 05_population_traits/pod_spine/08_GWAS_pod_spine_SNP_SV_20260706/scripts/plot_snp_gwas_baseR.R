args <- commandArgs(trailingOnly=TRUE)
run <- ifelse(length(args) >= 1, args[1], ".")
figdir <- file.path(run, "results", "figures")
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)
chr23997_start <- 88979818
chr23997_end <- 88982331
read_gwas <- function(path) {
  d <- read.table(path, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE, stringsAsFactors=FALSE)
  d <- d[d$TEST == "ADD" & !is.na(d$P) & d$P > 0, , drop=FALSE]
  names(d)[names(d)=="#CHROM"] <- "CHROM"
  d$CHROM <- as.character(d$CHROM)
  d$chr_num <- suppressWarnings(as.integer(gsub("^Chr", "", d$CHROM)))
  d <- d[!is.na(d$chr_num) & d$chr_num >= 1 & d$chr_num <= 8, , drop=FALSE]
  d$POS <- as.numeric(d$POS)
  d$mlog10p <- -log10(d$P)
  d
}
plot_manhattan <- function(d, prefix, title) {
  chrmax <- tapply(d$POS, d$chr_num, max, na.rm=TRUE)
  chrmax <- chrmax[as.character(1:8)]
  chrmax[is.na(chrmax)] <- 0
  offsets <- c(0, cumsum(chrmax))[1:8]
  names(offsets) <- as.character(1:8)
  d$x <- d$POS + offsets[as.character(d$chr_num)]
  centers <- offsets + chrmax/2
  gene_mid_x <- offsets["4"] + mean(c(chr23997_start, chr23997_end))
  ymax <- max(ceiling(max(d$mlog10p, na.rm=TRUE)), 10)
  for(ext in c("png", "pdf")) {
    out <- file.path(figdir, paste0(prefix, ".manhattan.", ext))
    if(ext=="png") png(out, width=1800, height=850, res=180) else pdf(out, width=10, height=4.8)
    par(mar=c(4.2,5,3,1), family="serif")
    plot(NA, xlim=range(d$x), ylim=c(0, ymax), xaxt="n", xlab="Chromosome", ylab=expression(-log[10](italic(P))), main=title, bty="l")
    cols <- c("#9da3a6", "#5d6670")
    for(chr in 1:8) {
      dd <- d[d$chr_num==chr,]
      points(dd$x, dd$mlog10p, pch=20, cex=0.35, col=adjustcolor(cols[(chr %% 2)+1], 0.65))
    }
    abline(h=-log10(5e-8), col="#d95f02", lty=2, lwd=1.2)
    abline(v=gene_mid_x, col="#1b9e77", lwd=1.3)
    text(gene_mid_x, ymax*0.95, "Chr23997", col="#1b9e77", srt=90, adj=c(1,0.5), cex=0.8)
    axis(1, at=centers, labels=1:8, tick=FALSE)
    dev.off()
  }
}
plot_zoom <- function(d, prefix, title) {
  z <- d[d$chr_num==4 & d$POS >= 88000000 & d$POS <= 90000000, , drop=FALSE]
  ymax <- max(ceiling(max(z$mlog10p, na.rm=TRUE)), 5)
  top <- z[order(z$P), ][1,]
  for(ext in c("png", "pdf")) {
    out <- file.path(figdir, paste0(prefix, ".Chr4_88_90Mb_zoom.", ext))
    if(ext=="png") png(out, width=1500, height=850, res=180) else pdf(out, width=8.5, height=4.8)
    par(mar=c(4.5,5,3,1), family="serif")
    plot(z$POS/1e6, z$mlog10p, pch=20, cex=0.55, col=adjustcolor("#5d6670",0.65), xlab="Chr4 position (Mb)", ylab=expression(-log[10](italic(P))), main=title, bty="l", ylim=c(0,ymax))
    rect(chr23997_start/1e6, 0, chr23997_end/1e6, ymax, col=adjustcolor("#1b9e77",0.15), border=NA)
    points(z$POS/1e6, z$mlog10p, pch=20, cex=0.55, col=adjustcolor("#5d6670",0.65))
    abline(h=-log10(5e-8), col="#d95f02", lty=2, lwd=1.2)
    points(top$POS/1e6, top$mlog10p, pch=21, bg="#e7298a", col="black", cex=1.2)
    text(mean(c(chr23997_start, chr23997_end))/1e6, ymax*0.92, "Chr23997", col="#1b9e77", cex=0.9)
    text(top$POS/1e6, top$mlog10p, labels=paste0("top ", top$POS, "\nP=", format(top$P, scientific=TRUE, digits=3)), pos=3, cex=0.7)
    dev.off()
  }
}
write_top <- function(d, prefix) {
  top <- d[order(d$P), ][1:50, c("CHROM","POS","ID","REF","ALT","A1","OBS_CT","OR","Z_STAT","P","mlog10p")]
  write.table(top, file=file.path(run,"results",paste0(prefix,".top50.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
}
no_file <- file.path(run,"results","sativa182_pod_spine_snp_gwas.corrected.no_covar.pod_spine.glm.logistic")
pc_file <- file.path(run,"results","sativa182_pod_spine_snp_gwas.corrected.PC5.pod_spine.glm.logistic")
no <- read_gwas(no_file); pc <- read_gwas(pc_file)
plot_manhattan(no, "sativa182_pod_spine_snp_gwas.corrected.no_covar", "Pod spine SNP GWAS, no covariate")
plot_zoom(no, "sativa182_pod_spine_snp_gwas.corrected.no_covar", "Chr4 88-90 Mb, no covariate")
write_top(no, "sativa182_pod_spine_snp_gwas.corrected.no_covar")
plot_manhattan(pc, "sativa182_pod_spine_snp_gwas.corrected.PC5", "Pod spine SNP GWAS, PC1-PC5")
plot_zoom(pc, "sativa182_pod_spine_snp_gwas.corrected.PC5", "Chr4 88-90 Mb, PC1-PC5")
write_top(pc, "sativa182_pod_spine_snp_gwas.corrected.PC5")
