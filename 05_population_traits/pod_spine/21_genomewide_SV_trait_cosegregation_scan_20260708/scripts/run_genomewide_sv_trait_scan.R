run_dir <- "path/to/project/N_4.pod_spiny/21_genomewide_SV_trait_cosegregation_scan_20260708"
mat_file <- file.path(run_dir, "source_data", "paper_style_unified_panSV.Msa_ref.PAV.matrix.tsv")
fisher_file <- file.path(run_dir, "source_data", "Chr23997.truth_calibration_summary.tsv")
fig_dir <- file.path(run_dir, "figures")
res_dir <- file.path(run_dir, "results")
sum_dir <- file.path(run_dir, "summary")

d <- read.delim(mat_file, stringsAsFactors=FALSE, check.names=FALSE)
spiny <- c("genome_410","genome_474","genome_Mpo","genome_R108","genome_436","genome_457","genome_454")
spineless <- c("genome_395","genome_461","genome_468","genome_472","genome_M22","genome_Mar","genome_Mru")
# ZM4 is treated as a known exception and is not included in the main scan.
# Msa self, M46 low-quality, and 482 uncertain are excluded.
spiny <- spiny[spiny %in% names(d)]
spineless <- spineless[spineless %in% names(d)]

parse01 <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.integer(x))
  out[!(x %in% c("0","1"))] <- NA_integer_
  out
}
A <- as.data.frame(lapply(d[spiny], parse01))
B <- as.data.frame(lapply(d[spineless], parse01))
a <- rowSums(A == 1, na.rm=TRUE); b <- rowSums(A == 0, na.rm=TRUE)
c <- rowSums(B == 1, na.rm=TRUE); e <- rowSums(B == 0, na.rm=TRUE)
nA <- a + b; nB <- c + e
key <- paste(a,b,c,e,sep="_")
uk <- unique(key)
pmap <- setNames(rep(NA_real_, length(uk)), uk)
for (kk in uk) {
  parts <- as.integer(strsplit(kk, "_", fixed=TRUE)[[1]])
  aa <- parts[1]; bb <- parts[2]; cc <- parts[3]; ee <- parts[4]
  if ((aa+bb) >= 5 && (cc+ee) >= 5) {
    pmap[kk] <- fisher.test(matrix(c(aa,bb,cc,ee), nrow=2, byrow=TRUE), alternative="two.sided")$p.value
  }
}
p <- unname(pmap[key])
neglog <- -log10(p)
neglog[!is.finite(neglog)] <- NA_real_
freqA <- a / nA; freqB <- c / nB
freqDiff <- freqB - freqA

out <- data.frame(
  FinalSV_ID=d$FinalSV_ID, source_methods=d$source_methods, reference=d$reference,
  CHROM=d$CHROM, POS=as.integer(d$POS), END=as.integer(d$END), SVTYPE=d$SVTYPE,
  SVLEN=as.integer(d$SVLEN), read_SV_ID=d$read_SV_ID, final_confidence=d$final_confidence,
  spiny_present=a, spiny_absent=b, spiny_valid=nA, spineless_present=c, spineless_absent=e, spineless_valid=nB,
  spiny_freq=freqA, spineless_freq=freqB, freq_diff_spineless_minus_spiny=freqDiff,
  fisher_p=p, neglog10_p=neglog,
  direction=ifelse(is.na(freqDiff), NA, ifelse(freqDiff>0, "spineless_enriched", ifelse(freqDiff<0, "spiny_enriched", "no_difference"))),
  stringsAsFactors=FALSE
)
out$tested <- !is.na(out$fisher_p)
out$complete_segregation <- out$tested & ((out$spiny_present==0 & out$spineless_absent==0) | (out$spiny_absent==0 & out$spineless_present==0))
out$group_specific_80 <- out$tested & ((out$spiny_freq<=0.2 & out$spineless_freq>=0.8) | (out$spiny_freq>=0.8 & out$spineless_freq<=0.2))
out$target_Chr23997 <- out$FinalSV_ID == "Msa_panSV_000339055" | out$read_SV_ID == "0_0_pbsv.DEL.130205"

fs <- read.delim(fisher_file, stringsAsFactors=FALSE, check.names=FALSE)
manual_p_core <- fs$fisher_p_DEL_presence[fs$model == "core_clear_no_ambiguous_no_A17"]
manual_p_validated <- fs$fisher_p_DEL_presence[fs$model == "validated_spiny_like_vs_core_spineless_with_A17"]
out$manual_IGV_corrected_p <- NA_real_
out$manual_IGV_corrected_neglog10_p <- NA_real_
out$manual_IGV_note <- ""
out$manual_IGV_corrected_p[out$target_Chr23997] <- manual_p_core
out$manual_IGV_corrected_neglog10_p[out$target_Chr23997] <- -log10(manual_p_core)
out$manual_IGV_note[out$target_Chr23997] <- sprintf("Raw matrix contains known false-positive presence in 436/457/474; manual IGV core model P=%.4g; validated+A17 P=%.4g", manual_p_core, manual_p_validated)

write.table(out, file.path(res_dir, "Msa_ref_panSV_trait_cosegregation_fisher.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
ord <- order(out$fisher_p, -abs(out$freq_diff_spineless_minus_spiny), na.last=NA)
write.table(out[ord[seq_len(min(5000,length(ord)))], ], file.path(res_dir, "Msa_ref_panSV_trait_cosegregation_top5000.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(out[out$complete_segregation & out$tested, ][order(out[out$complete_segregation & out$tested, ]$fisher_p), ], file.path(res_dir, "Msa_ref_panSV_complete_cosegregation_sites.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(out[out$group_specific_80 & out$tested, ][order(out[out$group_specific_80 & out$tested, ]$fisher_p), ], file.path(res_dir, "Msa_ref_panSV_group_specific80_sites.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Chromosome order and cumulative coordinates.
chr_num <- suppressWarnings(as.integer(sub("^Chr", "", out$CHROM)))
chr_levels <- unique(out$CHROM[order(ifelse(is.na(chr_num), 999, chr_num), out$CHROM)])
chr_levels <- chr_levels[!is.na(chr_levels) & chr_levels != ""]
chr_max <- tapply(out$END, out$CHROM, max, na.rm=TRUE)
chr_max <- chr_max[chr_levels]
offset <- c(0, cumsum(as.numeric(chr_max))[seq_len(length(chr_max)-1)])
names(offset) <- chr_levels
out$cum_pos <- out$POS + offset[out$CHROM]
chr_mid <- offset + as.numeric(chr_max)/2

# Summary.
summ <- data.frame(
  metric=c("total_SV", "tested_SV", "complete_cosegregation_SV", "group_specific80_SV", "p_lt_0.05", "p_lt_0.01", "Chr23997_raw_p", "Chr23997_manual_core_p"),
  value=c(nrow(out), sum(out$tested), sum(out$complete_segregation, na.rm=TRUE), sum(out$group_specific_80, na.rm=TRUE), sum(out$fisher_p < 0.05, na.rm=TRUE), sum(out$fisher_p < 0.01, na.rm=TRUE), out$fisher_p[out$target_Chr23997][1], manual_p_core)
)
write.table(summ, file.path(sum_dir, "Msa_ref_panSV_trait_cosegregation_scan_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

plot_scan <- function(outfile, device="pdf", ymode="neglog", chr23997_manual=TRUE, local=FALSE) {
  if (device=="pdf") pdf(outfile, width=8.8, height=3.6, useDingbats=FALSE, family="Helvetica")
  if (device=="png") png(outfile, width=8.8, height=3.6, units="in", res=600, type="cairo")
  op <- par(no.readonly=TRUE); on.exit({par(op); dev.off()}, add=TRUE)
  par(family="Helvetica", mar=c(3.2,3.6,0.8,0.8), xaxs="i", yaxs="i")
  dd <- out[out$tested & !is.na(out$cum_pos), ]
  if (local) dd <- dd[dd$CHROM=="Chr4" & dd$POS >= 86000000 & dd$POS <= 92000000, ]
  y <- if (ymode=="p") dd$fisher_p else dd$neglog10_p
  ymax <- if (ymode=="p") 1 else max(y, dd$manual_IGV_corrected_neglog10_p, na.rm=TRUE) * 1.12
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  if (ymode=="p") {
    plot(NA, xlim=range(dd$cum_pos, na.rm=TRUE), ylim=c(1,0), axes=FALSE, xlab="", ylab="Fisher exact P")
  } else {
    plot(NA, xlim=range(dd$cum_pos, na.rm=TRUE), ylim=c(0,ymax), axes=FALSE, xlab="", ylab=expression(-log[10](italic(P))))
  }
  ulevels <- if (local) "Chr4" else chr_levels
  for (i in seq_along(ulevels)) {
    ch <- ulevels[i]
    sub <- dd[dd$CHROM==ch, ]
    col <- ifelse(i %% 2 == 1, "#7ec7ee", "#e84a8a")
    points(sub$cum_pos, if (ymode=="p") sub$fisher_p else sub$neglog10_p, pch=16, cex=0.18, col=adjustcolor(col, 0.75))
  }
  # Highlight group-specific/complete co-segregating sites.
  hit <- dd[dd$group_specific_80, ]
  if (nrow(hit)) points(hit$cum_pos, if (ymode=="p") hit$fisher_p else hit$neglog10_p, pch=16, cex=0.28, col=adjustcolor("#333333", 0.75))
  comp <- dd[dd$complete_segregation, ]
  if (nrow(comp)) points(comp$cum_pos, if (ymode=="p") comp$fisher_p else comp$neglog10_p, pch=21, bg="#f2a900", col="#333333", cex=0.48, lwd=0.35)
  # Chr23997 target raw and manual-corrected marker.
  tar <- out[out$target_Chr23997, ][1, ]
  if (nrow(tar)==1 && (!local || (tar$CHROM=="Chr4" && tar$POS >= 86000000 && tar$POS <= 92000000))) {
    tx <- tar$cum_pos
    ty <- if (ymode=="p") tar$fisher_p else tar$neglog10_p
    points(tx, ty, pch=23, bg="white", col="black", cex=1.0, lwd=0.8)
    if (chr23997_manual && ymode != "p") {
      points(tx, tar$manual_IGV_corrected_neglog10_p, pch=23, bg="#d95f02", col="black", cex=1.15, lwd=0.8)
      arrows(tx, tar$manual_IGV_corrected_neglog10_p + 0.35, tx, tar$manual_IGV_corrected_neglog10_p + 0.05, length=0.08, lwd=0.8)
      text(tx, tar$manual_IGV_corrected_neglog10_p + 0.42, "Chr23997", cex=0.68, font=3, adj=c(0.5,0), xpd=NA)
    } else {
      text(tx, ty + ifelse(ymode=="p", -0.04, 0.25), "Chr23997", cex=0.68, font=3, adj=c(0.5,0), xpd=NA)
    }
  }
  axis(2, las=1, cex.axis=0.65, lwd=0.45, tck=-0.015)
  if (local) {
    at <- pretty(dd$POS, n=6)
    axis(1, at=at + offset["Chr4"], labels=sprintf("%.1f", at/1e6), cex.axis=0.65, lwd=0.45, tck=-0.015)
    mtext("Chr4 position (Mb)", side=1, line=2.0, cex=0.72)
  } else {
    axis(1, at=chr_mid[chr_levels], labels=sub("Chr", "", chr_levels), cex.axis=0.65, lwd=0.45, tck=-0.015)
    mtext("genome_Msa chromosomes", side=1, line=2.0, cex=0.72)
  }
  box(bty="l", lwd=0.6)
  if (ymode != "p") abline(h=-log10(0.05), lty=2, lwd=0.45, col="#777777")
  legend("topright", legend=c("All tested SVs", "80% group-specific", "Complete segregation", "Chr23997 manual IGV-corrected"),
         pch=c(16,16,21,23), pt.bg=c(NA,NA,"#f2a900","#d95f02"), col=c("#7ec7ee","#333333","#333333","black"), bty="n", cex=0.58, pt.cex=c(0.7,0.8,0.9,1.0))
}
plot_scan(file.path(fig_dir, "Msa_ref_SV_trait_cosegregation_scan.neglog10P.pdf"), "pdf", "neglog", TRUE, FALSE)
plot_scan(file.path(fig_dir, "Msa_ref_SV_trait_cosegregation_scan.neglog10P.png"), "png", "neglog", TRUE, FALSE)
plot_scan(file.path(fig_dir, "Msa_ref_SV_trait_cosegregation_scan.raw_P.pdf"), "pdf", "p", FALSE, FALSE)
plot_scan(file.path(fig_dir, "Msa_ref_SV_trait_cosegregation_scan.raw_P.png"), "png", "p", FALSE, FALSE)
plot_scan(file.path(fig_dir, "Chr23997_region_SV_trait_cosegregation_scan.neglog10P.pdf"), "pdf", "neglog", TRUE, TRUE)
plot_scan(file.path(fig_dir, "Chr23997_region_SV_trait_cosegregation_scan.neglog10P.png"), "png", "neglog", TRUE, TRUE)

cap <- c(
  "Genome-wide SV-phenotype co-segregation scan using the genome_Msa single-reference panSV matrix.",
  "Each point is a panSV event tested by two-tailed Fisher exact test comparing 7 spiny/spiny-like accessions (genome_410, genome_474, genome_Mpo, genome_R108, genome_436, genome_457, genome_454) with 7 core spineless accessions (genome_395, genome_461, genome_468, genome_472, genome_M22, genome_Mar, genome_Mru).",
  "Rows with fewer than five non-missing genotypes in either group were not tested. Msa self, genome_M46, genome_482 and ZM4 were excluded from the main scan.",
  "Chr23997 is marked twice in the -log10(P) plot: the white diamond is the raw panSV matrix result; the orange diamond is the manual IGV-corrected target deletion statistic. This distinction is required because the raw matrix contains known false-positive presence calls for genome_436/457/474 at this locus."
)
writeLines(cap, file.path(sum_dir, "figure_caption_and_methods_note.txt"))
