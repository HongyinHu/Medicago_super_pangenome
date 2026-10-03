# Base-R publication-style figure for the Chr23997 target intron deletion.
# Inputs are validation-calibrated tables produced under N_4.pod_spiny/19_...

outdir <- "path/to/project/N_4.pod_spiny/20_Chr23997_candidate_figure_20260708"
truth_file <- file.path(outdir, "source_data", "Chr23997.target_large_intron_DEL.truth_calibrated_per_species.tsv")
fisher_file <- file.path(outdir, "source_data", "Chr23997.truth_calibration_summary.tsv")
exact_file <- file.path(outdir, "source_data", "Chr23997.target_DEL_exactness_by_sample_caller.tsv")
figdir <- file.path(outdir, "figures")
sumdir <- file.path(outdir, "summary")

d <- read.delim(truth_file, stringsAsFactors=FALSE, check.names=FALSE)
fs <- read.delim(fisher_file, stringsAsFactors=FALSE, check.names=FALSE)
ex <- read.delim(exact_file, stringsAsFactors=FALSE, check.names=FALSE)

species_name <- c(
  genome_Msa="M. sativa T2T", genome_ZM4="M. sativa ZM4", genome_468="M. cretacea",
  genome_457="M. marina", genome_410="M. praecox", genome_R108="M. truncatula R108",
  genome_A17="M. truncatula A17", genome_Mpo="M. polymorpha", genome_M22="M. orbicularis",
  genome_461="M. secundiflora", genome_474="M. carstiensis", genome_472="M. suffruticosa",
  genome_395="M. lupulina", genome_482="M. edgeworthii", genome_M46="M. fischeriana",
  genome_436="M. radiata", genome_Mar="M. archiducis-nicolai", genome_Mru="M. ruthenica",
  genome_454="M. lanigera"
)
# Order follows the provided phylogeny screenshot, with A17 placed next to R108 as an extra accession.
phy_order <- c("genome_Msa","genome_ZM4","genome_468","genome_457","genome_410","genome_R108","genome_A17",
               "genome_Mpo","genome_M22","genome_461","genome_474","genome_472","genome_395",
               "genome_482","genome_M46","genome_436","genome_Mar","genome_Mru","genome_454")
d <- d[match(phy_order, d$sample), ]
d$latin <- unname(species_name[d$sample])
d$plot_label <- paste0(d$sample, "  ", d$latin)
d$trait_simple <- ifelse(grepl("spiny", d$phenotype_group_curated) & !grepl("spineless", d$phenotype_group_curated), "Spiny-like",
                         ifelse(grepl("spineless", d$phenotype_group_curated), "Spineless", "Excluded/uncertain"))
d$trait_simple[d$sample %in% c("genome_M46","genome_482","genome_Msa")] <- "Excluded/uncertain"
d$truth_simple <- ifelse(d$validated_truth_class == "DEL_present", "Target DEL",
                         ifelse(d$validated_truth_class == "INTACT_absent", "Intact",
                                ifelse(d$validated_truth_class == "UNCERTAIN", "Uncertain", "Excluded")))
d$truth_simple[d$sample == "genome_ZM4"] <- "Intact exception"
d$truth_simple[d$sample %in% c("genome_M46","genome_482","genome_Msa")] <- "Excluded/uncertain"

spiny_set <- c("genome_410","genome_474","genome_Mpo","genome_R108","genome_436","genome_457","genome_454","genome_A17")
spineless_core <- c("genome_395","genome_461","genome_468","genome_472","genome_M22","genome_Mar","genome_Mru")
spineless_plus_zm4 <- c(spineless_core, "genome_ZM4")

count_group <- function(samples) {
  sub <- d[d$sample %in% samples, ]
  del <- sum(sub$truth_DEL_binary == 1, na.rm=TRUE)
  intact <- sum(sub$truth_DEL_binary == 0, na.rm=TRUE)
  c("Target DEL"=del, "Intact"=intact)
}
bmat <- rbind(
  "Spiny-like\n(n=8)" = count_group(spiny_set),
  "Spineless core\n(n=7)" = count_group(spineless_core),
  "Spineless + ZM4\n(n=8)" = count_group(spineless_plus_zm4)
)
bprop <- bmat / rowSums(bmat) * 100

# Exact caller support matrix.
callers <- c("pbsv", "cutesv", "sniffles2")
val_samples <- c(spiny_set, spineless_core, "genome_ZM4")
val_samples <- val_samples[val_samples %in% d$sample]
heat <- matrix("No exact target record", nrow=length(val_samples), ncol=5,
               dimnames=list(val_samples, c("IGV truth", "Genotyped VCF", "pbsv exact", "cuteSV exact", "Sniffles2 exact")))
for (s in val_samples) {
  rr <- d[d$sample == s, ]
  heat[s, "IGV truth"] <- ifelse(rr$truth_DEL_binary == 1, "Target DEL", "Intact")
  gt <- rr$target_genotyped_GT
  heat[s, "Genotyped VCF"] <- ifelse(gt %in% c("0/1", "1/1", "1|1", "0|1"), "Target DEL", ifelse(gt %in% c("0/0", "0|0"), "Intact", "Missing"))
  for (cl in callers) {
    ee <- ex[ex$sample == s & ex$caller == cl, ]
    has <- if (nrow(ee)) any(ee$exact_target_like_DEL_records > 0, na.rm=TRUE) else FALSE
    colname <- if (cl == "cutesv") "cuteSV exact" else if (cl == "sniffles2") "Sniffles2 exact" else "pbsv exact"
    heat[s, colname] <- ifelse(has, "Target DEL", "No exact target record")
  }
}
write.table(d, file.path(sumdir, "Chr23997_figure_plotted_species_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(bmat, file.path(sumdir, "Chr23997_figure_group_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(heat, file.path(sumdir, "Chr23997_figure_validation_heatmap_matrix.tsv"), sep="\t", quote=FALSE, col.names=NA)

cols_trait <- c("Spiny-like"="#d95f5f", "Spineless"="#67aeb3", "Excluded/uncertain"="#bdbdbd")
cols_state <- c("Target DEL"="#d95f02", "Intact"="#2b6aa8", "Intact exception"="#7b7b7b", "Uncertain"="#bdbdbd", "Excluded/uncertain"="#e0e0e0", "Excluded"="#e0e0e0")
cols_heat <- c("Target DEL"="#d95f02", "Intact"="#2b6aa8", "No exact target record"="#eeeeee", "Missing"="#bdbdbd")

panel_label <- function(x, y, lab) text(x, y, lab, font=2, cex=1.25, adj=c(0,1), xpd=NA)

plot_panel_a <- function() {
  n <- nrow(d); yy <- rev(seq_len(n))
  plot(NA, xlim=c(0, 10), ylim=c(0.2, n+1.2), axes=FALSE, xlab="", ylab="")
  panel_label(0, n+1.1, "a")
  text(0.55, n+0.9, "Phylogeny-ordered validation status", adj=0, font=2, cex=0.86)
  text(5.8, n+0.2, "Pod trait", cex=0.62, font=2)
  text(7.5, n+0.2, "Chr23997\nintron-2 SV", cex=0.62, font=2)
  text(9.0, n+0.2, "Use", cex=0.62, font=2)
  for (i in seq_len(n)) {
    y <- yy[i]
    colline <- ifelse(i %% 2 == 0, "#f7f7f7", "white")
    rect(0.35, y-0.42, 9.75, y+0.42, border=NA, col=colline)
    text(0.55, y, d$plot_label[i], adj=0, cex=0.54)
    points(5.8, y, pch=21, bg=cols_trait[d$trait_simple[i]], col="white", cex=1.1)
    points(7.55, y, pch=22, bg=cols_state[d$truth_simple[i]], col="white", cex=1.15)
    use_lab <- ifelse(d$sample[i] %in% c("genome_M46","genome_482","genome_Msa"), "excluded", ifelse(d$sample[i]=="genome_ZM4", "exception", "tested"))
    text(9.0, y, use_lab, cex=0.48, col=ifelse(use_lab=="tested", "#333333", "#777777"))
  }
  legend(0.55, 1.2, legend=c("spiny-like", "spineless", "excluded/uncertain"), pt.bg=cols_trait, pch=21, pt.cex=1.1, bty="n", cex=0.50, x.intersp=0.75, y.intersp=0.8)
  legend(4.55, 1.2, legend=c("target DEL", "intact", "intact exception"), pt.bg=cols_state[c("Target DEL","Intact","Intact exception")], pch=22, pt.cex=1.1, bty="n", cex=0.50, x.intersp=0.75, y.intersp=0.8)
}

plot_panel_b <- function() {
  plot(NA, xlim=c(0, 100), ylim=c(0.4, 3.9), axes=FALSE, xlab="", ylab="")
  panel_label(0, 3.9, "b")
  text(3, 3.73, "Target DEL frequency by phenotype group", adj=0, font=2, cex=0.86)
  yv <- c(3,2,1)
  for (i in 1:3) {
    rect(0, yv[i]-0.28, bprop[i,"Target DEL"], yv[i]+0.28, col=cols_state["Target DEL"], border=NA)
    rect(bprop[i,"Target DEL"], yv[i]-0.28, 100, yv[i]+0.28, col=cols_state["Intact"], border=NA)
    text(-3, yv[i], rownames(bprop)[i], adj=1, cex=0.62, xpd=NA)
    text(102, yv[i], sprintf("%s/%s DEL", bmat[i,"Target DEL"], sum(bmat[i,])), adj=0, cex=0.58, xpd=NA)
  }
  axis(1, at=seq(0,100,25), labels=paste0(seq(0,100,25), "%"), cex.axis=0.58, lwd=0.5, tck=-0.025)
  mtext("Validated target-DEL frequency", side=1, line=1.7, cex=0.62)
  p1 <- fs$fisher_p_DEL_presence[fs$model == "validated_spiny_like_vs_core_spineless_with_A17"]
  p2 <- fs$fisher_p_DEL_presence[fs$model == "validated_spiny_like_vs_spineless_plus_ZM4_exception"]
  text(48, 0.42, sprintf("Fisher exact: spiny-like vs spineless core, P = %.2g\nIncluding ZM4 exception, P = %.2g", p1, p2), cex=0.56)
  legend(58, 3.1, legend=c("Target DEL", "Intact"), fill=cols_state[c("Target DEL","Intact")], bty="n", cex=0.58)
}

plot_panel_c <- function() {
  x0 <- 88979818; x1 <- 88982331
  del0 <- 88980831; del1 <- 88981052
  plot(NA, xlim=c(x0-80, x1+80), ylim=c(0, 1), axes=FALSE, xlab="", ylab="")
  panel_label(x0-80, 1, "c")
  text(x0+60, 0.92, "Chr23997 gene model and candidate intron deletion", adj=0, font=2, cex=0.86)
  segments(x0, 0.58, x1, 0.58, lwd=1.1, col="#333333")
  arrows(x1-130, 0.58, x1, 0.58, length=0.06, lwd=1.1, col="#333333")
  exons <- data.frame(start=c(88979818, 88980255, 88981335, 88981810), end=c(88980110, 88980375, 88981520, 88982331))
  for (i in seq_len(nrow(exons))) rect(exons$start[i], 0.49, exons$end[i], 0.67, col="#333333", border=NA)
  rect(del0, 0.20, del1, 0.82, col=adjustcolor("#d95f02", 0.20), border=NA)
  segments(c(del0,del1), 0.22, c(del0,del1), 0.84, col="#d95f02", lwd=1.1)
  text(mean(c(del0,del1)), 0.15, "221 bp DEL\n(second intron)", cex=0.60, col="#8a3c00")
  text(x0, 0.35, "Chr4:88,979,818", adj=0, cex=0.54)
  text(x1, 0.35, "88,982,331", adj=1, cex=0.54)
  text(mean(c(x0,x1)), 0.74, "Chr23997", cex=0.66, font=3)
}

plot_panel_d <- function() {
  h <- heat[rev(rownames(heat)), , drop=FALSE]
  nr <- nrow(h); nc <- ncol(h)
  plot(NA, xlim=c(0.5,nc+0.5), ylim=c(0.5,nr+1.5), axes=FALSE, xlab="", ylab="")
  panel_label(0.5, nr+1.5, "d")
  text(0.9, nr+1.32, "Manual truth versus automated target support", adj=0, font=2, cex=0.86)
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    rect(j-0.45, i-0.38, j+0.45, i+0.38, col=cols_heat[h[i,j]], border="white", lwd=0.6)
  }
  text(seq_len(nc), nr+0.55, colnames(h), srt=35, adj=0, cex=0.55, xpd=NA)
  text(0.35, seq_len(nr), rownames(h), adj=1, cex=0.50, xpd=NA)
  legend(3.1, 1.0, legend=names(cols_heat), fill=cols_heat, bty="n", cex=0.52, ncol=2)
  box(col="#dddddd")
}

make_fig <- function(file, device=c("pdf","png","svg")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width=7.2, height=6.2, useDingbats=FALSE, family="Helvetica")
  if (device == "png") png(file, width=7.2, height=6.2, units="in", res=600, type="cairo")
  if (device == "svg") svg(file, width=7.2, height=6.2, family="Helvetica")
  op <- par(no.readonly=TRUE)
  par(family="Helvetica", xaxs="i", yaxs="i", oma=c(0.1,0.1,0.1,0.1))
  layout(matrix(c(1,2,1,3,4,4), nrow=2, byrow=TRUE), widths=c(1.55,1), heights=c(1.25,1))
  par(mar=c(0.6,0.4,1.0,0.6)); plot_panel_a()
  par(mar=c(2.4,3.5,1.0,2.0)); plot_panel_b()
  par(mar=c(1.4,0.6,1.0,0.6)); plot_panel_c()
  par(mar=c(0.8,4.2,1.0,0.8)); plot_panel_d()
  par(op)
  dev.off()
}
make_fig(file.path(figdir, "Chr23997_candidate_SV_evidence_figure.pdf"), "pdf")
make_fig(file.path(figdir, "Chr23997_candidate_SV_evidence_figure.png"), "png")
make_fig(file.path(figdir, "Chr23997_candidate_SV_evidence_figure.svg"), "svg")

caption <- c(
"Figure contract: Chr23997 second-intron deletion is a validated candidate SV associated with pod spiny/spineless status.",
"Panel a: species/status map ordered according to the provided Medicago phylogeny screenshot; this panel is a status map, not a newly inferred tree.",
"Panel b: validated target DEL frequency in curated phenotype groups. Spiny-like accessions: 0/8 DEL. Core spineless accessions: 7/7 DEL. Fisher exact P=1.55e-4. Including the intact spineless ZM4 exception: 7/8 DEL, P=0.0014.",
"Panel c: schematic Chr23997 gene model and the candidate 221 bp DEL in the second intron on the genome_Msa reference coordinate system.",
"Panel d: comparison of manual IGV truth, genotyped VCF, and exact target-like DEL support from pbsv/cuteSV/Sniffles2. Broad nearby caller records are not counted as exact target support."
)
writeLines(caption, file.path(sumdir, "Chr23997_candidate_SV_evidence_figure_caption.txt"))
