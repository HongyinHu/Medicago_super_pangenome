args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: plot_local_gwas.R <local_association.tsv> <out.pdf>")

assoc <- read.table(args[[1]], header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "")
if (!nrow(assoc)) stop("local association table is empty")
assoc$pos_mb <- assoc$pos / 1e6
palette <- c(SNP = "#7A7A7A", INDEL = "#E07A3F", SV = "#6C5AAE", Candidate_SV = "#C23B22")
point_col <- palette[assoc$class]
point_col[is.na(point_col)] <- "#7A7A7A"

pdf(args[[2]], width = 8.2, height = 5.2, useDingbats = FALSE)
par(mar = c(4.5, 5.0, 1.2, 1.0), las = 1)
plot(assoc$pos_mb, assoc$minus_log10_p, type = "n", xlab = "Chr4 position (Mb)",
     ylab = expression(-log[10](italic(P))), bty = "l")
normal <- assoc$class != "Candidate_SV"
points(assoc$pos_mb[normal], assoc$minus_log10_p[normal], pch = 16, cex = 0.55,
       col = point_col[normal])
candidate <- assoc$class == "Candidate_SV"
points(assoc$pos_mb[candidate], assoc$minus_log10_p[candidate], pch = 23, cex = 1.6,
       bg = palette[["Candidate_SV"]], col = "black")
abline(h = -log10(0.05), lty = 2, col = "#8C8C8C")
abline(h = -log10(0.01), lty = 3, col = "#8C8C8C")
legend("topright", legend = c("SNP (global GEMMA)", "INDEL (global GEMMA)", "SV (global GEMMA)",
                              "Chr23997 221-bp PAV (targeted Fisher)"),
       pch = c(16, 16, 16, 23), pt.bg = c(NA, NA, NA, palette[["Candidate_SV"]]),
       col = c(palette[["SNP"]], palette[["INDEL"]], palette[["SV"]], "black"),
       bty = "n", cex = 0.8)
if (any(candidate)) {
  text(assoc$pos_mb[candidate], assoc$minus_log10_p[candidate], labels = "Chr23997", pos = 3, cex = 0.8)
}
mtext("Local association within +/-100 kb of Chr23997", side = 3, line = 0.2, cex = 0.9)
dev.off()
