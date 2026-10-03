#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_PAV_detection_frequency_by_size.R <pav.tsv> <out_table.tsv> <out_prefix> <exclude_sample>")
}
pav_file <- args[1]
out_table <- args[2]
out_prefix <- args[3]
exclude_sample <- args[4]

meta_cols <- c("SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "read_SV_ID", "svgap_ids", "final_confidence")
dat <- read.delim(pav_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
sample_cols <- setdiff(names(dat), c(meta_cols, exclude_sample))
# Coerce PAV values robustly. Missing/non-numeric values are treated as absent.
pav_mat <- as.matrix(dat[, sample_cols, drop = FALSE])
pav_num <- suppressWarnings(matrix(as.numeric(pav_mat), nrow = nrow(pav_mat), dimnames = dimnames(pav_mat)))
pav_num[is.na(pav_num)] <- 0
freq <- rowSums(pav_num > 0)
svlen <- suppressWarnings(abs(as.numeric(dat$SVLEN)))
svtype <- dat$SVTYPE

# Keep comparative events present in at least one non-reference genome and having an interpretable size.
# TRA records in this catalogue have SVLEN=0 and are excluded from size-bin plots.
keep <- freq > 0 & !is.na(svlen) & svlen >= 50 & svtype != "TRA"
plot_dat <- data.frame(
  SV_ID = dat$SV_ID[keep],
  SVTYPE = svtype[keep],
  SVLEN_abs = svlen[keep],
  detection_frequency = freq[keep],
  stringsAsFactors = FALSE
)
plot_dat$size_bin <- cut(
  plot_dat$SVLEN_abs,
  breaks = c(49, 500, 5000, Inf),
  labels = c("50-500", "500-5000", ">5000"),
  right = TRUE,
  include.lowest = TRUE
)
size_levels <- c("50-500", "500-5000", ">5000")
freq_levels <- 1:length(sample_cols)
count_tab <- as.data.frame(as.table(table(
  size_bin = factor(plot_dat$size_bin, levels = size_levels),
  detection_frequency = factor(plot_dat$detection_frequency, levels = freq_levels)
)), stringsAsFactors = FALSE)
names(count_tab)[3] <- "event_count"
count_tab$event_count <- as.integer(count_tab$event_count)
count_tab$total_in_size_bin <- ave(count_tab$event_count, count_tab$size_bin, FUN = sum)
count_tab$percent_of_size_bin <- ifelse(count_tab$total_in_size_bin > 0, count_tab$event_count / count_tab$total_in_size_bin, 0)
count_tab$percent_label <- sprintf("%.4f", count_tab$percent_of_size_bin)
write.table(count_tab, out_table, sep = "\t", quote = FALSE, row.names = FALSE)

summary_file <- sub("\\.tsv$", ".summary.tsv", out_table)
summary_tab <- data.frame(
  metric = c("input_pav", "excluded_sample", "remaining_species_count", "events_total_in_pav", "events_present_non_reference", "events_used_size_defined", "events_excluded_no_non_reference_presence", "events_excluded_TRA_or_no_size", "size_bins"),
  value = c(
    pav_file,
    exclude_sample,
    length(sample_cols),
    nrow(dat),
    sum(freq > 0),
    nrow(plot_dat),
    sum(freq == 0),
    sum(freq > 0) - nrow(plot_dat),
    paste(size_levels, collapse = ",")
  )
)
write.table(summary_tab, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

mat <- xtabs(percent_of_size_bin ~ size_bin + detection_frequency, data = count_tab)
mat <- mat[size_levels, as.character(freq_levels), drop = FALSE]
cols <- c("50-500" = "#4F76B8", "500-5000" = "#E7A91A", ">5000" = "#63B34E")

plot_one <- function(device, file, width, height, res = 600) {
  if (device == "pdf") {
    pdf(file, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  } else if (device == "png") {
    png(file, width = width, height = height, units = "in", res = res, type = "cairo")
  } else if (device == "svg") {
    svg(file, width = width, height = height, family = "Helvetica")
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); dev.off() }, add = TRUE)
  par(mar = c(4.2, 4.3, 1.0, 1.0), xaxs = "i", yaxs = "i", cex = 0.9)
  ymax <- max(mat, na.rm = TRUE) * 1.12
  bp <- barplot(
    mat,
    beside = TRUE,
    col = cols[rownames(mat)],
    border = NA,
    ylim = c(0, max(0.8, ymax)),
    axes = FALSE,
    xlab = "Detection frequency",
    ylab = "Percent of variants",
    cex.lab = 1.05,
    cex.names = 0.85,
    names.arg = freq_levels
  )
  axis(2, las = 1, at = seq(0, max(0.8, ymax), by = 0.1), cex.axis = 0.85, lwd = 0.7)
  axis(1, at = colMeans(bp), labels = freq_levels, tick = FALSE, cex.axis = 0.85)
  box(bty = "l", lwd = 0.7)
  legend(
    "topright",
    legend = names(cols),
    title = "SV Size (bp)",
    fill = cols,
    border = NA,
    bty = "n",
    cex = 0.85,
    title.cex = 0.9
  )
}
plot_one("pdf", paste0(out_prefix, ".pdf"), 4.6, 2.8)
plot_one("png", paste0(out_prefix, ".png"), 4.6, 2.8)
plot_one("svg", paste0(out_prefix, ".svg"), 4.6, 2.8)
cat("DONE\n")
cat(out_table, "\n")
cat(summary_file, "\n")
cat(paste0(out_prefix, ".pdf"), "\n")
cat(paste0(out_prefix, ".png"), "\n")
cat(paste0(out_prefix, ".svg"), "\n")
