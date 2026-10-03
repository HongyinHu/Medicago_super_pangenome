#!/usr/bin/env Rscript
# Extended Data Fig. 6b (revised): PAV detection frequency by size, DEL/INS only.
# Size bins: 50-499, 500-4,999, >=5,000 bp (left-closed, non-overlapping).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_PAV_detection_frequency_by_size_DEL_INS.R <pav.tsv> <out_table.tsv> <unused_out_prefix> <exclude_sample>")
}
pav_file <- args[1]
out_table <- args[2]
exclude_sample <- args[4]
expected_total <- 609145

meta_cols <- c("SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "read_SV_ID", "svgap_ids", "final_confidence")
dat <- read.delim(pav_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
sample_cols <- setdiff(names(dat), c(meta_cols, exclude_sample))
pav_mat <- as.matrix(dat[, sample_cols, drop = FALSE])
pav_num <- suppressWarnings(matrix(as.numeric(pav_mat), nrow = nrow(pav_mat), dimnames = dimnames(pav_mat)))
pav_num[is.na(pav_num)] <- 0
freq <- rowSums(pav_num > 0)
svlen <- suppressWarnings(abs(as.numeric(dat$SVLEN)))
svtype <- dat$SVTYPE

# PAVs = DEL/INS present in at least one non-reference genome.
keep <- freq > 0 & svtype %in% c("DEL", "INS") & !is.na(svlen) & svlen >= 50
n_del <- sum(keep & svtype == "DEL")
n_ins <- sum(keep & svtype == "INS")
if (sum(keep) != expected_total) {
  stop(sprintf("PAV total %d != expected %d", sum(keep), expected_total))
}

size_levels <- c("50-499", "500-4999", ">=5000")
size_bin <- cut(svlen[keep], breaks = c(50, 500, 5000, Inf), labels = size_levels, right = FALSE)
freq_levels <- seq_along(sample_cols)
count_tab <- as.data.frame(table(
  size_bin = factor(size_bin, levels = size_levels),
  detection_frequency = factor(freq[keep], levels = freq_levels)
), stringsAsFactors = FALSE)
names(count_tab)[3] <- "event_count"
count_tab$total_in_size_bin <- ave(count_tab$event_count, count_tab$size_bin, FUN = sum)
count_tab$percent_of_size_bin <- count_tab$event_count / count_tab$total_in_size_bin
count_tab$percent_label <- sprintf("%.4f", count_tab$percent_of_size_bin)
write.table(count_tab, out_table, sep = "\t", quote = FALSE, row.names = FALSE)

low_freq <- sum(freq[keep] <= 2)
low_freq_pct <- 100 * low_freq / sum(keep)
bin_totals <- tapply(count_tab$event_count, count_tab$size_bin, sum)[size_levels]
low_by_bin <- tapply(count_tab$event_count[count_tab$detection_frequency %in% c("1", "2")],
                     count_tab$size_bin[count_tab$detection_frequency %in% c("1", "2")], sum)[size_levels]
summary_file <- sub("\\.tsv$", ".summary.tsv", out_table)
summary_tab <- data.frame(
  metric = c("input_pav", "excluded_sample", "remaining_species_count", "svtypes_used",
             "PAV_total", "PAV_DEL", "PAV_INS", "excluded_DUP_present", "excluded_INV_present",
             paste0("bin_total_", size_levels),
             "PAV_in_1_or_2_genomes", "PAV_in_1_or_2_genomes_percent",
             paste0("bin_", size_levels, "_in_1_or_2_percent")),
  value = c(pav_file, exclude_sample, length(sample_cols), "DEL,INS",
            sum(keep), n_del, n_ins,
            sum(freq > 0 & svtype == "DUP"), sum(freq > 0 & svtype == "INV"),
            bin_totals, low_freq, sprintf("%.2f", low_freq_pct),
            sprintf("%.2f", 100 * low_by_bin / bin_totals))
)
write.table(summary_tab, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
# Figure is drawn from out_table by plot_PAV_detection_frequency_by_size_DEL_INS_figure.R
# (needs Times New Roman with en dash / >= glyphs, unavailable on this server).
cat("DONE\n")
print(summary_tab, row.names = FALSE)
