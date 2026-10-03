#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

read_plink_unadjusted <- function(path) {
  x <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if ("TEST" %in% names(x)) x <- x[x$TEST == "ADD", , drop = FALSE]
  out <- data.frame(
    CHR = suppressWarnings(as.integer(x$CHR)),
    POS = suppressWarnings(as.numeric(x$BP)),
    P = suppressWarnings(as.numeric(x$P)),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$CHR) & is.finite(out$POS) & is.finite(out$P) & out$P > 0 & out$P <= 1, , drop = FALSE]
  out[order(out$CHR, out$POS), , drop = FALSE]
}

read_gmmat_adjusted <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  out <- data.frame(
    CHR = suppressWarnings(as.integer(x$CHR)),
    POS = suppressWarnings(as.numeric(x$POS)),
    P = suppressWarnings(as.numeric(x$PVAL)),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$CHR) & is.finite(out$POS) & is.finite(out$P) & out$P > 0 & out$P <= 1, , drop = FALSE]
  out[order(out$CHR, out$POS), , drop = FALSE]
}

prepare_manhattan_data <- function(data) {
  data <- data[order(data$CHR, data$POS), , drop = FALSE]
  chr_levels <- sort(unique(data$CHR))
  chr_max <- vapply(chr_levels, function(chr) max(data$POS[data$CHR == chr]), numeric(1))
  offsets <- c(0, cumsum(chr_max[-length(chr_max)]))
  names(offsets) <- as.character(chr_levels)
  data$BPcum <- data$POS + offsets[as.character(data$CHR)]
  data$chr_label <- factor(as.character(data$CHR), levels = as.character(chr_levels))
  axis <- data.frame(
    CHR = chr_levels,
    center = vapply(chr_levels, function(chr) median(data$BPcum[data$CHR == chr]), numeric(1)),
    stringsAsFactors = FALSE
  )
  list(data = data, axis = axis)
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

theme_gwas <- function() {
  theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7.5),
      plot.caption = element_text(size = 6.5, hjust = 0),
      panel.grid = element_blank()
    )
}

manhattan_plot <- function(data, title, subtitle) {
  prepared <- prepare_manhattan_data(data)
  d <- prepared$data
  d$neglog10p <- -log10(d$P)
  bonferroni <- -log10(0.05 / nrow(d))
  ggplot(d, aes(x = BPcum, y = neglog10p, colour = chr_label)) +
    geom_hline(yintercept = bonferroni, linewidth = 0.35, linetype = "dashed", colour = "#C23B22") +
    geom_point(size = 0.16, alpha = 0.70, stroke = 0) +
    scale_colour_manual(values = rep(c("#3B6EA5", "#9BAEC8"), length.out = nlevels(d$chr_label)), guide = "none") +
    scale_x_continuous(breaks = prepared$axis$center, labels = prepared$axis$CHR, expand = expansion(mult = c(0.01, 0.01))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Chromosome",
      y = expression(-log[10](italic(P))),
      caption = paste0("Dashed line: Bonferroni threshold (0.05/", format(nrow(d), big.mark = ","), ")")
    ) +
    theme_gwas()
}

qq_plot <- function(data, title, subtitle) {
  p <- sort(data$P)
  n <- length(p)
  qq <- data.frame(
    expected = -log10(ppoints(n)),
    observed = -log10(p)
  )
  lim <- max(c(qq$expected, qq$observed)) * 1.03
  ggplot(qq, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed", colour = "#767676") +
    geom_point(size = 0.28, alpha = 0.58, colour = "#3B6EA5") +
    coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(Expected~~-log[10](italic(P))),
      y = expression(Observed~~-log[10](italic(P))),
      caption = paste0("n = ", format(n, big.mark = ","), "; lambdaGC = ", sprintf("%.3f", lambda_gc(p)))
    ) +
    theme_gwas()
}

save_plot_set <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "sans")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 180)
  print(plot)
  grDevices::dev.off()
}

write_calibration <- function(path, raw, adjusted) {
  x <- data.frame(
    model = c("Unadjusted_PLINK_allelic_no_PC_no_GRM", "Adjusted_GMMAT_PC1_5_plus_GRM"),
    n_tested = c(nrow(raw), nrow(adjusted)),
    lambda_gc = c(lambda_gc(raw$P), lambda_gc(adjusted$P)),
    min_p = c(min(raw$P), min(adjusted$P)),
    stringsAsFactors = FALSE
  )
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3) {
    stop("Usage: plot_snp_unadjusted_vs_adjusted.R unadjusted.assoc adjusted.gmmat.tsv output_dir")
  }
  if (!requireNamespace("ragg", quietly = TRUE)) stop("R package ragg is required")
  unadjusted_path <- args[[1]]
  adjusted_path <- args[[2]]
  output_dir <- args[[3]]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  raw <- read_plink_unadjusted(unadjusted_path)
  adjusted <- read_gmmat_adjusted(adjusted_path)
  if (nrow(raw) == 0 || nrow(adjusted) == 0) stop("No valid association P values available for plotting")

  save_plot_set(
    manhattan_plot(raw, "Unadjusted SNP association Manhattan plot", "PLINK allelic test; no PC covariates and no GRM"),
    file.path(output_dir, "SNP_unadjusted_manhattan")
  )
  save_plot_set(
    qq_plot(raw, "Unadjusted SNP association QQ plot", "PLINK allelic test; no PC covariates and no GRM"),
    file.path(output_dir, "SNP_unadjusted_qq")
  )
  save_plot_set(
    manhattan_plot(adjusted, "Structure-adjusted SNP association Manhattan plot", "GMMAT score test; PC1-PC5 plus SNP GRM"),
    file.path(output_dir, "SNP_adjusted_manhattan")
  )
  save_plot_set(
    qq_plot(adjusted, "Structure-adjusted SNP association QQ plot", "GMMAT score test; PC1-PC5 plus SNP GRM"),
    file.path(output_dir, "SNP_adjusted_qq")
  )
  write_calibration(file.path(output_dir, "SNP_unadjusted_vs_adjusted_calibration.tsv"), raw, adjusted)
  cat("Wrote four SNP Manhattan/QQ figures to", output_dir, "\n")
}

if (sys.nframe() == 0) main()
