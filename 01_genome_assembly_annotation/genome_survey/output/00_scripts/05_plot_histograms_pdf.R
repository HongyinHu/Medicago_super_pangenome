args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: 05_plot_histograms_pdf.R plot_dir histo...")
plot_dir <- args[1]
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
for (h in args[-1]) {
  dat <- tryCatch(read.table(h, header=FALSE), error=function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) next
  colnames(dat) <- c("depth", "count")
  dat <- dat[dat$depth >= 1 & dat$count > 0, ]
  if (nrow(dat) == 0) next
  sub <- dat[dat$depth <= 300, ]
  out <- file.path(plot_dir, paste0(basename(h), ".pdf"))
  pdf(out, width=9, height=6, useDingbats=FALSE)
  plot(sub$depth, sub$count, type="l", lwd=1.1, col="#1f77b4",
       xlab="k-mer depth", ylab="number of distinct k-mers",
       main=basename(h))
  grid(col="#dddddd")
  dev.off()
}
