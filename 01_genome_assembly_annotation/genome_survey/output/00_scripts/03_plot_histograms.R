args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: 03_plot_histograms.R plot_dir histo...")
plot_dir <- args[1]
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
for (h in args[-1]) {
  dat <- tryCatch(read.table(h, header=FALSE), error=function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) next
  colnames(dat) <- c("depth", "count")
  dat <- dat[dat$depth >= 1, ]
  if (nrow(dat) == 0) next
  max_count <- max(dat$count, na.rm=TRUE)
  sub <- dat[dat$depth <= 300 & dat$count > 0, ]
  png(file.path(plot_dir, paste0(basename(h), ".png")), width=1200, height=800)
  plot(sub$depth, sub$count, type="l", lwd=1.2, col="#1f77b4", xlab="k-mer depth", ylab="number of distinct k-mers", main=basename(h))
  grid(col="#dddddd")
  dev.off()
}
