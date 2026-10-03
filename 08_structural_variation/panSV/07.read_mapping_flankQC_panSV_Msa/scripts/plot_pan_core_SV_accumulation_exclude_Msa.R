#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) != 1 || hit == length(args)) stop(paste('Missing', flag))
  args[[hit + 1]]
}
input <- get_arg('--input')
out_prefix <- get_arg('--out-prefix')

df <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
fmt_y <- function(x) ifelse(x >= 1000, paste0(round(x/1000), 'k'), as.character(x))
plot_one <- function() {
  old <- par(no.readonly=TRUE); on.exit(par(old), add=TRUE)
  par(family='serif', mar=c(4.2,4.6,1.0,1.1), xaxs='i', yaxs='i')
  x <- df$sample_number
  pan <- df$pan_SVs_mean
  core <- df$core_SVs_mean
  ymax <- ceiling(max(pan) * 1.08 / 100000) * 100000
  plot(NA, xlim=c(0, max(x)+1), ylim=c(0,ymax), axes=FALSE, xlab='', ylab='SV number')
  axis(1, at=seq(0, max(x), by=5), tck=-0.015)
  yt <- seq(0, ymax, by=100000)
  axis(2, at=yt, labels=fmt_y(yt), las=1, tck=-0.015)
  abline(h=yt, col=adjustcolor('grey70',0.18), lwd=0.6)
  barw <- 0.78
  rect(x - barw/2, core, x + barw/2, pan, col='#5B8DB8', border=NA)
  rect(x - barw/2, 0, x + barw/2, core, col='#D94245', border=NA)
  box(bty='l', lwd=1.1)
  mtext('Sample numbers', side=1, line=2.7, cex=1.0)
  legend('topleft', legend=c('pan-SVs','core-SVs'), fill=c('#5B8DB8','#D94245'), border=NA, bty='n', cex=0.9)
}

dir.create(dirname(out_prefix), recursive=TRUE, showWarnings=FALSE)
for (ext in c('pdf','svg','png')) {
  out <- paste0(out_prefix, '.', ext)
  if (ext == 'pdf') pdf(out, width=4.1, height=4.3, useDingbats=FALSE)
  if (ext == 'svg') svg(out, width=4.1, height=4.3, pointsize=10)
  if (ext == 'png') png(out, width=4.1, height=4.3, units='in', res=500, type='cairo-png')
  plot_one()
  dev.off()
  cat('WROTE\t', out, '\n', sep='')
}