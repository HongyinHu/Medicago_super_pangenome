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
  x <- df$sample_number
  pan <- df$pan_SVs_expected
  core <- df$core_SVs_expected
  n <- max(x)
  xlim <- c(0, n + 1)
  barw <- 0.78
  col_pan <- '#5B8DB8'
  col_core <- '#D94245'

  lower_step <- 5000
  lower_max <- 20000
  lower_ticks <- seq(0, lower_max, by=lower_step)

  upper_step <- 100000
  upper_min <- 100000
  upper_max <- ceiling(max(pan) * 1.08 / upper_step) * upper_step
  upper_ticks <- seq(upper_min, upper_max, by=upper_step)

  layout(matrix(c(1,2), nrow=2), heights=c(2.55,1.35))

  par(family='serif', mar=c(0.15,4.8,0.9,1.1), xaxs='i', yaxs='i')
  plot(NA, xlim=xlim, ylim=c(upper_min, upper_max), axes=FALSE, xlab='', ylab='')
  rect(x - barw/2, core, x + barw/2, pan, col=col_pan, border=NA)
  rect(x - barw/2, 0, x + barw/2, core, col=col_core, border=NA)
  axis(2, at=upper_ticks, labels=fmt_y(upper_ticks), las=1, tck=-0.015)
  abline(h=upper_ticks, col=adjustcolor('grey70',0.18), lwd=0.6)
  box(bty='l', lwd=1.1)
  mtext('SV number', side=2, line=3.35, cex=1.0)
  legend('topleft', legend=c('Pan', 'Core'), fill=c(col_pan, col_core), border=NA, bty='n', cex=0.9)

  par(family='serif', mar=c(4.2,4.8,0.05,1.1), xaxs='i', yaxs='i')
  plot(NA, xlim=xlim, ylim=c(0, lower_max), axes=FALSE, xlab='', ylab='')
  rect(x - barw/2, core, x + barw/2, pan, col=col_pan, border=NA)
  rect(x - barw/2, 0, x + barw/2, core, col=col_core, border=NA)
  axis(2, at=lower_ticks, labels=fmt_y(lower_ticks), las=1, tck=-0.015, cex.axis=0.78)
  axis(1, at=seq(0, n, by=5), tck=-0.015)
  abline(h=lower_ticks, col=adjustcolor('grey70',0.16), lwd=0.5)
  box(bty='l', lwd=1.1)
  mtext('Sample numbers', side=1, line=2.7, cex=1.0)

  par(fig=c(0,1,0,1), new=TRUE, mar=c(0,0,0,0))
  plot.new()
  segments(0.085, 0.350, 0.103, 0.373, xpd=NA, lwd=1.0)
  segments(0.085, 0.368, 0.103, 0.391, xpd=NA, lwd=1.0)
}

dir.create(dirname(out_prefix), recursive=TRUE, showWarnings=FALSE)
for (ext in c('pdf','svg','png')) {
  out <- paste0(out_prefix, '.split.v2.', ext)
  if (ext == 'pdf') pdf(out, width=4.2, height=4.9, useDingbats=FALSE)
  if (ext == 'svg') svg(out, width=4.2, height=4.9, pointsize=10)
  if (ext == 'png') png(out, width=4.2, height=4.9, units='in', res=500, type='cairo-png')
  plot_one()
  dev.off()
  cat('WROTE\t', out, '\n', sep='')
}