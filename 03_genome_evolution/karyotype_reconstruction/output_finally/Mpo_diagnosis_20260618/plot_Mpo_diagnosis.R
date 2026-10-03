read_matrix <- function(file) {
  x <- read.delim(file, check.names = FALSE)
  rn <- x[[1]]
  m <- as.matrix(x[,-1, drop=FALSE])
  storage.mode(m) <- 'numeric'
  rownames(m) <- paste0('Mpo', rn)
  m
}
anc <- read_matrix('Mpo_km_ancestor_matrix_fraction.tsv')
r108 <- read_matrix('Mpo_vs_R108_besthit_matrix_fraction.tsv')
anc_counts <- read.delim('Mpo_km_ancestor_matrix_counts.tsv', check.names = FALSE)
r108_counts <- read.delim('Mpo_vs_R108_besthit_matrix_counts.tsv', check.names = FALSE)
anc_cols <- c('#4E79A7','#E15759','#59A14F','#76B7B2','#2F6F4E','#F1CE63','#B07AA1','#9C755F')
r108_cols <- anc_cols
plot_heat <- function(m, title, filebase) {
  pal <- colorRampPalette(c('#F8FAFC','#DCEBEE','#80B8B3','#2A7F78','#124C45'))(100)
  pdf(paste0(filebase,'.pdf'), width=7.2, height=4.8, useDingbats=FALSE)
  par(mar=c(5,5,3,5), family='serif')
  image(t(m[nrow(m):1,]), axes=FALSE, col=pal, zlim=c(0,1), xlab='', ylab='', main=title)
  axis(1, at=seq(0,1,length.out=ncol(m)), labels=colnames(m), las=2, cex.axis=0.8)
  axis(2, at=seq(0,1,length.out=nrow(m)), labels=rev(rownames(m)), las=1, cex.axis=0.85)
  for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
    val <- m[nrow(m)-i+1,j]
    if (val >= 0.05) text((j-1)/(ncol(m)-1), (i-1)/(nrow(m)-1), sprintf('%.0f', val*100), cex=0.65, col=ifelse(val>0.45,'white','#27313A'))
  }
  legend_x <- 1.08
  rect(legend_x, seq(0,0.88,length.out=10), legend_x+0.035, seq(0.08,0.96,length.out=10), col=pal[seq(1,100,length.out=10)], border=NA, xpd=NA)
  text(legend_x+0.08, c(0,0.5,0.96), c('0','50%','100%'), xpd=NA, cex=0.75, adj=0)
  dev.off()
  png(paste0(filebase,'.png'), width=2200, height=1450, res=300)
  par(mar=c(5,5,3,5), family='serif')
  image(t(m[nrow(m):1,]), axes=FALSE, col=pal, zlim=c(0,1), xlab='', ylab='', main=title)
  axis(1, at=seq(0,1,length.out=ncol(m)), labels=colnames(m), las=2, cex.axis=0.8)
  axis(2, at=seq(0,1,length.out=nrow(m)), labels=rev(rownames(m)), las=1, cex.axis=0.85)
  for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
    val <- m[nrow(m)-i+1,j]
    if (val >= 0.05) text((j-1)/(ncol(m)-1), (i-1)/(nrow(m)-1), sprintf('%.0f', val*100), cex=0.65, col=ifelse(val>0.45,'white','#27313A'))
  }
  legend_x <- 1.08
  rect(legend_x, seq(0,0.88,length.out=10), legend_x+0.035, seq(0.08,0.96,length.out=10), col=pal[seq(1,100,length.out=10)], border=NA, xpd=NA)
  text(legend_x+0.08, c(0,0.5,0.96), c('0','50%','100%'), xpd=NA, cex=0.75, adj=0)
  dev.off()
}
plot_stack <- function(m, cols, title, filebase) {
  pdf(paste0(filebase,'.pdf'), width=7.5, height=4.8, useDingbats=FALSE)
  par(mar=c(4.5,4.7,3,6.5), family='serif', xpd=FALSE)
  bp <- barplot(t(m), horiz=TRUE, col=cols, border='white', xlim=c(0,1), las=1, xlab='Fraction of genes', main=title, cex.names=0.85)
  abline(v=c(0.25,0.5,0.75), col='#E5E7EB', lty=3)
  legend('right', inset=c(-0.25,0), legend=colnames(m), fill=cols, border=NA, bty='n', cex=0.78, xpd=TRUE)
  dev.off()
  png(paste0(filebase,'.png'), width=2300, height=1450, res=300)
  par(mar=c(4.5,4.7,3,6.5), family='serif', xpd=FALSE)
  bp <- barplot(t(m), horiz=TRUE, col=cols, border='white', xlim=c(0,1), las=1, xlab='Fraction of genes', main=title, cex.names=0.85)
  abline(v=c(0.25,0.5,0.75), col='#E5E7EB', lty=3)
  legend('right', inset=c(-0.25,0), legend=colnames(m), fill=cols, border=NA, bty='n', cex=0.78, xpd=TRUE)
  dev.off()
}
plot_heat(anc, 'Mpo WGDI karyotype assignment vs aak ancestors (%)', 'Mpo_WGDI_aak_ancestor_heatmap')
plot_heat(r108, 'Mpo vs R108 best-hit chromosome matrix (%)', 'Mpo_R108_besthit_heatmap')
plot_stack(anc, anc_cols, 'Mpo WGDI karyotype assignment', 'Mpo_WGDI_aak_ancestor_stackedbar')
plot_stack(r108, r108_cols, 'Mpo vs R108 best-hit chromosome assignment', 'Mpo_R108_besthit_stackedbar')
