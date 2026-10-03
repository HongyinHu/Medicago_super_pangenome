args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
eigval <- args[2]
outprefix <- args[3]
d <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE)
varlab <- c("PC1", "PC2")
if (file.exists(eigval)) {
  ev <- scan(eigval, quiet=TRUE)
  if (length(ev) >= 2 && sum(ev) > 0) {
    varlab <- paste0(c("PC1", "PC2"), " (", round(ev[1:2] / sum(ev) * 100, 2), "%)")
  }
}
pdf(paste0(outprefix, ".pdf"), width=6.5, height=5.5)
plot(d$PC1, d$PC2, pch=21, bg="#2c7fb8", col="#08306b", xlab=varlab[1], ylab=varlab[2], main="PCA")
dev.off()
png(paste0(outprefix, ".png"), width=1800, height=1500, res=250)
plot(d$PC1, d$PC2, pch=21, bg="#2c7fb8", col="#08306b", xlab=varlab[1], ylab=varlab[2], main="PCA")
dev.off()
pdf(paste0(outprefix, ".labels.pdf"), width=8, height=7)
plot(d$PC1, d$PC2, pch=21, bg="#2c7fb8", col="#08306b", xlab=varlab[1], ylab=varlab[2], main="PCA with sample IDs")
text(d$PC1, d$PC2, labels=d$IID, cex=0.45, pos=3)
dev.off()
