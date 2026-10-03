args <- commandArgs(trailingOnly=TRUE)
fam <- args[1]
outdir <- args[2]
kmin <- as.integer(args[3]); kmax <- as.integer(args[4])
f <- read.table(fam, stringsAsFactors=FALSE)
ids <- f[,2]
cv <- data.frame(K=integer(), CV=numeric())
for (K in kmin:kmax) {
  qfile <- file.path(outdir, paste0("sativa182.admix.", K, ".Q"))
  if (!file.exists(qfile)) next
  q <- as.matrix(read.table(qfile))
  ord <- order(max.col(q, ties.method="first"), ids)
  pdf(file.path(outdir, paste0("admixture.K", K, ".pdf")), width=10, height=3.8)
  par(mar=c(5,4,2,1))
  barplot(t(q[ord,,drop=FALSE]), col=rainbow(K), border=NA, space=0, names.arg=ids[ord], las=2, cex.names=0.32, ylab="Ancestry proportion", main=paste0("ADMIXTURE K=", K))
  dev.off()
  png(file.path(outdir, paste0("admixture.K", K, ".png")), width=2600, height=900, res=220)
  par(mar=c(5,4,2,1))
  barplot(t(q[ord,,drop=FALSE]), col=rainbow(K), border=NA, space=0, names.arg=ids[ord], las=2, cex.names=0.32, ylab="Ancestry proportion", main=paste0("ADMIXTURE K=", K))
  dev.off()
}
logs <- list.files(outdir, pattern="^K[0-9]+\\.log$", full.names=TRUE)
if (length(logs)) {
  for (lf in logs) {
    txt <- readLines(lf, warn=FALSE)
    hit <- grep("CV error", txt, value=TRUE)
    if (length(hit)) {
      K <- as.integer(sub("^K([0-9]+)\\.log$", "\\1", basename(lf)))
      val <- as.numeric(sub(".*: ", "", hit[length(hit)]))
      cv <- rbind(cv, data.frame(K=K, CV=val))
    }
  }
  if (nrow(cv)) {
    cv <- cv[order(cv$K),]
    write.table(cv, file=file.path(outdir, "admixture_cv_errors.clean.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    pdf(file.path(outdir, "admixture_cv_errors.pdf"), width=5.5, height=4.5)
    plot(cv$K, cv$CV, type="b", pch=19, xlab="K", ylab="Cross-validation error", main="ADMIXTURE CV error")
    dev.off()
    png(file.path(outdir, "admixture_cv_errors.png"), width=1300, height=1000, res=220)
    plot(cv$K, cv$CV, type="b", pch=19, xlab="K", ylab="Cross-validation error", main="ADMIXTURE CV error")
    dev.off()
  }
}
