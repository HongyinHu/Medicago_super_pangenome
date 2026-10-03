#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 9) {
  stop("Usage: run_gmmat_gwas_maf.R bfile model.tsv grm.rel grm.rel.id out_file label n_pcs ncores min_maf")
}

bfile <- args[[1]]
model_file <- args[[2]]
grm_file <- args[[3]]
grm_id_file <- args[[4]]
out_file <- args[[5]]
label <- args[[6]]
n_pcs <- as.integer(args[[7]])
ncores <- as.integer(args[[8]])
min_maf <- as.numeric(args[[9]])
if (!is.finite(min_maf) || min_maf <= 0 || min_maf >= 0.5) {
  stop("min_maf must be strictly between 0 and 0.5")
}

root <- Sys.getenv("GWAS_ROOT")
if (!nzchar(root)) stop("GWAS_ROOT is not set")
library(GMMAT, lib.loc=file.path(root, "software", "Rlib"))

dir.create(dirname(out_file), recursive=TRUE, showWarnings=FALSE)
meta <- read.delim(model_file, check.names=FALSE, stringsAsFactors=FALSE)
fam <- read.table(paste0(bfile, ".fam"), stringsAsFactors=FALSE)
if (!identical(as.character(fam[[2]]), as.character(meta$id))) {
  stop("Model metadata order does not match PLINK FAM order")
}

grm_ids <- read.table(grm_id_file, stringsAsFactors=FALSE)
if (ncol(grm_ids) < 2) stop("Malformed GRM ID file")
grm_iids <- as.character(grm_ids[[2]])
grm_values <- scan(grm_file, quiet=TRUE)
n_grm <- length(grm_iids)
if (length(grm_values) != n_grm * n_grm) {
  stop("GRM dimensions do not match GRM ID file")
}
K <- matrix(grm_values, nrow=n_grm, ncol=n_grm, byrow=TRUE)
rownames(K) <- grm_iids
colnames(K) <- grm_iids
if (!all(meta$id %in% grm_iids)) stop("Some model samples are absent from the SNP GRM")
K <- K[meta$id, meta$id, drop=FALSE]
K <- (K + t(K)) / 2
if (any(!is.finite(K))) stop("GRM contains non-finite values")

eigen_min <- min(eigen(K, symmetric=TRUE, only.values=TRUE)$values)
ridge <- 0
if (eigen_min < 1e-8) {
  ridge <- 1e-8 - eigen_min
  diag(K) <- diag(K) + ridge
}

pc_names <- paste0("PC", seq_len(n_pcs))
if (!all(pc_names %in% names(meta))) stop("Missing requested PC covariates")
fixed <- as.formula(paste("y ~", paste(pc_names, collapse=" + ")))
null_model <- glmmkin(
  fixed,
  data=meta,
  kins=K,
  id="id",
  family=binomial(link="logit"),
  method="REML",
  method.optim="AI",
  maxiter=500,
  tol=1e-5,
  verbose=TRUE
)
saveRDS(null_model, paste0(out_file, ".null.rds"))

glmm.score(
  null_model,
  infile=bfile,
  outfile=out_file,
  MAF.range=c(min_maf, 0.5),
  miss.cutoff=0.20,
  missing.method="impute2mean",
  nperbatch=100,
  ncores=1,
  verbose=TRUE
)
if (!file.exists(out_file) || file.info(out_file)$size == 0) {
  stop("GMMAT score output was not created")
}

audit <- data.frame(
  metric=c("label", "samples", "spiral_cases", "nonspiral_controls", "PCs", "min_maf", "GRM_min_eigen_before_ridge", "GRM_ridge_added", "requested_cores", "score_scan_cores", "theta"),
  value=c(label, nrow(meta), sum(meta$y == 1), sum(meta$y == 0), n_pcs,
          format(min_maf, digits=12), format(eigen_min, digits=12), format(ridge, digits=12), ncores, 1,
          paste(null_model$theta, collapse=","))
)
write.table(audit, paste0(out_file, ".model_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
