args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(paste("usage: run_candidate_gmmat.R <private_lib> <strict.tsv> <sensitivity.tsv>",
             "<sample_order.keep> <pc.txt> <kinship.txt> <out.tsv> <log.tsv>"))
}

private_lib <- normalizePath(args[[1]], mustWork = FALSE)
strict_path <- args[[2]]
sensitivity_path <- args[[3]]
keep_path <- args[[4]]
pc_path <- args[[5]]
kinship_path <- args[[6]]
out_path <- args[[7]]
log_path <- args[[8]]

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(private_lib, .libPaths()))

empty_result <- function(callset, status, n_total = NA_integer_, n_called = NA_integer_) {
  data.frame(callset = callset, model = "GMMAT_logistic_PC1_PC5_GRM", n_total = n_total,
             n_called = n_called, beta = NA_real_, se = NA_real_, odds_ratio = NA_real_,
             ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_, status = status,
             stringsAsFactors = FALSE)
}

if (!requireNamespace("GMMAT", quietly = TRUE)) {
  result <- rbind(empty_result("strict", "unavailable: GMMAT is not installed"),
                  empty_result("sensitivity", "unavailable: GMMAT is not installed"))
  write.table(result, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("GMMAT not installed", log_path)
  quit(status = 0L)
}

keep <- read.table(keep_path, sep = "\t", stringsAsFactors = FALSE)
sample_order <- keep[[2]]
pcs <- read.table(pc_path, header = FALSE, stringsAsFactors = FALSE)
if (nrow(pcs) != length(sample_order) || ncol(pcs) < 6L) stop("PC file does not match keep file")
pc_df <- data.frame(sample = sample_order, PC1 = pcs[[2]], PC2 = pcs[[3]], PC3 = pcs[[4]],
                    PC4 = pcs[[5]], PC5 = pcs[[6]], stringsAsFactors = FALSE)
kinship <- as.matrix(read.table(kinship_path, header = FALSE, check.names = FALSE))
if (nrow(kinship) != length(sample_order) || ncol(kinship) != length(sample_order)) stop("kinship dimensions do not match keep file")
rownames(kinship) <- sample_order
colnames(kinship) <- sample_order

fit_one <- function(path, callset) {
  typed <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE, comment.char = "")
  typed$phenotype <- as.numeric(typed$phenotype)
  typed$genotype_num <- suppressWarnings(as.numeric(typed$genotype))
  called <- typed[!is.na(typed$genotype_num), c("sample", "phenotype", "genotype_num")]
  df <- merge(called, pc_df, by = "sample", sort = FALSE)
  df <- df[match(called$sample, df$sample), ]
  if (nrow(df) != nrow(called)) return(empty_result(callset, "not_estimable: PC join failed", nrow(typed), nrow(called)))
  if (length(unique(df$phenotype)) < 2L || length(unique(df$genotype_num)) < 2L) {
    return(empty_result(callset, "not_estimable: no phenotype or genotype variation", nrow(typed), nrow(called)))
  }
  k_sub <- kinship[df$sample, df$sample, drop = FALSE]
  rownames(k_sub) <- df$sample
  colnames(k_sub) <- df$sample
  tryCatch({
    fit_warnings <- character()
    fit <- withCallingHandlers(
      GMMAT::glmmkin(phenotype ~ genotype_num + PC1 + PC2 + PC3 + PC4 + PC5,
                      data = df, kins = k_sub, id = "sample", family = binomial(link = "logit")),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    coef_table <- summary(fit)$coefficients
    if (is.null(coef_table) || !("genotype_num" %in% rownames(coef_table))) stop("genotype coefficient absent")
    estimate <- coef_table["genotype_num", 1]
    se <- coef_table["genotype_num", 2]
    p_column <- grep("Pr", colnames(coef_table), value = TRUE)[1]
    pvalue <- if (is.na(p_column)) NA_real_ else coef_table["genotype_num", p_column]
    unstable <- !is.finite(estimate) || !is.finite(se) || abs(estimate) > 20 || se > 10 || length(fit_warnings) > 0
    if (unstable) {
      return(data.frame(callset = callset, model = "GMMAT_logistic_PC1_PC5_GRM", n_total = nrow(typed),
                        n_called = nrow(called), beta = estimate, se = se, odds_ratio = NA_real_,
                        ci_low = NA_real_, ci_high = NA_real_, pvalue = NA_real_,
                        status = paste0("unstable: ", paste(unique(fit_warnings), collapse = "; ")),
                        stringsAsFactors = FALSE))
    }
    data.frame(callset = callset, model = "GMMAT_logistic_PC1_PC5_GRM", n_total = nrow(typed),
               n_called = nrow(called), beta = estimate, se = se, odds_ratio = exp(estimate),
               ci_low = exp(estimate - 1.96 * se), ci_high = exp(estimate + 1.96 * se),
               pvalue = pvalue, status = "ok", stringsAsFactors = FALSE)
  }, error = function(err) {
    empty_result(callset, paste0("not_estimable: ", conditionMessage(err)), nrow(typed), nrow(called))
  })
}

sink(log_path, split = TRUE)
cat("GMMAT version:", as.character(packageVersion("GMMAT")), "\n")
result <- rbind(fit_one(strict_path, "strict"), fit_one(sensitivity_path, "sensitivity"))
print(result)
sink()
write.table(result, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
