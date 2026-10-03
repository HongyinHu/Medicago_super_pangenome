args <- commandArgs(trailingOnly=TRUE)
out <- ifelse(length(args) >= 1, args[1], "path/to/project/N_4.pod_spiny/24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709")
label_file <- file.path(out, "data/Chr23997_strict_RBH_plus_outgroup.tip_labels.tsv")
label_map <- read.delim(label_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
plot_labels <- setNames(label_map$plot_label, label_map$tree_id)
read_newick <- function(file) paste(readLines(file, warn=FALSE), collapse="")
tokenize <- function(s) {
  m <- gregexpr("\\(|\\)|,|:|;|[^\\(\\),:;\\s]+", s, perl=TRUE)[[1]]
  regmatches(s, list(m))[[1]]
}
parse_newick <- function(s) {
  tok <- tokenize(s); i <- 1
  parse_sub <- function() {
    if (tok[i] == "(") {
      i <<- i + 1
      children <- list()
      repeat {
        children[[length(children)+1]] <- parse_sub()
        if (tok[i] == ",") { i <<- i + 1; next }
        if (tok[i] == ")") { i <<- i + 1; break }
      }
      label <- ""
      if (!tok[i] %in% c(":", ",", ")", ";")) { label <- tok[i]; i <<- i + 1 }
      bl <- 0
      if (tok[i] == ":") { i <<- i + 1; bl <- suppressWarnings(as.numeric(tok[i])); i <<- i + 1 }
      return(list(label=label, length=ifelse(is.na(bl),0,bl), children=children, x=NA, y=NA))
    } else {
      label <- tok[i]; i <<- i + 1
      bl <- 0
      if (tok[i] == ":") { i <<- i + 1; bl <- suppressWarnings(as.numeric(tok[i])); i <<- i + 1 }
      return(list(label=label, length=ifelse(is.na(bl),0,bl), children=list(), x=NA, y=NA))
    }
  }
  parse_sub()
}
assign_xy <- function(node, x0=0, next_y_env) {
  node$x <- x0 + node$length
  if (length(node$children)==0) {
    node$y <- next_y_env$y
    next_y_env$y <- next_y_env$y + 1
  } else {
    ys <- numeric(length(node$children))
    for (j in seq_along(node$children)) {
      node$children[[j]] <- assign_xy(node$children[[j]], node$x, next_y_env)
      ys[j] <- node$children[[j]]$y
    }
    node$y <- mean(ys)
  }
  node
}
collect_nodes <- function(node, parent=NULL) {
  rows <- list(list(label=node$label, x=node$x, y=node$y,
                    parent_x=if (is.null(parent)) NA else parent$x,
                    parent_y=if (is.null(parent)) NA else parent$y,
                    is_leaf=length(node$children)==0))
  for (ch in node$children) rows <- c(rows, collect_nodes(ch, node))
  rows
}
leaf_labels <- function(ids) {
  out <- ifelse(ids %in% names(plot_labels), plot_labels[ids], ids)
  unname(out)
}
plot_tree <- function(treefile, outfile_prefix, main) {
  e <- new.env(); e$y <- 1
  tr <- assign_xy(parse_newick(read_newick(treefile)), 0, e)
  rows <- collect_nodes(tr)
  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors=FALSE))
  df$x <- as.numeric(df$x); df$y <- as.numeric(df$y)
  df$parent_x <- as.numeric(df$parent_x); df$parent_y <- as.numeric(df$parent_y)
  leaves <- df[df$is_leaf==TRUE,]
  leaves$plot_label <- leaf_labels(leaves$label)
  xmax <- max(df$x, na.rm=TRUE)
  n <- nrow(leaves)
  draw <- function() {
    par(mar=c(3.1,1.0,2.0,8.9), xaxs="i", yaxs="i")
    plot(NA, xlim=c(0, xmax*1.55), ylim=c(0.5,n+0.5), axes=FALSE,
         xlab="substitutions/site", ylab="", main=main, cex.main=0.90)
    axis(1, cex.axis=0.72, lwd=0.6, lwd.ticks=0.6)
    for (k in seq_len(nrow(df))) {
      if (!is.na(df$parent_x[k])) segments(df$parent_x[k], df$y[k], df$x[k], df$y[k], lwd=0.8)
    }
    draw_v <- function(node) {
      if (length(node$children)>0) {
        ys <- vapply(node$children, function(z) z$y, numeric(1))
        segments(node$x, min(ys), node$x, max(ys), lwd=0.8)
        if (!is.null(node$label) && nzchar(node$label)) text(node$x, node$y+0.22, node$label, cex=0.40, col="#444444")
        for (ch in node$children) draw_v(ch)
      }
    }
    draw_v(tr)
    text(leaves$x + xmax*0.018, leaves$y, leaves$plot_label, adj=0, cex=0.58)
    box(bty="n")
  }
  grDevices::pdf(paste0(outfile_prefix,".pdf"), width=7.4, height=4.9, useDingbats=FALSE); draw(); dev.off()
  grDevices::png(paste0(outfile_prefix,".png"), width=7.4, height=4.9, units="in", res=600); draw(); dev.off()
  grDevices::svg(paste0(outfile_prefix,".svg"), width=7.4, height=4.9); draw(); dev.off()
}
plot_tree(file.path(out,"trees/Chr23997_pep_ML_rooted.treefile"), file.path(out,"figures/Chr23997_pep_ML_rooted"), "Chr23997 strict RBH + Chr23998 outgroup - PEP ML")
plot_tree(file.path(out,"trees/Chr23997_cds_ML_rooted.treefile"), file.path(out,"figures/Chr23997_cds_ML_rooted"), "Chr23997 strict RBH + Chr23998 outgroup - CDS ML")
plot_tree(file.path(out,"trees/Chr23997_codon12_ML_rooted.treefile"), file.path(out,"figures/Chr23997_codon12_ML_rooted"), "Chr23997 strict RBH + Chr23998 outgroup - codon12 ML")
