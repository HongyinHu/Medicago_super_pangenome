#!/usr/bin/env Rscript
# ED5a radial redraw in base R. All drawing and visual QA stay in R.
OUT <- "path/to/project/37.karyotype_reconstruction/output_ED5a_redo_20260929"
W <- file.path(OUT, "02_wgdi")
FIG <- file.path(OUT, "03_figure")
CONNECTED <- "--connected" %in% commandArgs(trailingOnly=TRUE)
SUFFIX <- if (CONNECTED) "_connected" else ""
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

AMK <- c(royalblue="AMK1", red="AMK2", "#99cc00"="AMK3",
         deepskyblue="AMK4", "#339966"="AMK5", "#ffcc00"="AMK6",
         fuchsia="AMK7", "#aa6e20"="AMK8")
AMK_COLOR <- c(AMK1="#3F5FA8", AMK2="#E3262A", AMK3="#7FB53A",
               AMK4="#29A3D8", AMK5="#1E7B4F", AMK6="#E0A11B",
               AMK7="#9B4FB0", AMK8="#A5561A")
BACKBONE <- "#E4E4E2"
EDGE <- "#8A8A88"
INK <- "#222222"
MUTED <- "#6B6B6B"
HILITE <- "#C0392B"

FIG_W <- 7.2
FIG_H <- 8.7
CX <- FIG_W / 2
CY <- FIG_H / 2
RX <- 3.08
RY <- 3.72
PW <- 0.70
PH <- 0.44
FS <- 5.5
LH <- 0.088
PAD <- 0.045

read_lens <- function(path) {
  x <- read.delim(path, header=FALSE, stringsAsFactors=FALSE,
                  col.names=c("chr", "bp", "genes"))
  x$chr <- as.character(x$chr)
  x$genes <- as.integer(x$genes)
  x
}

read_gff <- function(path) {
  x <- read.delim(path, header=FALSE, stringsAsFactors=FALSE)
  data.frame(chr=as.character(x[[1]]), gene=as.character(x[[2]]),
             start=as.integer(x[[3]]), end=as.integer(x[[4]]),
             order=as.integer(x[[6]]), stringsAsFactors=FALSE)
}

read_km <- function(path) {
  x <- read.delim(path, header=FALSE, stringsAsFactors=FALSE,
                  col.names=c("chr", "start", "end", "color", "cls"))
  x$chr <- as.character(x$chr)
  x$start <- as.integer(x$start)
  x$end <- as.integer(x$end)
  x$AMK <- unname(AMK[tolower(trimws(x$color))])
  if (anyNA(x$AMK)) stop("Unknown WGDI ancestor colour in ", path)
  x
}

read_ancestor <- function(path) read_km(path)

# Fill internal gaps at the midpoint between adjacent retained segments.
# Adjacent segments with the same AMK are merged after the fill.
connect_segments <- function(lens, segs) {
  out <- vector("list", 0)
  for (cc in as.character(lens$chr)) {
    z <- segs[segs$chr == cc, , drop=FALSE]
    if (!nrow(z)) next
    z <- z[order(z$start, z$end), , drop=FALSE]
    if (nrow(z) > 1L) for (i in seq_len(nrow(z)-1L)) {
      cut <- floor((z$end[i] + z$start[i+1L]) / 2)
      z$end[i] <- cut
      z$start[i+1L] <- cut + 1L
    }
    merged <- vector("list", 0)
    cur <- z[1, , drop=FALSE]
    if (nrow(z) > 1L) for (i in 2:nrow(z)) {
      if (identical(cur$AMK[1], z$AMK[i]) && cur$end[1] + 1L == z$start[i]) {
        cur$end[1] <- z$end[i]
      } else {
        merged[[length(merged)+1L]] <- cur
        cur <- z[i, , drop=FALSE]
      }
    }
    merged[[length(merged)+1L]] <- cur
    out[[length(out)+1L]] <- do.call(rbind, merged)
  }
  if (!length(out)) return(segs[0, , drop=FALSE])
  do.call(rbind, out)
}
label_lines <- function(full_name) {
  parts <- strsplit(full_name, " +")[[1]]
  genus <- parts[1]
  if (length(parts) == 3 && grepl("[0-9]", parts[3])) {
    return(list(c(genus, TRUE), c(paste(parts[2], parts[3]), FALSE)))
  }
  out <- list(c(genus, TRUE), c(parts[2], TRUE))
  if (length(parts) > 2) out[[length(out)+1]] <- c(paste(parts[3:length(parts)], collapse=" "), FALSE)
  out
}

draw_karyotype <- function(x0, y0, w, h, lens, segs, bar_frac=0.55) {
  slot <- w / nrow(lens)
  bw <- slot * bar_frac
  gmax <- max(lens$genes)
  for (i in seq_len(nrow(lens))) {
    c <- lens[i, ]
    bx <- x0 + (i - 1) * slot + (slot - bw) / 2
    bh <- h * c$genes / gmax
    top <- y0 + h
    rect(bx, top - bh, bx + bw, top, col=BACKBONE, border=NA)
    ss <- segs[segs$chr == c$chr, , drop=FALSE]
    if (nrow(ss)) for (j in seq_len(nrow(ss))) {
      s <- ss[j, ]
      y_hi <- top - bh * (s$start - 1) / c$genes
      y_lo <- top - bh * s$end / c$genes
      rect(bx, y_lo, bx + bw, y_hi, col=AMK_COLOR[s$AMK], border="white", lwd=0.2)
    }
    rect(bx, top - bh, bx + bw, top, col=NA, border=EDGE, lwd=0.3)
  }
}

prepare <- function() {
  sp <- read.delim(file.path(OUT, "species_table.tsv"), sep="\t", stringsAsFactors=FALSE)
  sp <- sp[order(sp$order), ]
  stopifnot(nrow(sp) == 19L)
  panels <- vector("list", nrow(sp))
  src <- vector("list", 0)
  boxes <- vector("list", nrow(sp))
  for (k in seq_len(nrow(sp))) {
    s <- sp[k, ]
    d <- file.path(W, s$label)
    lens <- read_lens(file.path(d, "in.lens1"))
    raw_segs <- read_km(file.path(d, "km_result.txt"))
    segs <- if (CONNECTED) connect_segments(lens, raw_segs) else raw_segs
    gff <- read_gff(file.path(d, "in.gff1"))
    stopifnot(nrow(lens) == as.integer(s$basic_x),
              length(unique(segs$chr)) == as.integer(s$basic_x))
    panels[[k]] <- list(s=s, lens=lens, segs=segs, gff=gff)
    theta <- pi / 2 - 2 * pi * (k - 1) / nrow(sp)
    px <- CX + RX * cos(theta)
    py <- CY + RY * sin(theta)
    lines <- label_lines(s$full_name)
    label_h <- (length(lines) + 1) * LH
    fx0 <- px - PW / 2 - PAD
    fy0 <- py - PH / 2 - PAD
    fw <- PW + 2 * PAD
    fh <- PH + 2 * PAD
    above <- sin(theta) > 0.3
    if (above) {
      boxes[[k]] <- c(label=s$label, x1=fx0, y1=fy0, x2=fx0+fw, y2=fy0+fh+0.03+label_h)
    } else {
      boxes[[k]] <- c(label=s$label, x1=fx0, y1=fy0-0.03-label_h, x2=fx0+fw, y2=fy0+fh)
    }
    chrom_genes <- setNames(lens$genes, lens$chr)
    key <- paste(gff$chr, gff$order, sep=":")
    for (j in seq_len(nrow(segs))) {
      g <- segs[j, ]
      a <- gff[match(paste(g$chr, g$start, sep=":"), key), ]
      b <- gff[match(paste(g$chr, g$end, sep=":"), key), ]
      if (nrow(a) != 1L || nrow(b) != 1L) stop("GFF index lookup failed for ", s$label)
      src[[length(src)+1L]] <- data.frame(
        order=as.integer(s$order), label=s$label, species=s$full_name,
        basic_x=as.integer(s$basic_x), analysis_dir=d, chromosome=g$chr,
        AMK=g$AMK, display_color=unname(AMK_COLOR[g$AMK]),
        start_gene_index=g$start, end_gene_index=g$end,
        segment_gene_count=g$end-g$start+1L,
        chromosome_gene_count=as.integer(chrom_genes[[g$chr]]),
        start_gene_id=a$gene, start_bp=a$start,
        end_gene_id=b$gene, end_bp=b$end, stringsAsFactors=FALSE)
    }
  }
  list(sp=sp, panels=panels, src=do.call(rbind, src), boxes=boxes)
}

draw_figure <- function(dat) {
  par(mar=c(0,0,0,0), oma=c(0,0,0,0), xpd=NA, family="sans")
  plot.new()
  plot.window(xlim=c(0, FIG_W), ylim=c(0, FIG_H), asp=1)
  # The central AMK disc and its 8 ancestral chromosomes
  symbols(CX, CY, circles=1.0, inches=FALSE, add=TRUE,
          bg="#F6F6F4", fg="#D5D5D2", lwd=0.6)
  anc <- read_ancestor(file.path(W, "Mfis", "in.ancestor"))
  alens <- read_lens(file.path(W, "Mfis", "in.lens2"))
  aw <- 1.24; ah <- 0.74
  draw_karyotype(CX-aw/2, CY-ah/2+0.06, aw, ah, alens, anc, bar_frac=0.6)
  text(CX-aw/2 + ((seq_len(8)-0.5)*aw/8), CY-ah/2+0.01,
       labels=seq_len(8), cex=0.65, col=INK, adj=c(0.5,1))
  text(CX, CY+ah/2+0.16, "AMK", cex=1.0, font=2, col=INK, adj=c(0.5,0))
  text(CX-0.02, CY-ah/2-0.16, "Ancestral", cex=0.55, col=MUTED, adj=c(1,1))
  text(CX+0.02, CY-ah/2-0.16, "Medicago", cex=0.55, col=MUTED, font=3, adj=c(0,1))
  text(CX, CY-ah/2-0.25, "karyotype (x = 8)", cex=0.55, col=MUTED, adj=c(0.5,1))

  for (k in seq_along(dat$panels)) {
    p <- dat$panels[[k]]
    s <- p$s; lens <- p$lens; segs <- p$segs
    theta <- pi/2 - 2*pi*(k-1)/nrow(dat$sp)
    px <- CX + RX*cos(theta); py <- CY + RY*sin(theta)
    fx0 <- px-PW/2-PAD; fy0 <- py-PH/2-PAD; fw <- PW+2*PAD; fh <- PH+2*PAD
    above <- sin(theta) > 0.3
    phi <- atan2(py-CY, px-CX)
    sx <- CX + 1.02*cos(phi); sy <- CY + 1.02*sin(phi)
    tx <- (fw/2+0.03)/max(abs(cos(phi)), 1e-9)
    ty <- (fh/2+0.03)/max(abs(sin(phi)), 1e-9)
    tt <- min(tx, ty)
    segments(sx, sy, px-tt*cos(phi), py-tt*sin(phi), col="#C4C4C1", lwd=0.5)
    x7 <- as.integer(s$basic_x) == 7L
    rect(fx0, fy0, fx0+fw, fy0+fh, col="white",
         border=if (x7) HILITE else "#CFCFCC", lwd=if (x7) 0.8 else 0.5,
         lty=if (x7) 2 else 1)
    draw_karyotype(px-PW/2, py-PH/2, PW, PH, lens, segs)
    lines <- label_lines(s$full_name)
    label_h <- (length(lines)+1)*LH
    y <- if (above) fy0+fh+0.03+label_h else fy0-0.03
    for (z in lines) {
      text(px, y, z[1], cex=FS/10, col=INK, font=if (as.logical(z[2])) 3 else 1, adj=c(0.5,1))
      y <- y-LH
    }
    text(px, y, paste0("x = ", s$basic_x), cex=FS/10,
         col=if (x7) HILITE else MUTED, font=if (x7) 2 else 1, adj=c(0.5,1))
  }
  problems <- character()
  for (i in seq_len(length(dat$boxes)-1L)) for (j in (i+1L):length(dat$boxes)) {
    a <- dat$boxes[[i]]; b <- dat$boxes[[j]]
    if (a["x1"] < b["x2"] && b["x1"] < a["x2"] && a["y1"] < b["y2"] && b["y1"] < a["y2"])
      problems <- c(problems, paste("overlap", a["label"], "/", b["label"]))
  }
  for (b in dat$boxes) if (b["x1"] < 0.05 || b["y1"] < 0.05 || b["x2"] > FIG_W-0.05 || b["y2"] > FIG_H-0.05)
    problems <- c(problems, paste("off-canvas", b["label"]))
  invisible(problems)
}

dat <- prepare()
pdf(file.path(FIG, paste0("ED5a_radial", SUFFIX, ".pdf")), width=FIG_W, height=FIG_H, useDingbats=FALSE)
problems <- draw_figure(dat)
dev.off()
png(file.path(FIG, paste0("ED5a_radial", SUFFIX, ".png")), width=FIG_W, height=FIG_H, units="in", res=600, type="cairo-png")
problems2 <- draw_figure(dat)
dev.off()
svg(file.path(FIG, paste0("ED5a_radial", SUFFIX, ".svg")), width=FIG_W, height=FIG_H, onefile=TRUE, family="sans")
problems3 <- draw_figure(dat)
dev.off()
write.table(dat$src, file.path(FIG, paste0("ED5a_source_data", SUFFIX, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE)
if (length(problems) || length(problems2) || length(problems3))
  stop(paste(unique(c(problems, problems2, problems3)), collapse="; "))
cat("R redraw complete\nsegments:", nrow(dat$src), "\nlayout QC: PASS\nMfis chromosomes:", nrow(dat$panels[[which(dat$sp$label == "Mfis")]]$lens), "\n")



