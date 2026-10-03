#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript qa_pan_genome_outputs.R <postprocess_output_directory>")
}

output_dir <- normalizePath(args[[1]], mustWork = TRUE)
required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

read_tsv <- function(name) {
  path <- file.path(output_dir, name)
  if (grepl("\\.gz$", path)) {
    return(data.table::as.data.table(utils::read.delim(
      gzfile(path),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )))
  }
  data.table::fread(path)
}

gene_table <- read_tsv("Supplementary_Table_gene_category_counts.tsv")
frequency <- read_tsv("gene_family_frequency_distribution.tsv")
categories <- read_tsv("gene_family_category_summary.tsv")
private <- read_tsv("private_families_by_genome.tsv")
pan_summary <- read_tsv("pan_core_summary.tsv")
pan_all <- read_tsv("pan_core_all_combinations.tsv.gz")

stopifnot(
  nrow(gene_table) == 18L,
  !anyDuplicated(gene_table$Entry),
  !anyDuplicated(gene_table$sample),
  all(rowSums(gene_table[, .(
    Core_genes_n,
    Softcore_genes_n,
    Dispensable_genes_n,
    Private_genes_n
  )]) == gene_table$Genes_in_orthogroups_n),
  all(
    gene_table$Genes_in_orthogroups_n + gene_table$Unassigned_genes_n ==
      gene_table$Total_genes_n
  ),
  all(gene_table$Unassigned_genes_n >= 0L),
  nrow(frequency) == 18L,
  sum(frequency$family_count) == sum(categories$family_count),
  frequency[genome_frequency == 1L, family_count] ==
    sum(private$private_family_count),
  frequency[genome_frequency == 18L, family_count] ==
    categories[as.character(category) == "Core", family_count],
  nrow(private) == 18L,
  nrow(pan_all) == 2^18 - 1,
  nrow(pan_summary) == 18L,
  pan_summary[sample_number == 18L, combinations] == 1L,
  pan_summary[sample_number == 18L, pan_mean] == sum(frequency$family_count),
  pan_summary[sample_number == 18L, core_mean] ==
    frequency[genome_frequency == 18L, family_count]
)

plot_stubs <- c(
  "Panel_a_pan_core_accumulation",
  "Panel_b_gene_family_frequency",
  "Pan_genome_two_panel_figure"
)
plot_extensions <- c("svg", "pdf", "png", "tiff")
plot_files <- as.vector(outer(
  plot_stubs,
  plot_extensions,
  function(stub, ext) file.path(output_dir, paste0(stub, ".", ext))
))
plot_info <- file.info(plot_files)
stopifnot(
  all(file.exists(plot_files)),
  all(!is.na(plot_info$size)),
  all(plot_info$size > 1000)
)

read_signature <- function(path, n) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  readBin(connection, what = "raw", n = n)
}

png_dimensions <- function(path) {
  bytes <- read_signature(path, 24L)
  stopifnot(identical(as.integer(bytes[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
  decode_uint32_be <- function(x) {
    sum(as.numeric(as.integer(x)) * 256^(3:0))
  }
  c(
    width = decode_uint32_be(bytes[17:20]),
    height = decode_uint32_be(bytes[21:24])
  )
}

expected_png <- list(
  Panel_a_pan_core_accumulation = c(width = 2126, height = 1984),
  Panel_b_gene_family_frequency = c(width = 2126, height = 1984),
  Pan_genome_two_panel_figure = c(width = 4323, height = 2031)
)
for (stub in names(expected_png)) {
  observed <- png_dimensions(file.path(output_dir, paste0(stub, ".png")))
  stopifnot(all(abs(observed - expected_png[[stub]]) <= 1L))
}

for (stub in plot_stubs) {
  svg_text <- paste(
    readLines(file.path(output_dir, paste0(stub, ".svg")), warn = FALSE),
    collapse = "\n"
  )
  stopifnot(
    grepl("<text", svg_text, fixed = TRUE),
    grepl("Arial", svg_text, fixed = TRUE)
  )
  pdf_signature <- rawToChar(
    read_signature(file.path(output_dir, paste0(stub, ".pdf")), 4L)
  )
  stopifnot(identical(pdf_signature, "%PDF"))
}

summary_lines <- readLines(file.path(output_dir, "postprocess_summary.txt"))
stopifnot("status=PASS" %in% summary_lines)

qa_lines <- c(
  "status=PASS",
  paste0("genomes=", nrow(gene_table)),
  paste0("orthogroups=", sum(frequency$family_count)),
  paste0("core_families=", categories[as.character(category) == "Core", family_count]),
  paste0("softcore_families=", categories[as.character(category) == "Softcore", family_count]),
  paste0("dispensable_families=", categories[as.character(category) == "Dispensable", family_count]),
  paste0("private_families=", categories[as.character(category) == "Private", family_count]),
  paste0("assigned_genes=", sum(gene_table$Genes_in_orthogroups_n)),
  paste0("unassigned_genes=", sum(gene_table$Unassigned_genes_n)),
  paste0("total_genes=", sum(gene_table$Total_genes_n)),
  "category_closure=PASS",
  "orthofinder_total_closure=PASS",
  "pan_core_combinations=262143",
  "frequency_endpoints=PASS",
  "figure_bundle=SVG,PDF,PNG,TIFF",
  "png_resolution=600_dpi_at_90_or_183_mm",
  "svg_editable_text=PASS",
  "pdf_signature=PASS",
  "visual_inspection=REQUIRED_MANUALLY"
)
writeLines(qa_lines, file.path(output_dir, "QA_report.txt"))
cat(qa_lines, sep = "\n")
cat("\n")
