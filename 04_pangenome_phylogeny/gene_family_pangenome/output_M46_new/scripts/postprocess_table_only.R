#!/usr/bin/env Rscript

# Post-process an 18-genome OrthoFinder run and draw manuscript figures.
# Category thresholds: Core=18, Softcore=15-17, Dispensable=2-14, Private=1.
# Unassigned genes are reported separately and excluded from all four categories.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop(
    "Usage: Rscript postprocess_pan_genome.R ",
    "<Orthogroups.GeneCount.tsv> <Orthogroups_UnassignedGenes.tsv> ",
    "<Statistics_PerSpecies.tsv> <output_directory>"
  )
}

gene_count_path <- normalizePath(args[[1]], mustWork = TRUE)
unassigned_path <- normalizePath(args[[2]], mustWork = TRUE)
statistics_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- normalizePath(args[[4]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

workspace_library <- file.path(
  getwd(), "figure_assembly_quality_matrix", "R_library"
)
if (dir.exists(workspace_library)) {
  .libPaths(c(workspace_library, .libPaths()))
}

  library(data.table)

gene_counts <- data.table::fread(gene_count_path)
stopifnot(all(c("Orthogroup", "Total") %in% names(gene_counts)))
sample_names <- setdiff(names(gene_counts), c("Orthogroup", "Total"))
stopifnot(length(sample_names) == 18L, !anyDuplicated(gene_counts$Orthogroup))

count_matrix <- as.matrix(gene_counts[, ..sample_names])
storage.mode(count_matrix) <- "integer"
presence_matrix <- count_matrix > 0L
genome_frequency <- rowSums(presence_matrix)

stopifnot(
  all(genome_frequency >= 1L & genome_frequency <= 18L),
  all(rowSums(count_matrix) == gene_counts$Total)
)

category_by_family <- fifelse(
  genome_frequency == 18L,
  "Core",
  fifelse(
    genome_frequency >= 15L,
    "Softcore",
    fifelse(genome_frequency >= 2L, "Dispensable", "Private")
  )
)
category_levels <- c("Private", "Dispensable", "Softcore", "Core")

frequency_table <- data.table(
  genome_frequency = 1:18,
  family_count = tabulate(genome_frequency, nbins = 18L)
)
frequency_table[
  ,
  category := factor(
    fifelse(
      genome_frequency == 18L,
      "Core",
      fifelse(
        genome_frequency >= 15L,
        "Softcore",
        fifelse(genome_frequency >= 2L, "Dispensable", "Private")
      )
    ),
    levels = category_levels
  )
]

category_summary <- frequency_table[
  ,
  .(family_count = sum(family_count)),
  by = category
]
category_summary[, percentage := family_count / sum(family_count) * 100]

private_matrix <- count_matrix[genome_frequency == 1L, , drop = FALSE]
private_by_genome <- data.table(
  sample = sample_names,
  private_family_count = as.integer(colSums(private_matrix > 0L)),
  private_gene_count = as.integer(colSums(private_matrix))
)

category_gene_matrix <- vapply(
  c("Core", "Softcore", "Dispensable", "Private"),
  function(category_name) {
    as.integer(colSums(count_matrix[category_by_family == category_name, , drop = FALSE]))
  },
  integer(length(sample_names))
)
rownames(category_gene_matrix) <- sample_names

unassigned <- data.table::fread(unassigned_path, fill = TRUE)
stopifnot("Orthogroup" %in% names(unassigned))
missing_unassigned_columns <- setdiff(sample_names, names(unassigned))
stopifnot(length(missing_unassigned_columns) == 0L)

count_gene_cells <- function(values) {
  values[is.na(values)] <- ""
  values <- trimws(values)
  sum(lengths(strsplit(values, ",", fixed = TRUE))[nzchar(values)])
}
unassigned_counts <- vapply(
  sample_names,
  function(sample_name) count_gene_cells(unassigned[[sample_name]]),
  integer(1)
)

# Statistics_PerSpecies.tsv contains two tabular sections separated by a blank
# line.  Only the first section has the per-species totals used below.
statistics_lines <- readLines(statistics_path, warn = FALSE)
first_blank_line <- which(!nzchar(trimws(statistics_lines)))[1L]
stopifnot(!is.na(first_blank_line), first_blank_line > 2L)
statistics <- data.table::fread(
  text = paste(statistics_lines[seq_len(first_blank_line - 1L)], collapse = "\n"),
  check.names = FALSE
)
metric_column <- names(statistics)[1]
extract_statistic <- function(metric) {
  row <- statistics[get(metric_column) == metric]
  stopifnot(nrow(row) == 1L, all(sample_names %in% names(row)))
  as.numeric(row[, ..sample_names])
}
statistics_total <- extract_statistic("Number of genes")
statistics_assigned <- extract_statistic("Number of genes in orthogroups")
statistics_unassigned <- extract_statistic("Number of unassigned genes")

assigned_counts <- colSums(count_matrix)
stopifnot(
  identical(as.integer(assigned_counts), as.integer(statistics_assigned)),
  identical(as.integer(unassigned_counts), as.integer(statistics_unassigned)),
  identical(
    as.integer(assigned_counts + unassigned_counts),
    as.integer(statistics_total)
  ),
  identical(
    as.integer(rowSums(category_gene_matrix)),
    as.integer(assigned_counts)
  )
)

mtr_sample <- intersect(c("Mtr_R108", "Mtr_500"), sample_names)
stopifnot(length(mtr_sample) == 1L)

species_map <- data.table(
  Entry = 1:18,
  sample = c(
    "Mca_474", "Msu_472", "Mla_454", "Mcr_468", "Mma_457", "Mse_461",
    "Msa_T2T", "Mra_436", "Mpr_410", "Med_482", "Mfi_M46", "Mlu_395",
    "Mpo_200", mtr_sample, "Mor_22", "Mar_100", "Mru_300", "Msa_zm4"
  ),
  Sample_ID = c(
    "Mca-474", "Msu-472", "Mla-454", "Mcr-468", "Mma-457", "Mse-461",
    "Msa-T2T", "Mra-436", "Mpr-410", "Med-482", "Mfi-M46", "Mlu-395",
    "Mpo-200", ifelse(mtr_sample == "Mtr_R108", "Mtr-R108", "Mtr-500"),
    "Mor-M22", "Mar-100", "Mru-300", "Msa-zm4"
  ),
  Species = c(
    "Medicago carstiensis",
    "Medicago suffruticosa",
    "Medicago lanigera",
    "Medicago cretacea",
    "Medicago marina",
    "Medicago secundiflora",
    "Medicago sativa subsp. caerulea",
    "Medicago radiata",
    "Medicago praecox",
    "Medicago edgeworthii",
    "Medicago fischeriana",
    "Medicago lupulina",
    "Medicago polymorpha",
    "Medicago truncatula R108",
    "Medicago orbicularis",
    "Medicago archiducis-nicolai",
    "Medicago ruthenica",
    "Medicago sativa cv. Zhongmu No. 4"
  )
)
stopifnot(setequal(species_map$sample, sample_names))

gene_category_table <- data.table(
  sample = sample_names,
  Core_genes_n = category_gene_matrix[, "Core"],
  Softcore_genes_n = category_gene_matrix[, "Softcore"],
  Dispensable_genes_n = category_gene_matrix[, "Dispensable"],
  Private_genes_n = category_gene_matrix[, "Private"],
  Genes_in_orthogroups_n = as.integer(assigned_counts),
  Unassigned_genes_n = as.integer(unassigned_counts),
  Total_genes_n = as.integer(assigned_counts + unassigned_counts)
)
gene_category_table <- merge(
  species_map,
  gene_category_table,
  by = "sample",
  sort = FALSE
)
setorder(gene_category_table, Entry)
setcolorder(
  gene_category_table,
  c(
    "Entry", "Sample_ID", "sample", "Species", "Core_genes_n",
    "Softcore_genes_n", "Dispensable_genes_n", "Private_genes_n",
    "Genes_in_orthogroups_n", "Unassigned_genes_n", "Total_genes_n"
  )
)

gene_category_percentages <- copy(gene_category_table)
for (column_name in c(
  "Core_genes_n", "Softcore_genes_n", "Dispensable_genes_n", "Private_genes_n"
)) {
  percentage_name <- sub("_genes_n$", "_genes_percent_of_assigned", column_name)
  gene_category_percentages[
    ,
    (percentage_name) := get(column_name) / Genes_in_orthogroups_n * 100
  ]
}

# Exact exhaustive pan/core accumulation across all 2^18 - 1 genome subsets.
n_species <- length(sample_names)
n_masks <- bitwShiftL(1L, n_species)
all_masks <- 0:(n_masks - 1L)
bit_values <- bitwShiftL(1L, 0:(n_species - 1L))
family_masks <- as.integer(rowSums(sweep(presence_matrix, 2L, bit_values, `*`)))
exact_pattern_counts <- as.numeric(tabulate(family_masks + 1L, nbins = n_masks))

subset_zeta <- exact_pattern_counts
superset_zeta <- exact_pattern_counts
for (bit_value in bit_values) {
  # which() already returns the 1-based array positions for mask + 1.
  upper_index <- which(bitwAnd(all_masks, bit_value) != 0L)
  subset_zeta[upper_index] <- subset_zeta[upper_index] +
    subset_zeta[upper_index - bit_value]

  lower_index <- which(bitwAnd(all_masks, bit_value) == 0L)
  superset_zeta[lower_index] <- superset_zeta[lower_index] +
    superset_zeta[lower_index + bit_value]
}

popcount <- integer(n_masks)
for (mask in 1:(n_masks - 1L)) {
  popcount[mask + 1L] <- popcount[bitwShiftR(mask, 1L) + 1L] +
    bitwAnd(mask, 1L)
}

subset_masks <- 1:(n_masks - 1L)
full_mask <- n_masks - 1L
complement_masks <- bitwXor(subset_masks, full_mask)
pan_core <- data.table(
  subset_mask = subset_masks,
  sample_number = popcount[subset_masks + 1L],
  pan_family_number = as.integer(
    nrow(gene_counts) - subset_zeta[complement_masks + 1L]
  ),
  core_family_number = as.integer(superset_zeta[subset_masks + 1L])
)

pan_core_summary <- pan_core[
  ,
  .(
    combinations = .N,
    pan_mean = mean(pan_family_number),
    pan_median = median(pan_family_number),
    pan_min = min(pan_family_number),
    pan_max = max(pan_family_number),
    core_mean = mean(core_family_number),
    core_median = median(core_family_number),
    core_min = min(core_family_number),
    core_max = max(core_family_number)
  ),
  by = sample_number
]
setorder(pan_core_summary, sample_number)
stopifnot(
  nrow(pan_core) == 2^18 - 1,
  pan_core_summary[sample_number == 18L, combinations] == 1L,
  pan_core_summary[sample_number == 18L, pan_mean] == nrow(gene_counts),
  pan_core_summary[sample_number == 18L, core_mean] ==
    frequency_table[genome_frequency == 18L, family_count]
)

data.table::fwrite(
  gene_category_table,
  file.path(output_dir, "Supplementary_Table_gene_category_counts.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  gene_category_percentages,
  file.path(output_dir, "Supplementary_Table_gene_category_counts_with_percentages.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  frequency_table,
  file.path(output_dir, "gene_family_frequency_distribution.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  category_summary,
  file.path(output_dir, "gene_family_category_summary.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  private_by_genome,
  file.path(output_dir, "private_families_by_genome.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  pan_core,
  file.path(output_dir, "pan_core_all_combinations.tsv.gz"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  pan_core_summary,
  file.path(output_dir, "pan_core_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

writeLines(
  c(
    paste0("orthogroups=", nrow(gene_counts)),
    paste0("core_families=", category_summary[category == "Core", family_count]),
    paste0("softcore_families=", category_summary[category == "Softcore", family_count]),
    paste0("dispensable_families=", category_summary[category == "Dispensable", family_count]),
    paste0("private_families=", category_summary[category == "Private", family_count]),
    paste0("assigned_genes=", sum(assigned_counts)),
    paste0("unassigned_genes=", sum(unassigned_counts)),
    paste0("total_genes=", sum(statistics_total)),
    "pan_core_method=exact exhaustive enumeration of all 262143 non-empty genome subsets",
    "unassigned_in_four_categories=0",
    "status=PASS"
  ),
  file.path(output_dir, "postprocess_summary.txt")
)

cat(readLines(file.path(output_dir, "postprocess_summary.txt")), sep = "\n")
