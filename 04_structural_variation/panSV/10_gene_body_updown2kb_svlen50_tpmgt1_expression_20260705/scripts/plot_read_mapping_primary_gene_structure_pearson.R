suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

work <- 'path/to/project/N_3.call_SV/10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705'
tabdir <- file.path(work, 'tables')
figdir <- file.path(work, 'figures')
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

targets <- c('genome_410','genome_461','genome_472','genome_474',
             'genome_Mar','genome_Mru','genome_R108','genome_ZM4')
feature_levels <- c('exon', 'intron', '2kb downstream', '2kb upstream')
priority_levels <- c('exon', 'intron', '2kb upstream', '2kb downstream')

pearson_distance <- function(x, y) {
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  1 - suppressWarnings(cor(x, y, method = 'pearson'))
}

format_p <- function(p) {
  if (is.na(p)) return('NA')
  if (p < 2.2e-16) return('P < 2.2e-16')
  if (p < 0.001) return(paste0('P = ', formatC(p, format = 'e', digits = 1)))
  paste0('P = ', signif(p, 2))
}

expr <- read.delim(file.path(tabdir, 'one_to_one_orthologs_refTarget_TPMgt1.tsv'),
                   check.names = FALSE, stringsAsFactors = FALSE)
events <- read.delim(file.path(tabdir, 'SVLENgt50_geneBody_updown2kb_candidate_events_with_PAV.tsv'),
                     check.names = FALSE, stringsAsFactors = FALSE)

events$feature_class <- NA_character_
events$feature_class[events$category4 == 'Gene body' & events$gene_body_detail == 'exon_overlap'] <- 'exon'
events$feature_class[events$category4 == 'Gene body' & events$gene_body_detail == 'intron_only'] <- 'intron'
events$feature_class[events$category4 == 'Downstream'] <- '2kb downstream'
events$feature_class[events$category4 == 'Upstream'] <- '2kb upstream'
events <- events[!is.na(events$feature_class), ]
events$feature_priority <- match(events$feature_class, priority_levels)

gene_class_raw <- list()
k <- 1
for (tg in targets) {
  sub <- events[events[[tg]] == 1, c('SV_ID', 'nearest_gene_id', 'feature_class', 'feature_priority',
                                    'SVTYPE', 'SVLEN_abs')]
  if (nrow(sub) == 0) next
  sub$target_group <- tg
  gene_class_raw[[k]] <- sub
  k <- k + 1
}
gene_class_raw <- bind_rows(gene_class_raw)

# Mutually exclusive assignment: if one gene has several SV categories in one target,
# use exon > intron > upstream 2 kb > downstream 2 kb.
gene_class <- gene_class_raw %>%
  group_by(target_group, gene_id = nearest_gene_id) %>%
  arrange(feature_priority, desc(SVLEN_abs), .by_group = TRUE) %>%
  summarise(feature_class = first(feature_class),
            representative_SV_ID = first(SV_ID),
            representative_SVTYPE = first(SVTYPE),
            representative_SVLEN_abs = first(SVLEN_abs),
            n_candidate_events = n(),
            .groups = 'drop')

write.table(gene_class_raw,
            file.path(tabdir, 'read_mapping_primary_gene_structure_all_geneSV_pairs.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)
write.table(gene_class,
            file.path(tabdir, 'read_mapping_primary_gene_structure_exclusive_gene_category.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

plotdat <- expr %>%
  inner_join(gene_class, by = c('target_group', 'gene_id')) %>%
  mutate(target_group = factor(target_group, levels = targets),
         feature_class = factor(feature_class, levels = feature_levels),
         ref_log2 = log2(mean_TPM_ref + 1),
         target_log2 = log2(mean_TPM_target + 1))

category_counts <- plotdat %>%
  count(target_group, feature_class, name = 'n_genes') %>%
  arrange(target_group, feature_class)
write.table(category_counts,
            file.path(tabdir, 'read_mapping_primary_gene_structure_TPMgt1_gene_counts.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

set.seed(1)
B <- 10000
boots <- list()
sums <- list()
k <- 1
for (tg in targets) {
  for (fc in feature_levels) {
    sub <- plotdat %>% filter(target_group == tg, feature_class == fc)
    n <- nrow(sub)
    observed <- pearson_distance(sub$ref_log2, sub$target_log2)
    out <- rep(NA_real_, B)
    if (n >= 2) {
      for (i in seq_len(B)) {
        idx <- sample.int(n, n, replace = TRUE)
        out[i] <- pearson_distance(sub$ref_log2[idx], sub$target_log2[idx])
      }
    }
    bt <- data.frame(target_group = tg, feature_class = fc,
                     iter = seq_len(B), pearson_distance = out,
                     n_genes = n, observed_pearson_distance = observed)
    boots[[k]] <- bt
    sums[[k]] <- data.frame(target_group = tg, feature_class = fc, n_genes = n,
                            observed_pearson_distance = observed,
                            boot_mean = mean(out, na.rm = TRUE),
                            boot_median = median(out, na.rm = TRUE),
                            boot_ci95_low = quantile(out, 0.025, na.rm = TRUE),
                            boot_ci95_high = quantile(out, 0.975, na.rm = TRUE))
    k <- k + 1
  }
}
boot <- bind_rows(boots) %>%
  mutate(target_group = factor(target_group, levels = targets),
         feature_class = factor(feature_class, levels = feature_levels))
sumtab <- bind_rows(sums) %>%
  mutate(target_group = factor(target_group, levels = targets),
         feature_class = factor(feature_class, levels = feature_levels))

write.table(boot,
            file.path(tabdir, 'read_mapping_primary_gene_structure_pearson_distance.bootstrap.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)
write.table(sumtab,
            file.path(tabdir, 'read_mapping_primary_gene_structure_pearson_distance.summary.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

pair_defs <- data.frame(group1 = 'exon',
                        group2 = c('intron', '2kb downstream', '2kb upstream'),
                        stringsAsFactors = FALSE)
pvals <- list()
k <- 1
for (tg in targets) {
  ymax <- max(boot$pearson_distance[boot$target_group == tg], na.rm = TRUE)
  yrng <- diff(range(boot$pearson_distance[boot$target_group == tg], na.rm = TRUE))
  if (!is.finite(yrng) || yrng == 0) yrng <- 0.02
  for (j in seq_len(nrow(pair_defs))) {
    g1 <- pair_defs$group1[j]
    g2 <- pair_defs$group2[j]
    x <- boot$pearson_distance[boot$target_group == tg & boot$feature_class == g1]
    y <- boot$pearson_distance[boot$target_group == tg & boot$feature_class == g2]
    pv <- suppressWarnings(wilcox.test(x, y)$p.value)
    pvals[[k]] <- data.frame(target_group = tg, group1 = g1, group2 = g2,
                             p_value = pv, label = format_p(pv),
                             x1 = match(g1, feature_levels),
                             x2 = match(g2, feature_levels),
                             xmid = mean(c(match(g1, feature_levels), match(g2, feature_levels))),
                             y = ymax + yrng * (0.10 + 0.10 * (j - 1)),
                             y_tick = ymax + yrng * (0.08 + 0.10 * (j - 1)),
                             y_text = ymax + yrng * (0.115 + 0.10 * (j - 1)))
    k <- k + 1
  }
}
pvals <- bind_rows(pvals) %>%
  mutate(target_group = factor(target_group, levels = targets))
write.table(pvals,
            file.path(tabdir, 'read_mapping_primary_gene_structure_pearson_distance.wilcox_exon_vs_others.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

pal <- c('exon' = '#f26b5b',
         'intron' = '#4db6ac',
         '2kb downstream' = '#7e9ed8',
         '2kb upstream' = '#d9a441')

p <- ggplot(boot, aes(x = feature_class, y = pearson_distance, color = feature_class, fill = feature_class)) +
  geom_violin(trim = FALSE, linewidth = 0.34, alpha = 0.14, scale = 'width', width = 0.88, na.rm = TRUE) +
  geom_boxplot(width = 0.20, outlier.size = 0.12, fill = 'white', linewidth = 0.42, na.rm = TRUE) +
  geom_point(data = sumtab, aes(x = feature_class, y = boot_mean),
             inherit.aes = FALSE, shape = 95, size = 5.2, color = '#d62728') +
  geom_segment(data = pvals, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, color = 'black', linewidth = 0.28) +
  geom_segment(data = pvals, aes(x = x1, xend = x1, y = y, yend = y_tick),
               inherit.aes = FALSE, color = 'black', linewidth = 0.28) +
  geom_segment(data = pvals, aes(x = x2, xend = x2, y = y, yend = y_tick),
               inherit.aes = FALSE, color = 'black', linewidth = 0.28) +
  geom_text(data = pvals, aes(x = xmid, y = y_text, label = label),
            inherit.aes = FALSE, size = 2.1, color = 'black') +
  facet_wrap(~ target_group, ncol = 4, scales = 'free_y') +
  scale_color_manual(values = pal, guide = 'none') +
  scale_fill_manual(values = pal, guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.25))) +
  coord_cartesian(clip = 'off') +
  theme_classic(base_size = 9.5) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        panel.spacing.x = unit(1.0, 'lines'),
        panel.spacing.y = unit(1.2, 'lines'),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        plot.title = element_text(face = 'bold'),
        plot.margin = margin(5.5, 14, 5.5, 5.5)) +
  labs(x = NULL, y = 'Pearson distance',
       title = 'Expression divergence by SV position relative to gene structure',
       subtitle = 'Each target species is compared with genome_Msa; red bars indicate bootstrap means')

ggsave(file.path(figdir, 'read_mapping_primary_gene_structure_pearson_distance.violin_wilcox.pdf'),
       p, width = 10.5, height = 6.4, useDingbats = FALSE)
ggsave(file.path(figdir, 'read_mapping_primary_gene_structure_pearson_distance.violin_wilcox.png'),
       p, width = 10.5, height = 6.4, dpi = 450)
