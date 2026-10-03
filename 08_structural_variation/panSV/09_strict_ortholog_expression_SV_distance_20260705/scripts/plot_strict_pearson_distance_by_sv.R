suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
set.seed(20260705)
args <- commandArgs(trailingOnly=FALSE)
file_arg <- sub('--file=', '', args[grep('--file=', args)][1])
work <- normalizePath(file.path(dirname(file_arg), '..'), mustWork=TRUE)
tabdir <- file.path(work, 'tables')
figdir <- file.path(work, 'figures')
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

B <- 10000L
infile <- file.path(tabdir, 'strict_ortholog_expression_divergence.expressed_genes.tsv')
df <- read.delim(infile, check.names=FALSE, stringsAsFactors=FALSE)
df <- df %>%
  mutate(ref_log2 = log2(mean_TPM_ref + 1),
         target_log2 = log2(mean_TPM_target + 1),
         target_group = factor(target_group, levels=c('genome_410','genome_461','genome_472','genome_474','genome_Mar','genome_Mru','genome_R108','genome_ZM4')),
         sv_status = factor(sv_status, levels=c('Without SV','With SV')))

pearson_distance <- function(x, y) {
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  1 - suppressWarnings(cor(x, y, method='pearson'))
}

boot_one <- function(sub, B=10000L) {
  n <- nrow(sub)
  x <- sub$ref_log2
  y <- sub$target_log2
  observed <- pearson_distance(x, y)
  out <- numeric(B)
  for (i in seq_len(B)) {
    idx <- sample.int(n, n, replace=TRUE)
    out[i] <- pearson_distance(x[idx], y[idx])
  }
  data.frame(iter=seq_len(B), pearson_distance=out, observed=observed, n_genes=n)
}

pieces <- list()
summary_rows <- list()
k <- 1L
for (cmp in levels(droplevels(df$target_group))) {
  for (st in levels(df$sv_status)) {
    sub <- df %>% filter(target_group == cmp, sv_status == st)
    if (nrow(sub) < 20) next
    bt <- boot_one(sub, B)
    bt$target_group <- cmp
    bt$comparison <- unique(sub$comparison)[1]
    bt$sv_status <- st
    pieces[[k]] <- bt
    summary_rows[[k]] <- data.frame(
      target_group=cmp,
      comparison=unique(sub$comparison)[1],
      sv_status=st,
      n_genes=nrow(sub),
      observed_pearson_distance=bt$observed[1],
      boot_mean=mean(bt$pearson_distance, na.rm=TRUE),
      boot_median=median(bt$pearson_distance, na.rm=TRUE),
      boot_ci95_low=quantile(bt$pearson_distance, 0.025, na.rm=TRUE),
      boot_ci95_high=quantile(bt$pearson_distance, 0.975, na.rm=TRUE)
    )
    k <- k + 1L
  }
}
boot <- bind_rows(pieces)
sumtab <- bind_rows(summary_rows)
write.table(boot, file.path(tabdir, 'strict_pearson_distance_by_sv.bootstrap.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(sumtab, file.path(tabdir, 'strict_pearson_distance_by_sv.summary.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

pvals <- boot %>%
  group_by(target_group, comparison) %>%
  summarise(
    p_value = suppressWarnings(wilcox.test(pearson_distance[sv_status=='Without SV'], pearson_distance[sv_status=='With SV'])$p.value),
    y = max(pearson_distance, na.rm=TRUE) + 0.03,
    .groups='drop'
  ) %>%
  mutate(label = ifelse(p_value < 2.2e-16, 'p < 2.2e-16', paste0('p = ', formatC(p_value, format='e', digits=1))))
write.table(pvals, file.path(tabdir, 'strict_pearson_distance_by_sv.wilcox.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

boot$target_group <- factor(boot$target_group, levels=levels(df$target_group))
boot$sv_status <- factor(boot$sv_status, levels=c('Without SV','With SV'))
pvals$target_group <- factor(pvals$target_group, levels=levels(df$target_group))
pal <- c('Without SV'='#0b6fad', 'With SV'='#d8aa00')

p <- ggplot(boot, aes(x=target_group, y=pearson_distance, color=sv_status)) +
  geom_boxplot(width=0.58, outlier.size=0.18, outlier.alpha=0.20, fill='white',
               linewidth=0.55, position=position_dodge(width=0.72)) +
  stat_summary(fun=mean, geom='crossbar', width=0.42, fatten=0, color='#d62728',
               linewidth=0.35, position=position_dodge(width=0.72)) +
  geom_text(data=pvals, aes(x=target_group, y=y, label=label), inherit.aes=FALSE,
            angle=90, size=2.7, vjust=0.5, hjust=0) +
  scale_color_manual(values=pal) +
  coord_cartesian(clip='off') +
  labs(x=NULL, y='Pearson distance', color=NULL,
       title='Pearson distance between Msa and target species by nearby SV status') +
  theme_classic(base_size=8) +
  theme(legend.position='top',
        legend.key.width=unit(0.35, 'cm'),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        plot.title=element_text(hjust=0.5, size=9),
        plot.margin=margin(6, 8, 6, 8))

ggsave(file.path(figdir, 'strict_pearson_distance_by_sv.boxplot.pdf'), p, width=7.4, height=4.5, useDingbats=FALSE)
ggsave(file.path(figdir, 'strict_pearson_distance_by_sv.boxplot.png'), p, width=7.4, height=4.5, dpi=450)
cat('Pearson-distance figure written to', figdir, '\n')
