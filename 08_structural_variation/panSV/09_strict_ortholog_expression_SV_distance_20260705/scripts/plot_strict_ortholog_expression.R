suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
# Fallback for older R if pipe placeholder is unavailable
args <- commandArgs(trailingOnly=FALSE)
file_arg <- sub('--file=', '', args[grep('--file=', args)][1])
work <- normalizePath(file.path(dirname(file_arg), '..'), mustWork=TRUE)
tabdir <- file.path(work, 'tables')
figdir <- file.path(work, 'figures')
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)
strict <- read.delim(file.path(tabdir, 'strict_ortholog_expression_divergence.expressed_genes.tsv'), check.names=FALSE)
strict$target_group <- factor(strict$target_group, levels=c('genome_410','genome_461','genome_472','genome_474','genome_Mar','genome_Mru','genome_R108','genome_ZM4'))
strict$sv_status <- factor(strict$sv_status, levels=c('Without SV','With SV'))
pal <- c('Without SV'='#2b6aa0','With SV'='#d4a000')
p <- ggplot(strict, aes(x=target_group, y=expression_divergence, color=sv_status)) +
  geom_boxplot(width=0.62, outlier.size=0.25, outlier.alpha=0.18, fill='white', position=position_dodge(width=0.72), linewidth=0.55) +
  scale_color_manual(values=pal) +
  labs(x=NULL, y=expression('Expression divergence  |log'[2]*'(TPM + 1) difference|'), color=NULL,
       title='Strict one-to-one ortholog expression divergence') +
  theme_classic(base_size=8) +
  theme(legend.position='top', axis.text.x=element_text(angle=45, hjust=1, vjust=1), plot.title=element_text(hjust=0.5, size=9))
ggsave(file.path(figdir, 'strict_ortholog_expression_divergence.boxplot.pdf'), p, width=7.2, height=4.2, useDingbats=FALSE)
ggsave(file.path(figdir, 'strict_ortholog_expression_divergence.boxplot.png'), p, width=7.2, height=4.2, dpi=450)

sumtab <- read.delim(file.path(tabdir, 'strict_vs_reference_projected.expression_divergence.summary.tsv'), check.names=FALSE)
sumtab$method <- factor(sumtab$method, levels=c('08_reference_projected','09_strict_one_to_one_ortholog'), labels=c('08 reference-projected','09 strict ortholog'))
sumtab$target_group <- factor(sumtab$target_group, levels=levels(strict$target_group))
sumtab$sv_status <- factor(sumtab$sv_status, levels=c('Without SV','With SV'))
p2 <- ggplot(sumtab, aes(x=target_group, y=median_expression_divergence, group=method, color=method, shape=method)) +
  geom_line(linewidth=0.45) + geom_point(size=1.7) +
  facet_wrap(~sv_status, ncol=1, scales='free_y') +
  scale_color_manual(values=c('#777777','#c0392b')) +
  labs(x=NULL, y='Median expression divergence', color=NULL, shape=NULL,
       title='Reference-projected versus strict ortholog result') +
  theme_classic(base_size=8) +
  theme(legend.position='top', axis.text.x=element_text(angle=45, hjust=1, vjust=1), plot.title=element_text(hjust=0.5, size=9))
ggsave(file.path(figdir, 'strict_vs_reference_projected.median_comparison.pdf'), p2, width=7.2, height=5.2, useDingbats=FALSE)
ggsave(file.path(figdir, 'strict_vs_reference_projected.median_comparison.png'), p2, width=7.2, height=5.2, dpi=450)

cat('figures written to', figdir, '\n')
