suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
args <- commandArgs(trailingOnly=FALSE)
file_arg <- sub('--file=', '', args[grep('--file=', args)][1])
work <- normalizePath(file.path(dirname(file_arg), '..'), mustWork=TRUE)
tabdir <- file.path(work, 'tables')
figdir <- file.path(work, 'figures')
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

df <- read.delim(file.path(tabdir, 'strict_ortholog_expression_divergence.expressed_genes.tsv'), check.names=FALSE, stringsAsFactors=FALSE)
lev <- c('genome_410','genome_461','genome_472','genome_474','genome_Mar','genome_Mru','genome_R108','genome_ZM4')
df <- df %>%
  mutate(target_group=factor(target_group, levels=lev),
         sv_status=factor(sv_status, levels=c('Without SV','With SV')),
         ref_expr_log2=log2(mean_TPM_ref + 1),
         target_expr_log2=log2(mean_TPM_target + 1))

long <- df %>%
  select(gene_id, comparison, target_group, sv_status, ref_expr_log2, target_expr_log2) %>%
  pivot_longer(cols=c(ref_expr_log2, target_expr_log2), names_to='expression_source', values_to='expression_log2_TPM1') %>%
  mutate(expression_source=recode(expression_source,
                                  ref_expr_log2='Msa reference expression',
                                  target_expr_log2='Target ortholog expression'))
write.table(long, file.path(tabdir, 'strict_expression_level_by_sv.long.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

sumtab <- long %>%
  group_by(expression_source, target_group, comparison, sv_status) %>%
  summarise(n_genes=n(), median_expression=median(expression_log2_TPM1, na.rm=TRUE),
            mean_expression=mean(expression_log2_TPM1, na.rm=TRUE),
            q25=quantile(expression_log2_TPM1, 0.25, na.rm=TRUE),
            q75=quantile(expression_log2_TPM1, 0.75, na.rm=TRUE), .groups='drop')
write.table(sumtab, file.path(tabdir, 'strict_expression_level_by_sv.summary.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

pvals <- long %>%
  group_by(expression_source, target_group, comparison) %>%
  summarise(p_value=suppressWarnings(wilcox.test(expression_log2_TPM1[sv_status=='Without SV'], expression_log2_TPM1[sv_status=='With SV'])$p.value),
            y=max(expression_log2_TPM1, na.rm=TRUE) + 0.35, .groups='drop') %>%
  mutate(label=ifelse(p_value < 2.2e-16, 'p < 2.2e-16', paste0('p = ', formatC(p_value, format='e', digits=1))))
write.table(pvals, file.path(tabdir, 'strict_expression_level_by_sv.wilcox.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

target <- long %>% filter(expression_source == 'Target ortholog expression')
pvals_target <- pvals %>% filter(expression_source == 'Target ortholog expression')
pal <- c('Without SV'='#0b6fad', 'With SV'='#d8aa00')

theme_pub <- theme_classic(base_size=8) +
  theme(legend.position='top', legend.title=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        plot.title=element_text(hjust=0.5, size=9),
        strip.background=element_blank(), strip.text=element_text(face='bold'),
        plot.margin=margin(5, 8, 5, 8))

p1 <- ggplot(target, aes(x=target_group, y=expression_log2_TPM1, color=sv_status)) +
  geom_boxplot(width=0.58, outlier.size=0.15, outlier.alpha=0.18, fill='white',
               linewidth=0.5, position=position_dodge(width=0.72)) +
  geom_text(data=pvals_target, aes(x=target_group, y=y, label=label), inherit.aes=FALSE,
            angle=90, size=2.5, vjust=0.5, hjust=0) +
  scale_color_manual(values=pal) + coord_cartesian(clip='off') +
  labs(x=NULL, y=expression('Expression level (log'[2]*'(TPM + 1))'),
       title='Target-species ortholog expression by nearby SV status') +
  theme_pub

ggsave(file.path(figdir, 'strict_target_expression_level_by_sv.boxplot.pdf'), p1, width=7.4, height=4.5, useDingbats=FALSE)
ggsave(file.path(figdir, 'strict_target_expression_level_by_sv.boxplot.png'), p1, width=7.4, height=4.5, dpi=450)

p2 <- ggplot(long, aes(x=target_group, y=expression_log2_TPM1, color=sv_status)) +
  geom_boxplot(width=0.58, outlier.size=0.12, outlier.alpha=0.15, fill='white',
               linewidth=0.45, position=position_dodge(width=0.72)) +
  facet_wrap(~expression_source, ncol=1) +
  scale_color_manual(values=pal) +
  labs(x=NULL, y=expression('Expression level (log'[2]*'(TPM + 1))'),
       title='Expression level in Msa and target species by nearby SV status') +
  theme_pub

ggsave(file.path(figdir, 'strict_ref_target_expression_level_by_sv.boxplot.pdf'), p2, width=7.4, height=6.8, useDingbats=FALSE)
ggsave(file.path(figdir, 'strict_ref_target_expression_level_by_sv.boxplot.png'), p2, width=7.4, height=6.8, dpi=450)

cat('expression-level figures written to', figdir, '\n')
