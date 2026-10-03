suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
work <- 'path/to/project/N_3.call_SV/10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705'
tabdir <- file.path(work, 'tables')
figdir <- file.path(work, 'figures')
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)
df <- read.delim(file.path(tabdir, 'read_mapping_primary_exon_CDS_UTR_promoter2kb_expression_status.tsv'), check.names=FALSE, stringsAsFactors=FALSE)
lev <- c('genome_410','genome_461','genome_472','genome_474','genome_Mar','genome_Mru','genome_R108','genome_ZM4')
df <- df %>% mutate(target_group=factor(target_group, levels=lev),
                    sv_status=factor(sv_status, levels=c('Without SV','With SV')),
                    ref_log2=log2(mean_TPM_ref+1),
                    target_log2=log2(mean_TPM_target+1))
pearson_distance <- function(x, y) {
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  1 - suppressWarnings(cor(x, y, method='pearson'))
}
set.seed(1)
B <- 10000
boot_one <- function(sub) {
  x <- sub$ref_log2; y <- sub$target_log2; n <- length(x)
  observed <- pearson_distance(x, y)
  out <- rep(NA_real_, B)
  if (n >= 2) {
    for (i in seq_len(B)) {
      idx <- sample.int(n, n, replace=TRUE)
      out[i] <- pearson_distance(x[idx], y[idx])
    }
  }
  data.frame(iter=seq_len(B), pearson_distance=out, observed=observed, n_genes=n)
}
boots <- list(); sums <- list(); k <- 1
for (cmp in lev) {
  for (st in c('Without SV','With SV')) {
    sub <- df %>% filter(target_group == cmp, sv_status == st)
    bt <- boot_one(sub)
    bt$target_group <- cmp; bt$sv_status <- st
    boots[[k]] <- bt
    sums[[k]] <- data.frame(target_group=cmp, sv_status=st, n_genes=nrow(sub),
                            observed_pearson_distance=bt$observed[1],
                            boot_mean=mean(bt$pearson_distance, na.rm=TRUE),
                            boot_median=median(bt$pearson_distance, na.rm=TRUE),
                            boot_ci95_low=quantile(bt$pearson_distance, 0.025, na.rm=TRUE),
                            boot_ci95_high=quantile(bt$pearson_distance, 0.975, na.rm=TRUE))
    k <- k + 1
  }
}
boot <- bind_rows(boots)
sumtab <- bind_rows(sums)
write.table(boot, file.path(tabdir, 'read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson_distance.bootstrap.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(sumtab, file.path(tabdir, 'read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson_distance.summary.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
pvals <- boot %>% group_by(target_group) %>%
  summarise(p_value=suppressWarnings(wilcox.test(pearson_distance[sv_status=='Without SV'], pearson_distance[sv_status=='With SV'])$p.value),
            y=max(pearson_distance, na.rm=TRUE)+0.012, .groups='drop') %>%
  mutate(label=ifelse(is.na(p_value),'NA',ifelse(p_value<2.2e-16,'P < 2.2e-16',paste0('P = ', signif(p_value,2)))))
write.table(pvals, file.path(tabdir, 'read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson_distance.wilcox.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
boot$target_group <- factor(boot$target_group, levels=lev)
boot$sv_status <- factor(boot$sv_status, levels=c('Without SV','With SV'))
sumtab$target_group <- factor(sumtab$target_group, levels=lev)
sumtab$sv_status <- factor(sumtab$sv_status, levels=c("Without SV","With SV"))
pvals$target_group <- factor(pvals$target_group, levels=lev)
pvals <- pvals %>% mutate(x=as.numeric(target_group), xmin=x-0.18, xmax=x+0.18, y_tick=y-0.004, y_text=y+0.003)
pal <- c('Without SV'='#0b6fad', 'With SV'='#d8aa00')
p <- ggplot(boot, aes(x=target_group, y=pearson_distance, color=sv_status)) +
  geom_boxplot(position=position_dodge(width=0.72), width=0.55, outlier.size=0.15, fill='white', linewidth=0.65) +
  geom_point(data=sumtab, aes(x=target_group, y=boot_mean, group=sv_status), position=position_dodge(width=0.72), shape=95, size=6, color='#d62728', inherit.aes=FALSE) +
  geom_segment(data=pvals, aes(x=xmin, xend=xmax, y=y, yend=y), inherit.aes=FALSE, color="black", linewidth=0.32) +
  geom_segment(data=pvals, aes(x=xmin, xend=xmin, y=y, yend=y_tick), inherit.aes=FALSE, color="black", linewidth=0.32) +
  geom_segment(data=pvals, aes(x=xmax, xend=xmax, y=y, yend=y_tick), inherit.aes=FALSE, color="black", linewidth=0.32) +
  geom_text(data=pvals, aes(x=x, y=y_text, label=label), inherit.aes=FALSE, color="black", size=2.7, angle=90, hjust=0, vjust=0.5) +
  scale_y_continuous(expand=expansion(mult=c(0.04, 0.22))) +
  scale_color_manual(values=pal) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=40, hjust=1), legend.position='top') +
  labs(x=NULL, y='Pearson distance', color=NULL, title='Read-mapping primary SVs in exon/CDS/UTR/promoter 2 kb')
ggsave(file.path(figdir, 'read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson_distance.boxplot.pdf'), p, width=7.4, height=4.5, useDingbats=FALSE)
ggsave(file.path(figdir, 'read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson_distance.boxplot.png'), p, width=7.4, height=4.5, dpi=450)
