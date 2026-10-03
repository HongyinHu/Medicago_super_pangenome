suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})
out_dir <- "path/to/project/39.TE_type_soloLTR/output/publication_ltr_insertion_time"
fig_dir <- file.path(out_dir, "figures")
records <- read_tsv(file.path(out_dir, "tables", "ltr_insertion_time_passlist_records.tsv"), show_col_types = FALSE) %>%
  mutate(insertion_time_mya=as.numeric(insertion_time_mya), species_order=as.integer(species_order), label=if_else(is.na(label)|label=="", species, label),
         ltr_group=recode(ltr_class, Gypsy="LTRs-Gypsy", Copia="LTRs-Copia", Unknown="LTRs-unknown")) %>%
  filter(!is.na(insertion_time_mya), insertion_time_mya > 0, insertion_time_mya <= 8,
         !species %in% c("genome_474_T2T", "genome_A17", "genome_Msa_T2T"))

bin_width <- 0.10
panel_levels <- c("LTRs-Gypsy", "LTRs-Copia", "LTRs-unknown", "LTRs")
panel_labels <- c("a  LTRs-Gypsy", "b  LTRs-Copia", "c  LTRs-unknown", "d  LTRs")
labels <- records %>% distinct(species,label,species_order) %>% arrange(species_order,label)
time_grid <- seq(bin_width/2, 8 - bin_width/2, by=bin_width)
base_grid <- expand_grid(species=labels$species, time_bin=time_grid) %>% left_join(labels, by="species")
by_class <- records %>% mutate(time_bin=floor(insertion_time_mya/bin_width)*bin_width + bin_width/2) %>% count(species,label,species_order,ltr_group,time_bin,name="count") %>% rename(panel=ltr_group)
all_ltr <- records %>% mutate(time_bin=floor(insertion_time_mya/bin_width)*bin_width + bin_width/2) %>% count(species,label,species_order,time_bin,name="count") %>% mutate(panel="LTRs")
binned0 <- bind_rows(by_class, all_ltr) %>% mutate(panel=factor(panel, levels=panel_levels, labels=panel_labels))
# Fill missing bins with zero per species/panel.
binned <- expand_grid(labels, panel=factor(panel_labels, levels=panel_labels), time_bin=time_grid) %>%
  left_join(binned0 %>% select(species,panel,time_bin,count), by=c("species","panel","time_bin")) %>%
  mutate(count=replace_na(count,0)) %>%
  group_by(species,label,species_order,panel) %>% arrange(time_bin, .by_group=TRUE) %>%
  mutate(count_smooth=as.numeric(stats::filter(count, rep(1/3,3), sides=2)), count_smooth=if_else(is.na(count_smooth), as.numeric(count), count_smooth)) %>%
  ungroup()

pal <- c("genome_395"="#2C7FB8","genome_436"="#F28E2B","genome_457"="#59A14F","genome_468"="#D62728","genome_474"="#8F63B0","genome_482"="#8C564B","genome_M22"="#E377C2","Marc"="#7F7F7F","Mrut"="#B5BD00","Msat_Cae"="#17BECF","Msat_ZM4"="#111111","genome_410"="#74A9CF","genome_454"="#FDB27A","genome_461"="#74C476","genome_472"="#FB6A4A","genome_M46"="#9E9AC8","Mpol"="#C49C94","Mtru_R108"="#F2A6C8")
levels_label <- labels$label
palm <- pal[levels_label]; miss <- is.na(palm); if(any(miss)) palm[miss] <- grDevices::hcl.colors(sum(miss), "Dark 3"); names(palm) <- levels_label
binned <- binned %>% mutate(label=factor(label, levels=levels_label))
size_map <- rep(0.48, length(levels_label)); names(size_map)<-levels_label; size_map["Msat_ZM4"] <- 0.95
alpha_map <- rep(0.78, length(levels_label)); names(alpha_map)<-levels_label; alpha_map["Msat_ZM4"] <- 1.0
p <- ggplot(binned, aes(time_bin, count_smooth, colour=label, group=label)) +
  geom_point(aes(y=count, alpha=label), size=0.28, stroke=0, show.legend=FALSE) +
  geom_line(aes(linewidth=label, alpha=label), lineend="round") +
  facet_wrap(~panel, ncol=2, scales="free_y") +
  scale_colour_manual(values=palm, name=NULL) + scale_linewidth_manual(values=size_map, guide="none") + scale_alpha_manual(values=alpha_map, guide="none") +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), expand=expansion(mult=c(0,0.01))) +
  labs(x="Insertion time (Mya)", y="Number of intact LTR-RTs per 0.1 Mya") +
  guides(colour=guide_legend(nrow=3, byrow=TRUE, override.aes=list(linewidth=1.0, alpha=1))) +
  theme_classic(base_size=7, base_family="Arial") +
  theme(axis.line=element_line(linewidth=0.32, colour="#262626"), axis.ticks=element_line(linewidth=0.28, colour="#262626"),
        axis.text=element_text(size=6.2, colour="#222222"), axis.title=element_text(size=7.0, colour="#111111"),
        strip.background=element_blank(), strip.text=element_text(size=7.2, colour="#111111", margin=margin(b=3)),
        panel.spacing.x=unit(10,"mm"), panel.spacing.y=unit(8,"mm"), legend.position="bottom",
        legend.text=element_text(size=5.4), legend.key.width=unit(6,"mm"), legend.key.height=unit(2.8,"mm"),
        panel.grid.major.y=element_line(linewidth=0.14, colour="#ECECEC"), panel.grid.minor=element_blank(), plot.margin=margin(4,5,3,5))
prefix <- file.path(fig_dir, "ltr_insertion_time.passlist.main18_with_Msat_ZM4.no_zero.v2_bin0.1")
w <- 178/25.4; h <- 145/25.4
cairo_pdf(paste0(prefix,".pdf"), width=w, height=h, family="Arial"); print(p); dev.off()
ragg::agg_png(paste0(prefix,".png"), width=w, height=h, units="in", res=700, background="white"); print(p); dev.off()
ragg::agg_tiff(paste0(prefix,".tiff"), width=w, height=h, units="in", res=700, background="white", compression="lzw"); print(p); dev.off()
svg(paste0(prefix,".svg"), width=w, height=h, family="Arial"); print(p); dev.off()
write_tsv(binned %>% arrange(panel,species_order,label,time_bin), paste0(prefix,".binned_counts.tsv"))
