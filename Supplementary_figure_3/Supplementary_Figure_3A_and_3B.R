library(data.table)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggbeeswarm)
library(cowplot)
library(ggalluvial)
library(ggrepel)
library(GenomicRanges)
library(rtracklayer)

source("plot_parameters_TRiX.R")

############ GTEx ############
suppl_table_10 <- fread("Supplementary_tables/supplementary_table10.tsv")

# calculate metrics for plotting.
plot_metrics <- suppl_table_10 %>% dplyr::group_by(participant, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")
plot_metrics$tissue_id <- factor(plot_metrics$tissue_id)

#### Figure 3E, showing distribution of skewing across tissues ####
plot_metrics$skewed <- ifelse(plot_metrics$median >= 0.25, yes = "skewed", no = "not_skewed")

# count and make wide. Add total count column.
skewed_tissue_samples_count <- plot_metrics %>% dplyr::group_by(tissue_id, skewed) %>% dplyr::count(name = "skewed_count")

# Export source data
suppl_fig_3A_source_data <- skewed_tissue_samples_count
write.table(suppl_fig_3A_source_data, "source_data/Suppl_Fig_3A_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
suppl_fig_3A_source_data <- fread("source_data/Suppl_Fig_3A_source_data.tsv")

# add percentage of total skewed
skewed_tissue_samples_count_dcast <- dcast(setDT(suppl_fig_3A_source_data), tissue_id ~ skewed, value.var = c("skewed_count"))
skewed_tissue_samples_count_dcast$total <- skewed_tissue_samples_count_dcast$not_skewed + skewed_tissue_samples_count_dcast$skewed
skewed_tissue_samples_count_dcast$perc_skewed <- skewed_tissue_samples_count_dcast$skewed / skewed_tissue_samples_count_dcast$total

# make plot order
skewed_tissue_samples_count_dcast_order <- skewed_tissue_samples_count_dcast[order(skewed_tissue_samples_count_dcast$perc_skewed),]$tissue_id



distribution_skewing_across_tissues_top <- 
  ggplot(suppl_fig_3A_source_data,
         aes(x=factor(tissue_id, levels = skewed_tissue_samples_count_dcast_order), y = skewed_count, fill = skewed))+
  geom_bar(stat="identity", position = "fill") +
  theme_AL_box_rotX() +
  labs(y="fraction skewed (%)", x="")+
  geom_text(aes(label=skewed_count), size = 2.5, position = "fill")+
  theme(axis.text.y = element_blank(), legend.position = "top") #, axis.text.x = element_blank())


# calculate metrics for plotting.
plot_metrics <- plot_metrics
plot_metrics$tissue_id <- factor(plot_metrics$tissue_id)

plot_metrics$tissue_id <- factor(plot_metrics$tissue_id, levels = skewed_tissue_samples_count_dcast_order)
plot_metrics <- plot_metrics[order(plot_metrics$tissue_id),]

# Export source data
suppl_fig_3B_source_data <- plot_metrics
write.table(suppl_fig_3B_source_data, "source_data/Suppl_Fig_3B_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
suppl_fig_3B_source_data <- fread("source_data/Suppl_Fig_3B_source_data.tsv")

bottom_plot <- 
  ggplot(suppl_fig_3B_source_data, 
                      aes(y = median, x = factor(tissue_id, levels = skewed_tissue_samples_count_dcast_order))) +
  geom_quasirandom(aes(col=tissue_id), alpha = .5, fill = "black", width = 0.25) +
  theme_AL_box_rotX() +
  theme(legend.position = "none")+
  geom_hline(yintercept = c(0.25), lty=2)+
  theme(axis.text.y = element_blank())

ggsave2(filename = "01_plots/Supplementary_Figure_3A_3B.pdf", width = 8, height = 7,
        plot_grid(ncol=1, rel_heights = c(0.5,0.5),
                  
                  distribution_skewing_across_tissues_top,
                  bottom_plot)
)

