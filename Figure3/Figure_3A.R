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


##### read in data #####
df_miniTRiX <- fread("Supplementary_tables/supplementary_table7.tsv")[-c(1),]
df_TRiX <- fread("Supplementary_tables/supplementary_table6.tsv")[-c(1),]

colnames(df_miniTRiX) <- as.character(df_miniTRiX[1,])
colnames(df_TRiX) <- as.character(df_TRiX[1,])

df_miniTRiX <- df_miniTRiX[-1,]
df_TRiX <- df_TRiX[-1,]

melt_df_miniTRiX <- melt(data.table(df_miniTRiX), id.vars = c("sample"), measure.vars = c("mean_skewing","TRiXi1_skewing", "TRiXi3_skewing"))
melt_df_miniTRiX$value <- as.numeric(melt_df_miniTRiX$value)

melt_df_TRiX <- melt(data.table(df_TRiX), id.vars = c("sample"), measure.vars = c("mean_skewing", "TRiXi1_skewing", "TRiXi2_skewing", "TRiXi3_skewing", "TRiXi4_skewing", "TRiXi5_skewing"))
melt_df_TRiX$value <- as.numeric(melt_df_TRiX$value)

# rbind
TRiXi_data <- rbind(melt_df_miniTRiX,
                    melt_df_TRiX)

# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample) %>% rstatix::get_summary_stats(value, type = "common")

# calculate descriptive statistics of descriptive statistics, to plot lines mean, and 1sd and 2sd.
stats_stats_TRiXi_data <- stats_TRiXi_data %>% rstatix::get_summary_stats(mean, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "non-skewed"
stats_TRiXi_data[stats_TRiXi_data$mean >= 75,]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean >= 90,]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "non-skewed"))

############ GTEx ############
suppl_table_10 <- fread("Supplementary_tables/supplementary_table10.tsv")

# calculate metrics for plotting.
suppl_table_10_metrics <- suppl_table_10 %>% dplyr::group_by(participant, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")
suppl_table_10_metrics$tissue_id <- factor(suppl_table_10_metrics$tissue_id)

# get stats
stats_suppl_table_10_metrics_order <- suppl_table_10_metrics %>% dplyr::group_by(participant) %>% rstatix::get_summary_stats(median, type = "common")

# Export source data
source_data_fig_3a <- stats_suppl_table_10_metrics_order
write.table(source_data_fig_3a, "source_data/Fig_3b_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_3a <- fread("source_data/Fig_3A_source_data.tsv")

# calculate descriptive statistics of descriptive statistics.
stats_stats_suppl_table_10_metrics_order <- source_data_fig_3a %>% rstatix::get_summary_stats(median, type = "common")

#### Figure 3a, showing distribution of skewing in GTEx ####
fig_3a <- 
  ggplot(source_data_fig_3a, 
         aes(x=1, y=median)) +
  geom_quasirandom(size = 1, alpha = .5) +
  geom_hline(yintercept = stats_stats_suppl_table_10_metrics_order$median, col = "black") +
  geom_hline(yintercept = 0.25, col = "pink")+
  geom_hline(yintercept = 0.40, col = "red")+
  annotate(geom = "text", x=0.7, y = c(stats_stats_suppl_table_10_metrics_order$median, 0.25, 0.40), label = c("median", ">=75", ">=90"))+
  coord_cartesian(ylim=c(0,0.5)) +
  theme_AL_box() + 
  theme(axis.text.x = element_blank()) + 
  labs(x="")+#, y="skewing (median nonPAR allele-specific expression)")+
  theme(legend.position = "top")+
  scale_y_continuous(breaks = c(0,0.1, 0.2,0.3,0.4, 0.5))+
  geom_text_repel(data=source_data_fig_3a[source_data_fig_3a$median > (stats_stats_suppl_table_10_metrics_order$median+2*stats_stats_suppl_table_10_metrics_order$sd) & source_data_fig_3a$participant %in% c("GTEX-UPIC", "GTEX-13PLJ", "GTEX-ZZPU"),], 
                  aes(label=participant), min.segment.length = 0.001, size = 3, box.padding = 0.1, nudge_x = 0.25)#, nudge_y = -0.05)

# get metrics
stats_source_data_fig_3a <- source_data_fig_3a

stats_source_data_fig_3a$'XCI pattern' <- "non-skewed"
stats_source_data_fig_3a[stats_source_data_fig_3a$median >= 0.25,]$'XCI pattern' <- "skewed"
stats_source_data_fig_3a[stats_source_data_fig_3a$median >= 0.4,]$'XCI pattern' <- "extremely skewed"
stats_source_data_fig_3a$'XCI pattern' <- factor(stats_source_data_fig_3a$'XCI pattern', levels = c("extremely skewed", "skewed", "non-skewed"))

stats_source_data_fig_3a %>% dplyr::group_by(`XCI pattern`) %>% dplyr::count()


ggsave2(filename = paste0("01_plots/Figure_3A.pdf"), height = 116, width = 120, units = "mm",
        plot_grid(fig_3a
)
)

