library(data.table)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggbeeswarm)
library(cowplot)
library(ggalluvial)
library(ggrepel)

source("plot_parameters_TRiX.R")

##### read in data #####
df_miniTRiX <- fread("Supplementary_tables/supplementary_table5.tsv")[-c(1),]
df_TRiX <- fread("Supplementary_tables/supplementary_table4.tsv")[-c(1),]

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


##### Figure 2A ####
# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample) %>% rstatix::get_summary_stats(value, type = "common")

# calculate descriptive statistics of descriptive statistics, to plot lines mean, and 1sd and 2sd.
stats_stats_TRiXi_data <- stats_TRiXi_data %>% rstatix::get_summary_stats(mean, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "mosaic"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+1*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "mosaic"))

# Export source data
source_data_a <- stats_TRiXi_data
source_data_b <- stats_stats_TRiXi_data
write.table(source_data_a, "source_data/Fig_2A_source_data_a.tsv", quote = F, row.names = F, sep = "\t")
write.table(source_data_b, "source_data/Fig_2A_source_data_b.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_a <- fread("source_data/Fig_2A_source_data_a.tsv")
source_data_b <- fread("source_data/Fig_2A_source_data_b.tsv")

# plot skewing
distribution_of_skewing_per_individual <- 
  ggplot(source_data_a, 
         aes(x=reorder(sample, mean), y=mean, col = `XCI pattern`)) +
  coord_cartesian(ylim=c(45,100)) +
  geom_ribbon(aes(ymin=median-sd, ymax=median+sd, group = 1), col = "lightblue", fill = "lightblue", alpha = 0.1)+
  geom_point() +
  theme_AL_box() + 
  geom_hline(yintercept = source_data_b$mean, col = "black")+
  geom_hline(yintercept = c(source_data_b$mean+source_data_b$sd), col = "pink")+
  geom_hline(yintercept = c(source_data_b$mean+2*source_data_b$sd), col = "red")+
  geom_hline(yintercept = c(source_data_b$mean+3*source_data_b$sd), col = "black")+
  theme(axis.text.x = element_blank()) + 
  labs(x="", y="skewing (mean ± sd %)")+
  theme(legend.position = "right") +
  expand_limits(x= c(-10, length(source_data_a$sample) + 10))+
  scale_color_manual(values=colors)+
  annotate(geom = "text", x=15, y = c(source_data_b$mean+2*source_data_b$sd, source_data_b$mean+source_data_b$sd, source_data_b$mean+1.5, source_data_b$mean+3*source_data_b$sd), label = c("+2 SD", "+1 SD", "mean", "+3 SD"))+
  scale_y_continuous(breaks = c(50, 75, 100))+
  geom_text_repel(data=source_data_a[source_data_a$sample %in% c("320280"),], 
                  aes(label=sample), min.segment.length = 0.001, size = 3, box.padding = 0.1, col = "black")


ggsave2("01_plots/figure_2A.pdf", height = 0.5*116, width = 0.5*300, units = "mm",
        plot_grid(distribution_of_skewing_per_individual))

