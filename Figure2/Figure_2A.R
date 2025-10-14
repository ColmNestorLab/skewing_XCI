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
df_miniTRiX <- fread("data/minitrix_results.tsv")
df_TRiX <- fread("data/trix_results.tsv")
df_miniTRiX_age <- fread("data/minitrix_AGE_results.tsv")
df_HUMARA <- fread("data/HUMARA_results.tsv")
df_miniTRiX_HUMARA_matched <- fread("data/minitrix_HUMARA_matched_results.tsv")

df_miniTRiX$method <- "mini"
df_TRiX$method <- "full"

melt_df_miniTRiX <- melt(data.table(df_miniTRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing"))
melt_df_miniTRiX$value <- as.numeric(melt_df_miniTRiX$value)

melt_df_TRiX <- melt(data.table(df_TRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing", "CNKSR2_skewing", "TCAC1_skewing", "ZIC3_skewing"))
melt_df_TRiX$value <- as.numeric(melt_df_TRiX$value)

# rbind
TRiXi_data <- rbind(melt_df_miniTRiX,
                    melt_df_TRiX)



##### Figure 2A ####
# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample, method) %>% rstatix::get_summary_stats(value, type = "common")

# calculate descriptive statistics of descriptive statistics, to plot lines mean, and 1sd and 2sd.
stats_stats_TRiXi_data <- stats_TRiXi_data %>% rstatix::get_summary_stats(mean, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "mosaic"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+1*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "mosaic"))

# plot skewing
distribution_of_skewing_per_individual <- 
  ggplot(stats_TRiXi_data, 
         aes(x=reorder(sample, mean), y=mean, col = `XCI pattern`)) +
  coord_cartesian(ylim=c(45,100)) +
  geom_ribbon(aes(ymin=median-sd, ymax=median+sd, group = 1), col = "lightblue", fill = "lightblue", alpha = 0.1)+
  geom_point() +
  theme_AL_box() + 
  geom_hline(yintercept = stats_stats_TRiXi_data$mean, col = "black")+
  geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+stats_stats_TRiXi_data$sd), col = "pink")+
  geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd), col = "red")+
  geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+3*stats_stats_TRiXi_data$sd), col = "black")+
  theme(axis.text.x = element_blank()) + 
  labs(x="", y="skewing (mean ± sd %)")+
  theme(legend.position = "right") +
  expand_limits(x= c(-10, length(stats_TRiXi_data$sample) + 10))+
  scale_color_manual(values=colors)+
  annotate(geom = "text", x=15, y = c(stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd, stats_stats_TRiXi_data$mean+stats_stats_TRiXi_data$sd, stats_stats_TRiXi_data$mean+1.5, stats_stats_TRiXi_data$mean+3*stats_stats_TRiXi_data$sd), label = c("+2 SD", "+1 SD", "mean", "+3 SD"))+
  scale_y_continuous(breaks = c(50, 75, 100))+
  geom_text_repel(data=stats_TRiXi_data[stats_TRiXi_data$sample %in% c("320280"),], 
                  aes(label=sample), min.segment.length = 0.001, size = 3, box.padding = 0.1, col = "black")


ggsave2("01_plots/figure_2A.pdf", height = 0.5*116, width = 0.5*300, units = "mm",
        plot_grid(distribution_of_skewing_per_individual))

# get some summary stats
stats_stats_TRiXi_data <- stats_TRiXi_data %>% dplyr::group_by(`XCI pattern`) %>% dplyr::count()

# samples that are skewed or extremely skewed.
(53+107)/(869+(53+107))*100

# fraction extremely skewed
(53/1027)*100

# fraction skewed
(107/1027)*100
