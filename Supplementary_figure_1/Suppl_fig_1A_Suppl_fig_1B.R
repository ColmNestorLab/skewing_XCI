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


##### Figure 2F, showing skewing per repeat ####


# Get XCI pattern per sample (i.e. sample ID)
stats_TRiXi <- stats_TRiXi_data[,c("sample", "XCI pattern")]

# Make new df, merging the XCI pattern classification with it, Fix some variable names.
compare_distributions <- TRiXi_data[!is.na(TRiXi_data$value)]
compare_distributions <- merge(compare_distributions, stats_TRiXi, by = "sample")
compare_distributions$variable <- gsub("_skewing", "", compare_distributions$variable)
compare_distributions$variable <- gsub("_", "", compare_distributions$variable)

# plot the distribution of skewing per repeat, split on XCI pattern.
distribution_of_skewing_per_repeat <-
  ggplot(compare_distributions[which(!is.na(compare_distributions$value)),], 
         aes(x=variable, y=value, col=factor(`XCI pattern`, 
                                             levels = c("mosaic", "skewed", "extremely skewed")))) + 
  geom_quasirandom(dodge.width = 0.9, alpha = .5) + 
  stat_summary(fun.data = "mean_sdl", position = position_dodge(width=0.9), geom = "crossbar", width = .75)+
  theme_AL_box(legend.position="top")+
  scale_color_manual(values=colors)+
  #geom_hline(yintercept = stats_stats_melt_df_TRiX$mean, col = "black")+
  #geom_hline(yintercept = c(stats_stats_melt_df_TRiX$mean+stats_stats_melt_df_TRiX$sd), col = "pink")+
  geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+3*stats_stats_TRiXi_data$sd), col = "red")+
  #geom_hline(yintercept = c(stats_stats_melt_df_TRiX$mean+3*stats_stats_melt_df_TRiX$sd), col = "black")+
  labs(x="", y="skewing (%)")+
  scale_y_continuous(breaks = c(50, 75, 100)) +
  theme(legend.title = element_blank())+
  facet_grid(~variable=="mean", space = "free_x", scales = "free_x")

ggsave2("01_plots/Suppl_fig_1B.pdf", height = 6, width = 8,
        plot_grid(distribution_of_skewing_per_repeat))


# Show 
heehee <- TRiXi_data %>% dplyr::group_by(variable) %>% dplyr::count(value) %>% subset(is.na(value))
heehee

heehee2 <- TRiXi_data %>% dplyr::group_by(variable) %>% dplyr::count(name = "total_count")


heehee_fin <- merge(heehee, heehee2, by = c("variable"))
heehee_fin$frac <- 100-round((heehee_fin$n / heehee_fin$total_count)*100,1)

heehee_fin$variable <- gsub(heehee_fin$variable, pattern = "_skewing", replacement = "")

ggsave2(filename = "01_plots/Suppl_fig_1A.pdf",
ggplot(heehee_fin[heehee_fin$variable != "mean",], aes(x=factor(variable, levels=c("ARHGAP6", "RP2", "ZIC3", "CNKSR2", "TCAC1")), y=frac)) + geom_bar(stat="identity") + geom_text(aes(label=frac), vjust = -.15) + theme_AL_box_rotX() + labs(x="", y="fraction heterozygous (%)")
)
