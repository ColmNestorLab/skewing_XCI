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



# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample, method) %>% rstatix::get_summary_stats(value, type = "common")

# calculate descriptive statistics of descriptive statistics, to plot lines mean, and 1sd and 2sd.
stats_stats_TRiXi_data <- stats_TRiXi_data %>% rstatix::get_summary_stats(mean, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "mosaic"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+1*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "mosaic"))

##### SFigure 1C, showing skewing per repeat for the completely skewed females ####

# Get XCI pattern per sample (i.e. sample ID)
stats_TRiXi <- stats_TRiXi_data[,c("sample", "XCI pattern")]

# Make new df, merging the XCI pattern classification with it, Fix some variable names.
compare_distributions <- TRiXi_data[!is.na(TRiXi_data$value)]
compare_distributions <- merge(compare_distributions, stats_TRiXi, by = "sample")
compare_distributions$variable <- gsub("_skewing", "", compare_distributions$variable)
compare_distributions$variable <- gsub("_", "", compare_distributions$variable)


# Plot extremely skewed individuals for each repeat
extremely_skewed_femmes <- compare_distributions[compare_distributions$`XCI pattern` == "extremely skewed",]

eXCI_femme_per_repeat <- 
  plot_grid(rel_widths = c(1,0.5),
            ggplot(extremely_skewed_femmes[extremely_skewed_femmes$variable != "mean",], 
                   aes(x=reorder(factor(sample), -value, mean), y=value)) + 
              geom_point(aes(col=variable)) + 
              stat_summary(fun.data = mean_se, geom = "crossbar") +
              geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+3*stats_stats_TRiXi_data$sd), 
                         col = "red")+
              theme_AL_box_rotX() +
              labs(x="extremely skewed XCI samples", y="skewing(%)")
            ,
            ggplot(extremely_skewed_femmes[extremely_skewed_femmes$variable != "mean",], 
                   aes(x=reorder(variable, value, mean), y=value)) + 
              geom_quasirandom(aes(col=variable), width = .35) + 
              stat_summary(fun.data = mean_se, geom = "crossbar") +
              geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+3*stats_stats_TRiXi_data$sd), 
                         col = "red")+
              theme_AL_box_rotX() +
              labs(x="", y="skewing(%)") +
              stat_compare_means(method="anova")+
              theme(legend.title = element_blank())
  )


ggsave2("01_plots/Suppl_fig_1C.pdf", height = 4, width = 14,
        plot_grid(eXCI_femme_per_repeat))
