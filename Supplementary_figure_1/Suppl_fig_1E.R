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



# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample) %>% rstatix::get_summary_stats(value, type = "common")

# calculate descriptive statistics of descriptive statistics, to plot lines mean, and 1sd and 2sd.
stats_stats_TRiXi_data <- stats_TRiXi_data %>% rstatix::get_summary_stats(mean, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "mosaic"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+1*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "mosaic"))

##### SFigure 1E, showing skewing per repeat for the completely skewed females ####

# Get XCI pattern per sample (i.e. sample ID)
stats_TRiXi <- stats_TRiXi_data[,c("sample", "XCI pattern")]

# Make new df, merging the XCI pattern classification with it, Fix some variable names.
compare_distributions <- TRiXi_data[!is.na(TRiXi_data$value)]
compare_distributions <- merge(compare_distributions, stats_TRiXi, by = "sample")
compare_distributions$variable <- gsub("_skewing", "", compare_distributions$variable)
compare_distributions$variable <- gsub("_", "", compare_distributions$variable)


# Plot extremely skewed individuals for each repeat
extremely_skewed_femmes <- compare_distributions[compare_distributions$`XCI pattern` == "extremely skewed",]

# Export source data
source_data_suppl_fig_1E <- extremely_skewed_femmes[extremely_skewed_femmes$variable != "mean",]
write.table(source_data_suppl_fig_1E, "source_data/Suppl_Fig_1E_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_suppl_fig_1E <- fread("source_data/Suppl_Fig_1E_source_data.tsv")

eXCI_femme_per_repeat <- 
  plot_grid(rel_widths = c(1,0.5),
            ggplot(source_data_suppl_fig_1E, 
                   aes(x=reorder(factor(sample), -value, mean), y=value)) + 
              geom_point(aes(col=variable)) + 
              stat_summary(fun.data = mean_se, geom = "crossbar") +
              geom_hline(yintercept = c(stats_stats_TRiXi_data$mean+3*stats_stats_TRiXi_data$sd), 
                         col = "red")+
              theme_AL_box_rotX() +
              labs(x="extremely skewed XCI samples", y="skewing(%)")
            ,
            ggplot(source_data_suppl_fig_1E, 
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


ggsave2("01_plots/Suppl_fig_1E.pdf", height = 4, width = 14,
        plot_grid(eXCI_femme_per_repeat))
