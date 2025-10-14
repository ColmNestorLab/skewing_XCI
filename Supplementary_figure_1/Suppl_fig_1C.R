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


##### Suppl_Figure 1B, showing distribution of repeat lengths across the cohort #####

repeat_df_mini <- melt(df_miniTRiX, id.vars = c("sample", "method"), measure.vars = c("ARHGAP6_allele1", "ARHGAP6_allele2", "RP2_allele1", "RP2_allele2"), variable.name = "repeat_name")

repeat_df_TRiXi <- melt(df_TRiX, id.vars = c("sample", "method"), measure.vars = c("ARHGAP6_allele1", "ARHGAP6_allele2", "RP2_allele1", "RP2_allele2", "CNKSR2_allele1", "CNKSR2_allele2","TCAC1_allele1", "TCAC1_allele2", "ZIC3_allele1", "ZIC3_allele2"), variable.name = "repeat_name")

TRiXi_repeat_data <- rbind(repeat_df_mini,
                           repeat_df_TRiXi)


melt_df_allele_sizes <- merge(rbind(repeat_df_mini,
                                    repeat_df_TRiXi), 
                              stats_TRiXi_data[,c("sample", "XCI pattern")], 
                              by = "sample")


# Fix names
melt_df_allele_sizes$repeat_name <- str_sub(melt_df_allele_sizes$repeat_name, end = -8)
melt_df_allele_sizes$repeat_name <- gsub("_", "", melt_df_allele_sizes$repeat_name)

# set numeric
melt_df_allele_sizes$value <- as.numeric(melt_df_allele_sizes$value)

# Perform kruskal wallis test
kruskal_results <- melt_df_allele_sizes %>%
  group_by(repeat_name) %>%
  kruskal_test(value ~ `XCI pattern`)

countz <- unique(melt_df_allele_sizes[,c("sample","XCI pattern")]) %>% dplyr::group_by(`XCI pattern`) %>% dplyr::count()

# plot
repeat_distribution <- 
  gghistogram(data=melt_df_allele_sizes, x="value", fill = "XCI pattern", binwidth = 1) + 
  facet_grid(`XCI pattern`~repeat_name, scales = "free") + 
  theme_AL_box(legend.position="none") + 
  labs(x="repeat length", title = "extremely skewed_9, skewed_44, mosaic_972"  ) +
  geom_text(data = kruskal_results, 
            aes(x = Inf, y = Inf, label = paste("K-S p-value:\n", format(p, digits = 4))), 
            hjust = 1.1, vjust = 1.5, size = 3, color = "black") +
  scale_fill_manual(values=colors)


ggsave2("01_plots/Suppl_figure_1C.pdf", height = 4.5, width = 7.5,
        plot_grid(repeat_distribution))