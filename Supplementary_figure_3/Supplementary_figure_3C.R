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
suppl_table_10_metrics <- suppl_table_10 %>% dplyr::group_by(participant, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")
suppl_table_10_metrics$tissue_id <- factor(suppl_table_10_metrics$tissue_id)

# get metrics
metricx <- suppl_table_10_metrics %>%
  mutate(
    skewing = 50 + median * 100
  ) %>%
  mutate(
    skew_group = case_when(
      skewing >= 75 ~ "skewed",
      TRUE ~ "not_skewed"
    )
  )

# count
metricx_tissues <- metricx %>% dplyr::group_by(tissue_id, skew_group) %>% dplyr::count()

# keep only skewed
metricx_tissues_skew <- metricx_tissues[metricx_tissues$skew_group == "skewed",]

# count fraction of skewed tissues that are whole blood (8.2%)
sum(metricx_tissues_skew[metricx_tissues_skew$tissue_id == "whole blood",]$n) / (sum(metricx_tissues_skew$n))

# no. of females with skewed tissues, not including whole blood (207).
length(unique(metricx[metricx$skew_group == "skewed" & metricx$tissue_id != "whole blood",]$participant))


# plot XCI skewing variation per individual
variation_df <- suppl_table_10_metrics[suppl_table_10_metrics$tissue_id != "whole blood",] %>%
  mutate(
    skewing = 50 + median * 100
  ) %>%
  group_by(participant) %>%
  summarise(
    median_skewing = median(skewing),
    mad_skewing = mad(skewing, constant = 1),
    range_skewing = max(skewing) - min(skewing),
    n_tissues = n(),
    .groups = "drop"
  )

# add skewing threhsolds
variation_df <- variation_df %>%
  mutate(
    skew_group = case_when(
      median_skewing >= 90 ~ "extreme",
      median_skewing >= 75 ~ "skewed",
      TRUE ~ "not_skewed"
    )
  )

# Export source data
Sfig_3c <- variation_df
write.table(Sfig_3c, "source_data/Suppl_fig_3c_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
Sfig_3c <- fread("source_data/Suppl_fig_3c_source_data.tsv")


# plot
stat.test <- Sfig_3c %>%
  dunn_test(mad_skewing ~ skew_group, p.adjust.method = "holm") %>%
  add_y_position()

stat.test$p.adj <- round(stat.test$p.adj, 5)


ggsave2(filename = "01_plots/Suppl_fig_3c.pdf", width = 6, height = 4,
  ggplot(Sfig_3c,
               aes(x=factor(skew_group,levels=c("not_skewed", "skewed", "extreme")), y=mad_skewing, col = factor(skew_group,levels=c("not_skewed", "skewed", "extreme")))) +
          geom_boxplot(outlier.shape = NA) +
          geom_quasirandom() +
          theme_AL_box_rotX(legend.position="none") +
          labs(x="", y="Tissue-to-tissue XCI skewing variability (MAD)") + 
          stat_compare_means(label.y = 17) +
          stat_pvalue_manual(
            stat.test,
            label = "p.adj"
          )+
          coord_cartesian(ylim=c(0,20))
        )
