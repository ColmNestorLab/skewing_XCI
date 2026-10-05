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

############ GTEx ############
suppl_table_10 <- fread("Supplementary_tables/supplementary_table10.tsv")

# calculate metrics for plotting.
suppl_table_10_metrics <- suppl_table_10 %>% dplyr::group_by(participant, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")
suppl_table_10_metrics$tissue_id <- factor(suppl_table_10_metrics$tissue_id)

suppl_table_10_metrics$skewed <- ifelse(suppl_table_10_metrics$median >= 0.25, yes = "skewed", no = "not_skewed")


#### Figure 3c, showing no. of skewed tissues per individual ####
countz <- suppl_table_10_metrics[suppl_table_10_metrics$tissue_id != "whole blood",] %>% dplyr::group_by(participant, skewed) %>% dplyr::count()
large_df_frac <- dcast(setDT(countz), participant ~ skewed, value.var = "n")
large_df_frac[is.na(large_df_frac$not_skewed),]$not_skewed <- 0
large_df_frac[is.na(large_df_frac$skewed),]$skewed <- 0

large_df_frac$total <- large_df_frac$skewed + large_df_frac$not_skewed
large_df_frac$fraction_skewed_tissues <- large_df_frac$skewed / large_df_frac$total
large_df_frac <- large_df_frac[order(large_df_frac$fraction_skewed_tissues, decreasing = T),]

# get metrics
paste0("total femmes", ": ", length(large_df_frac$participant))
paste0("total samples", ": ", sum(large_df_frac$not_skewed) + sum(large_df_frac$skewed))
paste0("not skewed samples", ": ", sum(large_df_frac$not_skewed))
paste0("skewed samples", ": ", sum(large_df_frac$skewed))
paste0("fraction skewed samples", ": ", sum(large_df_frac$skewed) / (sum(large_df_frac$not_skewed) + sum(large_df_frac$skewed)))

paste0("fraction skewed samples", ": ", sum(large_df_frac$skewed) / (sum(large_df_frac$not_skewed) + sum(large_df_frac$skewed)))

paste0("no. of femmes with at least one skewed tissue", ": ", length(large_df_frac[large_df_frac$skewed > 0,]$participant))
paste0("fraction femmes with at least one skewed tissue", ": ", length(large_df_frac[large_df_frac$skewed > 0,]$participant) / length(large_df_frac$participant))

# show how common 1,2,3 etc skewed tissues are per female.
large_df_frac_count <- large_df_frac %>% dplyr::group_by(skewed) %>% dplyr::count()
large_df_frac_count_total <- sum(large_df_frac_count$n)
large_df_frac_count$frac <- large_df_frac_count$n / large_df_frac_count_total
large_df_frac_count$frac2 <- round(1 / large_df_frac_count$frac, 1)

colnames(large_df_frac_count) <- c("skewed_tissue_count", "count", "frac", "frac2")

# Export source data
source_data_fig_3c <- large_df_frac_count
write.table(source_data_fig_3c, "source_data/Fig_3c_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_3c <- fread("source_data/Fig_3C_source_data.tsv")

source("plot_parameters_TRiX.R")

frequencies_per_skewing_count <-
  ggplot(source_data_fig_3c, 
         aes(x=skewed_tissue_count, y=count)) + geom_bar(stat="identity") +
  geom_text(aes(label=count)) + 
  annotate(geom="text", 
           x=source_data_fig_3c$skewed_tissue_count, 
           y=-4.5, 
           label = paste0("1/", source_data_fig_3c$frac2))+
  theme_AL_box()+
  labs(x="skewed tissue (count)", y="GTEx females (count)")

ggsave2(filename = "01_plots/Figure_3c.pdf", width = 8, height = 4,
        plot_grid(frequencies_per_skewing_count))

# how many females have one, two and three or more skewed tissues?
paste0("total = ", sum(source_data_fig_3c$count))
paste0("1 or more skewed tissue = ", sum(source_data_fig_3c[source_data_fig_3c$skewed_tissue_count >= 1,]$count))
paste0("2 or more  skewed tissue = ", sum(source_data_fig_3c[source_data_fig_3c$skewed_tissue_count >= 2,]$count))
paste0("3 or more  skewed tissue = ", sum(source_data_fig_3c[source_data_fig_3c$skewed_tissue_count >= 3,]$count))

paste0("4 or more  skewed tissue = ", sum(source_data_fig_3c[source_data_fig_3c$skewed_tissue_count >= 4,]$count))
paste0("5 or more  skewed tissue = ", sum(source_data_fig_3c[source_data_fig_3c$skewed_tissue_count >= 5,]$count))
paste0("10 or more  skewed tissue = ", sum(source_data_fig_3c[source_data_fig_3c$skewed_tissue_count >= 10,]$count))
