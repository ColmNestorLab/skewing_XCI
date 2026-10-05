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
repeat_df <- fread("data/repeats_length_data_revision.tsv")

# set colnames
melt_repeat_df <- melt(data.table(repeat_df), id.vars = c("sample"), measure.vars = c("TRiXi1_allele1", "TRiXi1_allele2" ,"TRiXi2_allele1" ,"TRiXi2_allele2" ,"TRiXi3_allele1", "TRiXi3_allele2", "TRiXi4_allele1", "TRiXi4_allele2",  "TRiXi5_allele1", "TRiXi5_allele2"))

# set character
melt_repeat_df$sample <- as.character(melt_repeat_df$sample)


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

##### Suppl. Figure 1E ####
# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data[TRiXi_data$variable != "mean_skewing",] %>% dplyr::group_by(sample) %>% rstatix::get_summary_stats(value, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "non-skewed"
stats_TRiXi_data[stats_TRiXi_data$mean >= 75,]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean >= 90,]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "non-skewed"))

##### Suppl_Figure 1E, showing distribution of repeat lengths across the cohort #####
repeat_df_smash <- merge(melt_repeat_df, stats_TRiXi_data[,c("sample", "XCI pattern")], by = "sample")

# Fix names
repeat_df_smash$variable <- str_sub(repeat_df_smash$variable, end = -9)

# Export source data
# remove sample column
source_data_fig_2c <- repeat_df_smash[,-1]
write.table(source_data_fig_2c, "source_data/Suppl_Fig_2c_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_2c <- fread("source_data/Suppl_Fig_2c_source_data.tsv")

source_data_fig_2c$value <- as.numeric(source_data_fig_2c$value)

# Perform kruskal wallis test
kruskal_results <- source_data_fig_2c %>%
  group_by(variable) %>%
  kruskal_test(value ~ `XCI pattern`)

# plot
repeat_distribution <- 
  gghistogram(data=source_data_fig_2c, x="value", fill = "XCI pattern", binwidth = 1) + 
  facet_grid(`XCI pattern`~variable, scales = "free") + 
  theme_AL_box(legend.position="none") + 
  labs(x="repeat length") +
  geom_text(data = kruskal_results, 
            aes(x = Inf, y = Inf, label = paste("K-S p-value:\n", format(p, digits = 4))), 
            hjust = 1.1, vjust = 1.5, size = 3, color = "black") +
  scale_fill_manual(values=colors)


ggsave2("01_plots/Suppl_fig_2c.pdf", height = 4.5, width = 7.5,
        plot_grid(repeat_distribution))
