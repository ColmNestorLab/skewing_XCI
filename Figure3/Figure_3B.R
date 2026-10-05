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

#####  TRIXI  ###########
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

# keep only whole blood
suppl_table_10_WB <- suppl_table_10_metrics[suppl_table_10_metrics$tissue_id == "whole blood",c("participant","variable","median")]

# make TRiXi data frame thinner
stats_TRiXi_data_short <- stats_TRiXi_data[,c("sample", "variable",  "median")]

# change colnames to match
colnames(stats_TRiXi_data_short) <- c("sample", "variable",  "skewing")
colnames(suppl_table_10_WB) <- c("sample", "variable",  "skewing")

# rbind
smash_dat <- rbind(stats_TRiXi_data_short, suppl_table_10_WB)

smash_dat$variable <- gsub(smash_dat$variable, pattern = "value", replacement = "TRiXi - whole blood")
smash_dat$variable <- gsub(smash_dat$variable, pattern = "allelic_expression", replacement = "GTEx - whole blood")

# rescale function
scale_range <- function(x, new_min = 50, new_max = 100, old_min = 0, old_max = 0.5) {
  (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
}

# rename & rescale
test_scale <- smash_dat[smash_dat$variable %in% c("GTEx - whole blood"),]
test_scale$skewing <- scale_range(test_scale$skewing)

# rename & add skewed tag
fin_scaled <- rbind(test_scale, smash_dat[smash_dat$variable %in% c("TRiXi - whole blood"),])
fin_scaled$skewed <- ifelse(fin_scaled$skewing >= 90, yes = "skewed", no = "not_skewed")

# Export source data
source_data_fig_3b <- fin_scaled
write.table(source_data_fig_3b, "source_data/Fig_3B_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_3b <- fread("source_data/Fig_3b_source_data.tsv")

# calculate fraction skewed (to add to plot later)
fin_scaled_count_tmp <- source_data_fig_3b %>% dplyr::group_by(variable, skewed) %>% dplyr::count()
fin_scaled_count <- fin_scaled_count_tmp %>%
  group_by(variable) %>%
  mutate(pct = n / sum(n))


# plot density
pdf("01_plots/Figure_3b.pdf", height = 5)
plot_grid(rel_widths = c(0.9,0,1), ncol = 2,
          ggdensity(source_data_fig_3b, x="skewing", col = "variable", add = "median") + 
            geom_vline(xintercept = 90, lty=2) + 
            annotate(geom="text", x=c(95,95), y=c(0.0125,0.0015), 
                     label= round(fin_scaled_count[fin_scaled_count$skewed == "skewed",]$pct, 3)*100)
)
dev.off()
