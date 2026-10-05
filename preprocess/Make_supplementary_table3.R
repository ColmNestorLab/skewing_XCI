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
# make supplementary table 3
df_miniTRiX <- fread("data/minitrix_results.tsv")
df_TRiX <- fread("data/trix_results.tsv")
df_HUMARA <- fread("data/HUMARA_results.tsv")
df_miniTRiX$method <- "mini"
df_TRiX$method <- "full"

melt_df_miniTRiX <- melt(data.table(df_miniTRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing"))
melt_df_miniTRiX$value <- as.numeric(melt_df_miniTRiX$value)
melt_df_TRiX <- melt(data.table(df_TRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing", "CNKSR2_skewing", "TCAC1_skewing", "ZIC3_skewing"))
melt_df_TRiX$value <- as.numeric(melt_df_TRiX$value)

# rbind
TRiXi_data <- rbind(melt_df_miniTRiX, melt_df_TRiX)

# get TRiX runs for the HUMARA matched samples.
common_elements <- intersect(df_TRiX$sample, df_HUMARA$sample)
missing_elements <- setdiff(df_HUMARA$sample, df_TRiX$sample)

df_TRiX_HUMARA_samples <- df_TRiX[df_TRiX$sample %in% common_elements,]
df_HUMARA_filt <- df_HUMARA[df_HUMARA$sample %in% common_elements,]

# Export supplementary table 3
write.table(merge(df_TRiX_HUMARA_samples, df_HUMARA_filt, by = "sample"), row.names = F, quote = F, sep = "\t", file = "supplementary_table3.tsv")
