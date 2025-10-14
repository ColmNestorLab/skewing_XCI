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



# get TRiX runs for the HUMARA matched samples.
common_elements <- intersect(df_TRiX$sample, df_HUMARA$sample)
missing_elements <- setdiff(df_HUMARA$sample, df_TRiX$sample)

df_TRiX_HUMARA_samples <- df_TRiX[df_TRiX$sample %in% common_elements,]
df_HUMARA_filt <- df_HUMARA[df_HUMARA$sample %in% common_elements,]

# Export supplementary table 3
write.table(merge(df_TRiX_HUMARA_samples, df_HUMARA_filt, by = "sample"), row.names = F, quote = F, sep = "\t", file = "supplementary_table3.tsv")



trix_humara <- df_TRiX_HUMARA_samples[,c("sample", "mean_skewing")]
humara_humara <- df_HUMARA_filt[,c("sample", "AR_skewing")]

trix_humara$method <- "TRiX"
humara_humara$method <- "HUMARA"

colnames(trix_humara) <- c("sample", "skewing", "method")
colnames(humara_humara) <- c("sample", "skewing", "method")

humara_comparison <- rbind(trix_humara, humara_humara)

humara_comparison$tag <- ifelse(is.na(humara_comparison$skewing), yes = "homozygous samples", no = "heterozygous samples")

humara_comparison_order <- humara_comparison %>% dplyr::group_by(sample) %>% dplyr::summarise(meanz = median(skewing, na.rm = T))
humara_comparison_order <- humara_comparison_order[order(humara_comparison_order$meanz, decreasing = F),]

# calculate XCI skew classification with HUMARA
df_HUMARA_filt_stats <- df_HUMARA_filt[!is.na(df_HUMARA_filt$AR_skewing),] %>% rstatix::get_summary_stats(AR_skewing, type = "common")

df_TRiX_HUMARA_samples_stats <- df_TRiX_HUMARA_samples[!is.na(df_TRiX_HUMARA_samples$mean_skewing),] %>% rstatix::get_summary_stats(mean_skewing, type = "common")

df_wide <- humara_comparison[,c("sample", "skewing",   "method")] %>%
  pivot_wider(names_from = method, values_from = skewing)


## Figure 1E ##
# Scatter plot with regression line & R² value
HUMARA_comparison_correlation <- 
  plot_grid(ncol=1,
            ggscatter(df_wide, x = "TRiX", y = "HUMARA",
                      add = "reg.line",        # Add regression line
                      conf.int = TRUE,         # Confidence interval
                      cor.coef = TRUE,         # Show correlation coefficient (R)
                      cor.method = "spearman") + geom_point(aes(col=sample=="320280"), show.legend = F)
  )



homo_humara <- df_wide[is.na(df_wide$HUMARA),]$sample

humara_comparison$homo_humara <- ifelse(humara_comparison$sample %in% homo_humara, yes = "homozygous in HUMARA", no = NA)

humara_comparison_boxplot <- 
  ggplot(humara_comparison,
         aes(x=method, y=skewing))+
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(aes( col = homo_humara), show.legend = F)+
  geom_quasirandom(data=humara_comparison[humara_comparison$sample == "320280",], col= "blue") +
  geom_text_repel(data=humara_comparison[humara_comparison$sample == "320280",], aes(label=sample)) +
  theme_AL_box_rotX() +
  stat_compare_means(label.x = 1.5)

humara_comparison_plot <-
  plot_grid(ncol=2,
            HUMARA_comparison_correlation,
            humara_comparison_boxplot)


ggsave2(filename = "01_plots/Figure_1E_and_1F.pdf", units = "mm", #width = 400, height = 700,
        plot_grid(humara_comparison_plot))



