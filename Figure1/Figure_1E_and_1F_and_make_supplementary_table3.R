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
#df_miniTRiX <- fread("data/minitrix_results.tsv")
#df_TRiX <- fread("data/trix_results.tsv")
#df_HUMARA <- fread("data/HUMARA_results.tsv")

#df_miniTRiX$method <- "mini"
#df_TRiX$method <- "full"

#melt_df_miniTRiX <- melt(data.table(df_miniTRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing"))
#melt_df_miniTRiX$value <- as.numeric(melt_df_miniTRiX$value)

#melt_df_TRiX <- melt(data.table(df_TRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing", "CNKSR2_skewing", "TCAC1_skewing", "ZIC3_skewing"))
#melt_df_TRiX$value <- as.numeric(melt_df_TRiX$value)

# rbind
#TRiXi_data <- rbind(melt_df_miniTRiX, melt_df_TRiX)

# get TRiX runs for the HUMARA matched samples.
#common_elements <- intersect(df_TRiX$sample, df_HUMARA$sample)
#missing_elements <- setdiff(df_HUMARA$sample, df_TRiX$sample)

#df_TRiX_HUMARA_samples <- df_TRiX[df_TRiX$sample %in% common_elements,]
#df_HUMARA_filt <- df_HUMARA[df_HUMARA$sample %in% common_elements,]

# Export supplementary table 3
#write.table(merge(df_TRiX_HUMARA_samples, df_HUMARA_filt, by = "sample"), row.names = F, quote = F, sep = "\t", file = "supplementary_table3.tsv")


supplementary_table3 <- fread("Supplementary_tables/supplementary_table3.tsv")[-c(1:2),]
colnames(supplementary_table3) <- c("sample", "TRiXi_mean_skewing", "ARHGAP6_skewing", "CNKSR2_skewing", "TCAC1_skewing", "RP2_skewing", "ZIC3_skewing", "AR_skewing")


# make compatible by removing excess columns
trix_humara <- supplementary_table3[,c("sample", "TRiXi_mean_skewing")]
humara_humara <- supplementary_table3[,c("sample", "AR_skewing")]

# add method
trix_humara$method <- "TRiX"
humara_humara$method <- "HUMARA"

# change colnames
colnames(trix_humara) <- c("sample", "skewing", "method")
colnames(humara_humara) <- c("sample", "skewing", "method")

# rbind the two
humara_comparison <- rbind(trix_humara, humara_humara)

# add tag
humara_comparison$tag <- ifelse(is.na(humara_comparison$skewing), yes = "homozygous samples", no = "heterozygous samples")

# pivot wider
df_wide <- humara_comparison[,c("sample", "skewing",   "method")] %>%  pivot_wider(names_from = method, values_from = skewing)

# Export source data
source_data <- df_wide
write.table(source_data, "source_data/Fig_1E_and_Fig_1F_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data <- fread("source_data/Fig_1E_and_Fig_1F_source_data.tsv")

source_data$HUMARA <- as.numeric(source_data$HUMARA)
source_data$TRiX <- as.numeric(source_data$TRiX)

## Figure 1E ##
# Scatter plot with regression line & R² value
HUMARA_comparison_correlation <- 
  plot_grid(ncol=1,
            ggscatter(source_data, x = "TRiX", y = "HUMARA",
                      add = "reg.line",        # Add regression line
                      conf.int = TRUE,         # Confidence interval
                      cor.coef = TRUE,         # Show correlation coefficient (R)
                      cor.method = "spearman") + geom_point(aes(col=sample=="320280"), show.legend = F)
  )

# find samples with homozygous AR
homo_humara <- source_data[is.na(source_data$HUMARA),]$sample

source_data_long <- melt(source_data, value.name = "skewing", id.vars = "sample", variable.name = "method")

# add tag
source_data_long$homo_humara <- ifelse(source_data_long$sample %in% homo_humara, yes = "homozygous in HUMARA", no = NA)

humara_comparison_boxplot <- 
  ggplot(source_data_long,
         aes(x=method, y=skewing))+
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(aes( col = homo_humara), show.legend = F)+
  geom_quasirandom(data=source_data_long[source_data_long$sample == "320280",], col= "blue") +
  geom_text_repel(data=source_data_long[source_data_long$sample == "320280",], aes(label=sample)) +
  theme_AL_box_rotX() +
  stat_compare_means(label.x = 1.5)

# merge plots
humara_comparison_plot <-
  plot_grid(ncol=2,
            HUMARA_comparison_correlation,
            humara_comparison_boxplot)

# Plot
ggsave2(filename = "01_plots/Figure_1E_and_1F.pdf", units = "mm", #width = 400, height = 700,
        plot_grid(humara_comparison_plot))



