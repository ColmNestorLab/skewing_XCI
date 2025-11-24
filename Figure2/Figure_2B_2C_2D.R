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

# also read in the longitudinal data
df_miniTRiX_age <- fread("Supplementary_tables/supplementary_table6.tsv")[-1,]

colnames(df_miniTRiX_age) <- as.character(df_miniTRiX_age[1,])
df_miniTRiX_age <- df_miniTRiX_age[-1,]

#melt
melt_df_miniTRiX_age <- melt(df_miniTRiX_age, id.vars = "sample_ID", measure.vars = c("mean_skewing_age_0_to_1_year", "mean_skewing_age_14_to_16_years"))

# Classify XCI pattern of individuals.
melt_df_miniTRiX_age$'XCI pattern' <- "mosaic"
melt_df_miniTRiX_age[melt_df_miniTRiX_age$value > stats_stats_TRiXi_data$mean+stats_stats_TRiXi_data$sd,]$'XCI pattern' <- "skewed"
melt_df_miniTRiX_age[melt_df_miniTRiX_age$value > (stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "extremely skewed"
melt_df_miniTRiX_age$'XCI pattern' <- factor(melt_df_miniTRiX_age$'XCI pattern', levels = c("extremely skewed", "skewed", "mosaic"))


# make new df
df_for_alluvial <- melt_df_miniTRiX_age[,c("sample_ID", "variable", "XCI pattern")]

# make sure its all factors
df_for_alluvial$sample_ID <- factor(df_for_alluvial$sample_ID)
df_for_alluvial$variable <- factor(df_for_alluvial$variable)
df_for_alluvial$`XCI pattern` <- factor(df_for_alluvial$`XCI pattern`)

# change rows a bit
df_for_alluvial$variable <- gsub(df_for_alluvial$variable, pattern= "mean_skewing_age_0_to_1_year", replacement="1")
df_for_alluvial$variable <- gsub(df_for_alluvial$variable, pattern= "mean_skewing_age_14_to_16_years", replacement="2")



# Combine XCI pattern with sample identifier to ensure uniqueness
df_spread <- df_for_alluvial %>%
  dplyr::mutate(Sample = paste0("Sample_", variable)) %>%
  dplyr::select(sample_ID, Sample, `XCI pattern`) %>%
  spread(key = Sample, value = `XCI pattern`) %>%
  dplyr::mutate(Sample_1 = paste("Sample 1", Sample_1),
                Sample_2 = paste("Sample 2", Sample_2))

# spread data
df_spread <- df_spread %>%
  mutate(Sample_1 = factor(Sample_1, levels = c("Sample 1 extremely skewed", "Sample 1 skewed", "Sample 1 mosaic")),
         Sample_2 = factor(Sample_2, levels = c("Sample 2 extremely skewed", "Sample 2 skewed", "Sample 2 mosaic")))

# make df for plotting XCI pattern per age
df_spread_count <- df_for_alluvial %>% dplyr::group_by(variable, `XCI pattern`) %>% dplyr::count()

# Calculate the total for each category
df_spread_count <- df_spread_count %>%
  dplyr::group_by(variable) %>%
  mutate(total = sum(n),
         fraction = paste0(round((n / total)*100, 2), "%"))  # Calculate fractions

# Export source data
source_data_fig_2C <- df_spread_count
write.table(source_data_fig_2C, "source_data/Fig_2C_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_2C <- fread("source_data/Fig_2C_source_data.tsv")

# get statistics
# Create a contingency table
contingency_table_age <- xtabs(n ~ variable + `XCI pattern`, data = source_data_fig_2C)

# Perform the Chi-squared test
chi_squared_result_age <- chisq.test(contingency_table_age)

# Perform Fisher's exact test
fisher_result_age <- fisher.test(contingency_table_age, workspace = 2e8)

# plot
bar_distribution_of_XCI_patterns_across_age <-
  ggplot(source_data_fig_2C, aes(x=variable, y=n, fill=`XCI pattern`, label = fraction)) + 
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=colors) +
  geom_text(position = "fill", vjust = 1.25, col = "white")+
  theme_AL_box() +
  theme(plot.subtitle = element_text(size=8), legend.position = "top")+
  labs(x="age measurement", y="fraction (%)", subtitle = paste0("Chi-squared pval = ",round(chi_squared_result_age$p.value, 3),"\n", "Fisher's exact test pval = ", round(fisher_result_age$p.value, 3)))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),  expand = expansion(mult = c(0, 0.025)))

# Export age plots
ggsave2(filename = "01_plots/Figure_2C.pdf", units = "mm", #width = 400, height = 700,
        plot_grid(bar_distribution_of_XCI_patterns_across_age)
)



# Set order
HM_order <- unique(c(c(1022, 1081, 229, 1676, 864, 4308, 212, 315, 511, 1061, 5786), 
                     setdiff(df_for_alluvial$sample_ID, c(1022, 1081, 229, 1676, 864, 4308, 212, 315, 511, 1061, 5786))))


# add difference column. set factor.
df_miniTRiX_age$change <- as.numeric(df_miniTRiX_age$mean_skewing_age_14_to_16_years) - as.numeric(df_miniTRiX_age$mean_skewing_age_0_to_1_year)
df_miniTRiX_age$sample_ID <- factor(df_miniTRiX_age$sample_ID)

# make new df
df_Figure_2D <- df_for_alluvial

#set factor.
df_Figure_2D$sample_ID <- factor(df_Figure_2D$sample_ID)

# merge
df_Figure_2D <- merge(df_Figure_2D, df_miniTRiX_age[,c("change", "sample_ID"),], by.x = "sample_ID", by.y = "sample_ID")

# Export source data
source_data_fig_2D <- df_Figure_2D[df_Figure_2D$sample_ID %in% c(1022, 1081, 229, 1676, 864, 4308, 212, 315, 511, 1061, 5786),]
write.table(source_data_fig_2D, "source_data/Fig_2D_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_2D <- fread("source_data/Fig_2D_source_data.tsv")

# plot figure 2D
ggsave2(filename = "01_plots/Figure_2D.pdf", height = 4, width = 4,
ggplot(source_data_fig_2D, 
       aes(y=factor(sample_ID, levels = rev(HM_order)), x=variable, fill = `XCI pattern`)) + 
  geom_tile(color="white") + 
  theme_AL_box(legend.position = "none") +   
  geom_text(aes(label = round(change, 2)), color = "black") + 
  coord_fixed(1/1) + labs(x="", y="")
)


# Export source data
source_data_fig_2B <- melt_df_miniTRiX_age
write.table(source_data_fig_2B, "source_data/Fig_2B_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_fig_2B <- fread("source_data/Fig_2B_source_data.tsv")

# Make order
ordrrr2 <- source_data_fig_2B[source_data_fig_2B$variable == "mean_skewing_age_0_to_1_year",]
ordrrr2 <- ordrrr2[order(ordrrr2$value, decreasing = F),]
ordrrr2 <- ordrrr2$sample_ID

# plot figure 2B
ggsave2(filename = "01_plots/Figure_2B.pdf", height = 4, width = 10,
ggplot(source_data_fig_2B, 
       aes(x=factor(sample_ID, levels=ordrrr2), y=value, col = variable)) + 
  geom_point() + 
  theme_AL_box_rotX() + 
  geom_line(aes(group=sample_ID), col = "grey", lty=2) + 
  geom_hline(yintercept = c(62.4,71.09, 79.76), col = c("black","pink", "red"), lty=2) + 
  scale_color_manual(values=c("mean_skewing_age_0_to_1_year"="black","mean_skewing_age_14_to_16_years"="red"))+
  labs(x="", y="skewing (%)") +
   coord_cartesian(ylim=c(50,90))
)

