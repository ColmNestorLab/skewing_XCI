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

##### Suppl. Fig. 2B, showing skewing per repeat ####
# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample) %>% rstatix::get_summary_stats(value, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "non-skewed"
stats_TRiXi_data[stats_TRiXi_data$mean >= 75,]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean >= 90,]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "non-skewed"))

# Get XCI pattern per sample (i.e. sample ID)
stats_TRiXi <- stats_TRiXi_data[,c("sample", "XCI pattern")]

# Make new df, merging the XCI pattern classification with it, Fix some variable names.
compare_distributions <- TRiXi_data[!is.na(TRiXi_data$value),]
compare_distributions <- merge(compare_distributions, stats_TRiXi, by = "sample")
compare_distributions$variable <- gsub("_skewing", "", compare_distributions$variable)
compare_distributions$variable <- gsub("_", "", compare_distributions$variable)

# Export source data
source_data_suppl_fig_2B <- compare_distributions[which(!is.na(compare_distributions$value)),]
write.table(source_data_suppl_fig_2B, "source_data/Suppl_Fig_2B_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_suppl_fig_2b <- fread("source_data/Suppl_Fig_2B_source_data.tsv")

# plot the distribution of skewing per repeat, split on XCI pattern.
suppl_fig_2b <-
  ggplot(source_data_suppl_fig_2b, 
         aes(x=variable, y=value, col=factor(`XCI pattern`, 
                                             levels = c("non-skewed", "skewed", "extremely skewed")))) + 
  geom_quasirandom(dodge.width = 0.9, alpha = .5) + 
  stat_summary(fun.data = "mean_sdl", position = position_dodge(width=0.9), geom = "crossbar", width = .75)+
  theme_AL_box(legend.position="top")+
  scale_color_manual(values=colors)+
  geom_hline(yintercept = c(90), col = "red")+
  labs(x="", y="skewing (%)")+
  scale_y_continuous(breaks = c(50, 75, 100)) +
  theme(legend.title = element_blank())+
  facet_grid(~variable=="mean", space = "free_x", scales = "free_x")

ggsave2("01_plots/Suppl_fig_2b.pdf", height = 6, width = 8,
        plot_grid(suppl_fig_2b))


# Suppl. fig. 2a.
# Show fraction heterozygous females per repeat
# Find homozygous
df_heterozygous_part1 <- TRiXi_data %>% dplyr::group_by(variable) %>% dplyr::count(value) %>% subset(is.na(value))
df_heterozygous_part2 <- TRiXi_data %>% dplyr::group_by(variable) %>% dplyr::count(name = "total_count")

# merge
df_heterozygous_fin <- merge(df_heterozygous_part1, df_heterozygous_part2, by = c("variable"))

# calculate fraction
df_heterozygous_fin$frac <- 100-round((df_heterozygous_fin$n / df_heterozygous_fin$total_count)*100,1)

# make nice variable names
df_heterozygous_fin$variable <- gsub(df_heterozygous_fin$variable, pattern = "_skewing", replacement = "")

# Export source data
source_data_suppl_fig_2a <- df_heterozygous_fin[df_heterozygous_fin$variable != "mean", c("variable", "n", "total_count", "frac")]
write.table(source_data_suppl_fig_2a, "source_data/Suppl_Fig_2A_source_data.tsv", quote = F, row.names = F, sep = "\t")

#### Read in source data. ####
# Read in.
source_data_suppl_fig_2a <- fread("source_data/Suppl_Fig_2A_source_data.tsv")

# plot
ggsave2(filename = "01_plots/Suppl_fig_2a.pdf",
ggplot(source_data_suppl_fig_2a, aes(x=factor(variable, levels=c("TRiXi1", "TRiXi2", "TRiXi3", "TRiXi4", "TRiXi5")), y=frac)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label=frac), vjust = -.15) + 
  theme_AL_box_rotX() + 
  labs(x="", y="fraction heterozygous (%)")
)
