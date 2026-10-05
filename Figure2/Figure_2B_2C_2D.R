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

##### Figure 2B ####
# Read in the longitudinal data
df_miniTRiX_age <- fread("Supplementary_tables/supplementary_table8.tsv")[-1,]

colnames(df_miniTRiX_age) <- as.character(df_miniTRiX_age[1,])
df_miniTRiX_age <- df_miniTRiX_age[-1,]

#melt
melt_df_miniTRiX_age <- melt(df_miniTRiX_age, id.vars = "sample_ID", measure.vars = c("mean_skewing_age_0_to_1_year", "mean_skewing_age_14_to_16_years"))

# set numeric
melt_df_miniTRiX_age$value <- as.numeric(melt_df_miniTRiX_age$value)

# Classify XCI pattern of individuals.
melt_df_miniTRiX_age$'XCI pattern' <- "non-skewed"
melt_df_miniTRiX_age[melt_df_miniTRiX_age$value >= 75,]$'XCI pattern' <- "skewed"
melt_df_miniTRiX_age[melt_df_miniTRiX_age$value >= 90 ,]$'XCI pattern' <- "extremely skewed"
melt_df_miniTRiX_age$'XCI pattern' <- factor(melt_df_miniTRiX_age$'XCI pattern', levels = c("extremely skewed", "skewed", "non-skewed"))


#### Plot 2B ####
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
          geom_hline(yintercept = c(75, 90), col = c("pink", "red"), lty=2) + 
          scale_color_manual(values=c("mean_skewing_age_0_to_1_year"="black","mean_skewing_age_14_to_16_years"="red"))+
          labs(x="", y="skewing (%)") +
          coord_cartesian(ylim=c(50,90))
)


# plot 2D
# make new df
df_Figure_2D_preprocess <- source_data_fig_2B[,c("sample_ID", "variable", "XCI pattern")]

# make sure its all factors
df_Figure_2D_preprocess$sample_ID <- factor(df_Figure_2D_preprocess$sample_ID)
df_Figure_2D_preprocess$variable <- factor(df_Figure_2D_preprocess$variable)
df_Figure_2D_preprocess$`XCI pattern` <- factor(df_Figure_2D_preprocess$`XCI pattern`)

# change rows a bit
df_Figure_2D_preprocess$variable <- gsub(df_Figure_2D_preprocess$variable, pattern= "mean_skewing_age_0_to_1_year", replacement="1")
df_Figure_2D_preprocess$variable <- gsub(df_Figure_2D_preprocess$variable, pattern= "mean_skewing_age_14_to_16_years", replacement="2")

# add difference column. set factor.
df_miniTRiX_age$change <- as.numeric(df_miniTRiX_age$mean_skewing_age_14_to_16_years) - as.numeric(df_miniTRiX_age$mean_skewing_age_0_to_1_year)
df_miniTRiX_age$sample_ID <- factor(df_miniTRiX_age$sample_ID)

# make new df
df_Figure_2D <- df_Figure_2D_preprocess

#set factor.
df_Figure_2D$sample_ID <- factor(df_Figure_2D$sample_ID)

# merge
df_Figure_2D <- merge(df_Figure_2D, df_miniTRiX_age[,c("change", "sample_ID", "mean_skewing_age_0_to_1_year", "mean_skewing_age_14_to_16_years"),], by.x = "sample_ID", by.y = "sample_ID")

# Here I select samples that change 5% skewing or more and/or has skewing >=75 at any timepoint.
large_change <- as.character(unique(df_Figure_2D[(df_Figure_2D$change <= -5 | df_Figure_2D$change >= 5 | df_Figure_2D$`XCI pattern` == "skewed"),]$sample_ID))
source_data_fig_2C <- df_Figure_2D
source_data_fig_2D <- df_Figure_2D[df_Figure_2D$sample_ID %in% large_change,]

# export
write.table(source_data_fig_2C, "source_data/Fig_2C_source_data.tsv", quote = F, row.names = F, sep = "\t")
write.table(source_data_fig_2D, "source_data/Fig_2D_source_data.tsv", quote = F, row.names = F, sep = "\t")

# Read in.
source_data_fig_2C <- fread("source_data/Fig_2C_source_data.tsv")

ggsave2(filename = "01_plots/Figure_2C.pdf", height = 4, width = 4,
        ggplot(source_data_fig_2C[source_data_fig_2C$variable == "2",], 
         aes(x=reorder(sample_ID, change), 
             y=change, 
             fill = sample_ID%in%source_data_fig_2C[source_data_fig_2C$`XCI pattern` == "skewed",]$sample_ID)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept = c(-5,5), lty=2, col="red") + 
  theme_AL_box_rotX() + 
  coord_cartesian(ylim=c(-16,16)) + 
  scale_y_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_hline(yintercept = median(source_data_fig_2C[source_data_fig_2C$variable == "1",]$change), lty=1)+
  theme(legend.title = element_blank())
)

# Read in.
source_data_fig_2D <- fread("source_data/Fig_2D_source_data.tsv")

# Set order
HM_order <- unique(c(c(1022, 1081, 229, 1676, 864, 4308, 212, 315, 511, 1061, 5786), 
                     setdiff(df_Figure_2D_preprocess$sample_ID, c(1022, 1081, 229, 1676, 864, 4308, 212, 315, 511, 1061, 5786))))

# plot figure 2D
ggsave2(filename = "01_plots/Figure_2D.pdf", height = 4, width = 4,
        ggplot(source_data_fig_2D, 
               aes(y=factor(sample_ID, levels = rev(HM_order)), x=variable, fill = `XCI pattern`)) + 
          geom_tile(color="white") + 
          theme_AL_box(legend.position = "top") +   
          geom_text(aes(label = round(change, 2)), color = "black") + 
          coord_fixed(1/1) + labs(x="", y="")
)


