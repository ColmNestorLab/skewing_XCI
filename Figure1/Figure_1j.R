library(data.table)
library(tidyverse)
library(ggplot2)

source("plot_parameters_TRiX.R")

AML_samples <- melt(fread("data/revision_new_experiment/AML_samples/AML_samples_TRiXi.tsv"), id.vars = "sample", measure.vars = c("TRiXi1", "TRiXi2", "TRiXi3", "TRiXi4", "TRiXi5"))

AML_samples$value <- as.numeric(AML_samples$value) 

# export source data
write.table(AML_samples, "source_data/figure_1J.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# read in source data
figure_1J_source_data <- fread("source_data/figure_1J.tsv")

ggsave2(filename = "01_plots/AML_samples.pdf",
ggplot(figure_1J_source_data, aes(x=variable, y=value)) + 
  geom_bar(stat="identity") + 
  theme_AL_box_rotX() + 
  facet_grid(~sample, space = "free", scale = "free") + 
  geom_hline(yintercept = 90, lty=2, col = "red")
)
