# Load necessary library
library(dplyr)
library(ggplot2)
library(data.table)
library(biomaRt)
library(tidyverse)
library(karyoploteR)

# good to have.
source("plot_parameters_TRiX.R")

result_df_trip <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/result_df_trip.tsv")
result_df_quad <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/result_df_quad.tsv")
result_df_auto <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/result_df_auto.tsv")

trip_smash <- result_df_trip[result_df_trip$num_repeats > 6,]
quad_smash <- result_df_quad[result_df_quad$num_repeats > 6,]
auto_smash <- result_df_auto[result_df_auto$num_repeats > 6,]

trip_smash$tag <- "trip"
quad_smash$tag <- "quad"
auto_smash$tag <- "auto"

smash_allez <- rbind(trip_smash, quad_smash, auto_smash)
selected_repeats <- c("AR", "ARHGAP6", "CNKSR2", "RP2", "TCAC1", "ZIC3", "ATXN7", "HTT", "RPL14")

smash_allez_selected <- smash_allez[smash_allez$closest_gene %in% selected_repeats,]

smash_allez_selected$name_short <- gsub(smash_allez_selected$name, pattern = "\\.R1", replacement = "")
smash_allez_selected$name_short <- gsub(smash_allez_selected$name_short, pattern = "\\.R2", replacement = "")

pdf("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/01_plots/Fig_1B.pdf", height = 10, width = 12)
ggplot(smash_allez[smash_allez$closest_gene %in% selected_repeats,], 
       aes(x=num_repeats, fill = closest_gene%in%selected_repeats)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(factor(tag, levels = c("trip", "quad", "auto"))~closest_gene, scales = "free")+
  theme_AL_box()
dev.off()

