library(seqinr)
library(ggplot2)
library(tidyverse)

source("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/plot_parameters_TRiX.R")

#### Make electropherograms ####
fsa_directory <- "datas/"
fsa_files <- list.files(path = fsa_directory, pattern = "*.fsa$", full.names = TRUE)

read_fsa_to_df <- function(file) {
  seq_data <- unique(data.frame(seqinr::read.abif(file)$Data$DATA.1))
  seq_data$filename <- basename(file)
  seq_data$run_index <- rownames(seq_data)
  return(seq_data)
}

# Read all files into a list of data.frames
seq_list <- lapply(fsa_files, read_fsa_to_df)

# Combine into a single data.frame
combined_df <- do.call(rbind, seq_list)

# change colnames
colnames(combined_df) <- c("intensity", "filename", "run_index")

# separate file name
combined_df <- combined_df %>% tidyr::separate(filename, into = c("digest", "sample"))

# add index
combined_df$run_index <- as.numeric(combined_df$run_index)

# add gene symbol
combined_df$gene <- "LOL"
combined_df[which(combined_df$run_index > 2500 & combined_df$run_index < 3250),]$gene <- "ARHGAP6"
combined_df[which(combined_df$run_index > 3250 & combined_df$run_index < 3750),]$gene <- "CNKSR2"
combined_df[which(combined_df$run_index > 3750 & combined_df$run_index < 4500),]$gene <- "TCAC1"
combined_df[which(combined_df$run_index > 4500 & combined_df$run_index < 5500),]$gene <- "RP2"
combined_df[which(combined_df$run_index > 5500 & combined_df$run_index < 6000),]$gene <- "ZIC3"
 
# plot
pdf("01_plots/Figure_1D.pdf", width = 16, height = 8)
ggplot(combined_df[combined_df$run_index > 2500 & combined_df$run_index < 6001 & combined_df$sample == "320280",], 
       aes(x=run_index,
           y=intensity,
           col=digest))+
  geom_line()+
  facet_grid(factor(digest, levels = c("U", "D"))~.) + 
  theme_AL_box()
dev.off()


