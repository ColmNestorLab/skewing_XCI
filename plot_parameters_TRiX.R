# author of ggplot code below: Antonio Lentini(AL).
# modifies ggplots to look the same across plots.
require(ggplot2)
theme_AL_simple <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_simple_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_box_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 45), ... )
}
theme_AL_simple_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 45), ... )
}

chrom_colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF", "#F39B7FFF")

chrom_colors_gg <- c("1" = "#E64B35FF","2" = "#4DBBD5FF","3" ="#00A087FF","4" ="#3C5488FF","5" ="#F39B7FFF","6" ="#8491B4FF","7" ="#91D1C2FF","8" ="#7E6148FF","9" ="#B09C85FF","10" ="#E64B35FF","11" ="#4DBBD5FF","12" ="#00A087FF","13" ="#3C5488FF","14" ="#F39B7FFF","15" ="#8491B4FF","16" ="#91D1C2FF","17" ="#7E6148FF","18" ="#B09C85FF","19" ="#E64B35FF","20" ="#4DBBD5FF","21" ="#00A087FF","22" ="#3C5488FF", "X" ="#F39B7FFF")


# cell colors
colors = c("mosaic"="#1f77b4","skewed"="#ff7f0e","extremely skewed"="#2ca02c")

# factor level
XCI_pattern_factor = c("mosaic","skewed","extremely skewed")


# repeat colors
repeat_colors_TRiX <- c("ARHGAP6" = "#009E73", "CNKSR2" = "#56B4E9", "RP2" = "#CC79A7", "TCAC1" = "#000000", "ZIC3" = "#E69F00")
repeat_colors_miniTRiX <- c("ARHGAP6" = "#009E73", "CNKSR2" = "#56B4E9", "RP2" = "#CC79A7", "TCAC1" = "#000000", "ZIC3" = "#E69F00")

# Chromosome order
chrom_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX")


color_valzz <- c("brain-1" = "#0072B5FF",
                 "brain-2" ="#0072B5FF" ,
                 "brain-3" ="#0072B5FF",
                 "brain-4" = "#0072B5FF", 
                 "brain-5" = "#0072B5FF" ,
                 "brain-6" = "#0072B5FF" ,
                 "brain-7" = "#0072B5FF" ,
                 "brain-8" = "#0072B5FF" ,
                 "brain-9" = "#0072B5FF",
                 "brain-10" = "#0072B5FF",
                 "brain-11" = "#0072B5FF",
                 "brain-12" = "#0072B5FF",
                 "brain-13" ="#0072B5FF",
                 "esophagus-1" = "#4DBBD5FF" ,
                 "esophagus-2" ="#4DBBD5FF",
                 "esophagus-3" = "#4DBBD5FF",
                 "skin-1" ="#00A087FF",
                 "skin-2" = "#00A087FF",
                 "heart-1" = "#3C5488FF",
                 "heart-2" = "#3C5488FF",
                 "adipose-1" = "#F39B7FFF",
                 "adipose-2" ="#F39B7FFF",
                 "artery-1" = "#91D1C2FF",
                 "artery-2" = "#91D1C2FF" ,
                 "artery-3" = "#91D1C2FF",
                 "adrenal gland" = "#8491B4FF",
                 "kidney" = "#7E6148FF",
                 "bladder" = "#DC0000FF",
                 "nerve" = "#B09C85FF",
                 "salivary gland" = "#BC3C29FF",
                 "thyroid" = "#E64B35FF",
                 "pituitary" = "#E18727FF",
                 "pancreas" = "#20854EFF",
                 "spleen" = "#7876B1FF",
                 "liver" = "#6F99ADFF",
                 "lung" = "#FFDC91FF",
                 "breast" ="#EE4C97FF",
                 "muscle" = "#374E55FF",
                 "small intestine" = "#DF8F44FF",
                 "colon-1" = "#DF8F44FF",
                 "colon-2" = "#DF8F44FF",
                 "stomach" ="#DF8F44FF",
                 "lymphocytes" = "grey",
                 "fibroblasts" = "#B24745FF",
                 "whole blood" = "#79AF97FF",
                 "ovary" = "red",
                 "uterus" = "#6F99ADFF",
                 "vagina" = "#EE4C97FF")