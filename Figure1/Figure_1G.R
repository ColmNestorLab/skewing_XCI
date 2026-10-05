library(tidyverse)
library(NanoMethViz)
library(data.table)
library(cowplot)

#### Plot split reads using split bam files #####

# get exons
exon_tibble <- get_exons_hg38()

# Read in modbams
mbr <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("data/320280_BAM/hp1.bam",
      "data/320280_BAM/hp2.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
)

ARHGAP6_PCR_prod_plot <- plot_region(mbr, chr = "chrX", start = 11427654, end = 11427890, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 11427816) + geom_vline(xintercept = c(11427770, 11427772), col = "red")
CNKSR2_PCR_prod_plot <- plot_region(mbr, chr = "chrX", start = 21374408, end = 21374708, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 21374595) + geom_vline(xintercept = c(21374428,21374430, 21374576,21374578, 21374669, 21374671), col = "red")
RP2_PCR_prod_plot <- plot_region(mbr, chr = "chrX", start = 46836203, end = 46836692, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 46836331) + geom_vline(xintercept = c(46836650, 46836652), col = "red")
TCAC1_PCR_prod_plot <- plot_region(mbr, chr = "chrX", start = 137545307, end = 137545689, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 137545369) + geom_vline(xintercept = c(137545622, 137545624, 137545633, 137545635), col = "red")
ZIC3_PCR_prod_plot <- plot_region(mbr, chr = "chrX", start = 137568039, end = 137568606, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 137568340) + geom_vline(xintercept = c(137568262,137568264, 137568507,137568509, 137568525,137568527), col = "red")

genes_PCR_prod_plot <- 
  plot_grid(nrow=1,
            ARHGAP6_PCR_prod_plot,
            CNKSR2_PCR_prod_plot,
            RP2_PCR_prod_plot,
            TCAC1_PCR_prod_plot,
            ZIC3_PCR_prod_plot)

pdf("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/01_plots/Figure_1G.pdf", height = 6, width = 24)
genes_PCR_prod_plot
dev.off()
