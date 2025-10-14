# Load necessary library
library(ggplot2)
library(data.table)
library(biomaRt)
library(karyoploteR)
library(tidyverse)
library(cowplot)
library(Biostrings)

# good to have.
source("plot_parameters_TRiX.R")

# repeat data
trip_data <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/trip_repeats_sequences_to_grep.tsv", 
                   col.names = c("contig", "repeat_start", "repeat_start_20_before_100_after", "repeat_type", "sequence_directly_after_repeat"))

quad_data <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_sequences_to_grep.tsv", 
                   col.names = c("contig", "repeat_start", "repeat_start_20_before_100_after", "repeat_type", "sequence_directly_after_repeat"))

# read in repeat meta
trip_repeats_meta <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/trip_repeats_meta.tsv", 
                   col.names = c("contig","repeat_start","repeat_end","repeat_start_minus_300","repeat_end_plus_300"))

quad_repeats_meta <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_meta.tsv", 
                           col.names = c("contig","repeat_start","repeat_end","repeat_start_minus_300","repeat_end_plus_300"))

# calculate no. of repeats.
trip_repeats_meta$repeat_count <- round((trip_repeats_meta$repeat_end-trip_repeats_meta$repeat_start)/3, 0)
quad_repeats_meta$repeat_count <- round((quad_repeats_meta$repeat_end-quad_repeats_meta$repeat_start)/4, 0)

# add repeat info to the data frame
trip_data <- merge(trip_data, trip_repeats_meta, by = c("contig", "repeat_start"))
quad_data <- merge(quad_data, quad_repeats_meta, by = c("contig", "repeat_start"))

# Merge repeat type and count for nicer plotting further down.
trip_data$repeat_type_and_count <- paste0(trip_data$repeat_type,"(", trip_data$repeat_count, ")")
quad_data$repeat_type_and_count <- paste0(quad_data$repeat_type,"(", quad_data$repeat_count, ")")

# get ensembl database hg38
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# get gene info for chrX
genes <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
               filters = 'chromosome_name', values = 'X', mart = ensembl)


# Function to find the closest gene for a given position
find_closest_gene <- function(pos) {
  # Calculate the distance to each gene (start or end position)
  genes <- genes %>%
    mutate(distance = pmin(abs(start_position - pos), abs(end_position - pos))) %>%
    arrange(distance)
  
  # If the position is within a gene, return that gene; otherwise, return the closest one
  within_gene <- genes %>% filter(start_position <= pos & end_position >= pos)
  
  if (nrow(within_gene) > 0) {
    return(within_gene$hgnc_symbol[1])
  } else {
    return(genes$hgnc_symbol[1])
  }
}


# Apply the function to each position in the repeat data.
try(trip_data <- trip_data %>%
  rowwise() %>%
  mutate(closest_gene = find_closest_gene(repeat_start)), silent = T)

# workaround ensembl being down.
#trip_data$closest_gene <- "whaaat"

# add gene name. The repeats are not necessarily in the actual gene, here I just add the closest gene.
trip_data[trip_data$repeat_start %in% 37737518,]$closest_gene <- "XK"
trip_data[trip_data$repeat_start %in% 38562575,]$closest_gene <- "TSPAN7"
trip_data[trip_data$repeat_start %in% 48631673,]$closest_gene <- "WDR13"
trip_data[trip_data$repeat_start %in% 69082485,]$closest_gene <- "PJA1"
trip_data[trip_data$repeat_start %in% 69083812,]$closest_gene <- "PJA1_2"
trip_data[trip_data$repeat_start %in% 123767927,]$closest_gene <- "THOC"

# selected
trip_data[trip_data$repeat_start %in% 67545317,]$closest_gene <- "AR"
trip_data[trip_data$repeat_start %in% 11427816,]$closest_gene <- "ARHGAP6"
	

try(quad_data <- quad_data %>%
  rowwise() %>%
  mutate(closest_gene = find_closest_gene(repeat_start)), silent = T)

# workaround ensembl being down.
#quad_data$closest_gene <- "whaaat"

# add gene name. The repeats are not necessarily in the actual gene, here I just add the closest gene.
quad_data[quad_data$repeat_start %in% 833522,]$closest_gene <- "SHOX"
quad_data[quad_data$repeat_start %in% 1109474,]$closest_gene <- "CRLF2"
quad_data[quad_data$repeat_start %in% 68903924,]$closest_gene <- "EFNB1"
quad_data[quad_data$repeat_start %in% 127008981,]$closest_gene <- "PRR32"
quad_data[quad_data$repeat_start %in% 1658563,]$closest_gene <- "ASMT"
quad_data[quad_data$repeat_start %in% 4551340,]$closest_gene <- "FAM239B"
quad_data[quad_data$repeat_start %in% 39613271,]$closest_gene <- "MIR3937"
quad_data[quad_data$repeat_start %in% 48244717,]$closest_gene <- "SSX1"
quad_data[quad_data$repeat_start %in% 88417303,]$closest_gene <- "CPXCR1"
quad_data[quad_data$repeat_start %in% 137554666,]$closest_gene <- "LINC02931_2"

# selected
quad_data[quad_data$repeat_start %in% 137545369,]$closest_gene <- "TCAC1"
quad_data[quad_data$repeat_start %in% 137568340,]$closest_gene <- "ZIC3"
quad_data[quad_data$repeat_start %in% 46836337,]$closest_gene <- "RP2"

# Read in methylation data. Remove positions with a read depth lower than 4.
df_meth <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/50_meth_summarized.sorted.samples.merged.bed", col.names = c("contig", "CG_start", "CG_end", "read_depth", "methylation_value"))[read_depth > 4]
df_unmeth <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/not_50_meth_summarized.sorted.samples.merged.bed", col.names = c("contig", "CG_start", "CG_end", "read_depth", "methylation_value"))[read_depth > 4]

# Merge the methylation data.
df_smash_meth <- rbind(df_meth, df_unmeth)[methylation_value > 0]

# read in CG and CCGG positions. Keep only 50% methylated CCGG positions.
df_CG_pos <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/CGpos.txt", col.names = c("CG_position"))
df_CCGG <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/CCGGpos.txt", col.names = c("CCGG_position"))
df_CCGG_pos <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/CCGGpos.txt", col.names = c("CCGG_position"))[CCGG_position %in% df_meth$CG_start]


trip_gr2 <- makeGRangesFromDataFrame(trip_data[trip_data$repeat_start %in% c("11427816", "21374595", "46836337", "137545369", "137568340"),], keep.extra.columns = T, seqnames.field = "contig", start.field="repeat_start_minus_300", end.field = "repeat_end_plus_300")
quad_gr2 <- makeGRangesFromDataFrame(quad_data[quad_data$repeat_start %in% c("11427816", "21374595", "46836337", "137545369", "137568340"),], keep.extra.columns = T, seqnames.field = "contig", start.field="repeat_start_minus_300", end.field = "repeat_end_plus_300")

# get methylation, add tag whether or not its a CCGG site
df_smash_meth$CCGG <- ifelse(df_smash_meth$CG_start %in% df_CCGG$CCGG_position, yes = "CCGG", no = "CpG")

# make GRanges object for plotting.
ccgg_methylation_gr <- makeGRangesFromDataFrame(df_smash_meth[df_smash_meth$CCGG == "CCGG",], keep.extra.columns = T, seqnames.field = "contig", start.field="CG_start", end.field = "CG_end" )

cg_methylation_gr <- makeGRangesFromDataFrame(df_smash_meth[df_smash_meth$CCGG == "CpG",], keep.extra.columns = T, seqnames.field = "contig", start.field="CG_start", end.field = "CG_end" )


# selected repeats, show primer binding sites.
selected_repeats <- c("AR","ARHGAP6", "CNKSR2", "RP2", "TCAC1", "ZIC3")

smash_data <- rbind(trip_data[trip_data$closest_gene %in% selected_repeats,],
                    quad_data[quad_data$closest_gene %in% selected_repeats,])

primer_summary <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/primers_summary_TRiXi_corrected.tsv")

smash_data_with_primer_summary <- merge(smash_data, primer_summary[,c("gene", "PCR_product_start", "PCR_product_end", "PCR_product_length", "PCR_product_sequence")], by.x = "closest_gene", by.y="gene")


# add GC% to the plot.
# got sequences from IGV and added to primer_summary.tsv
# calc GC content
gc_content_ARHGAP6 <- letterFrequency(DNAString(smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "ARHGAP6",]$PCR_product_sequence), letters = c("G", "C"), as.prob = TRUE)
gc_content_ARHGAP6 <- sum(gc_content_ARHGAP6) * 100

gc_content_CNKSR2 <- letterFrequency(DNAString(smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "CNKSR2",]$PCR_product_sequence), letters = c("G", "C"), as.prob = TRUE)
gc_content_CNKSR2 <- sum(gc_content_CNKSR2) * 100

gc_content_RP2 <- letterFrequency(DNAString(smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "RP2",]$PCR_product_sequence), letters = c("G", "C"), as.prob = TRUE)
gc_content_RP2 <- sum(gc_content_RP2) * 100

gc_content_TCAC1 <- letterFrequency(DNAString(smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "TCAC1",]$PCR_product_sequence), letters = c("G", "C"), as.prob = TRUE)
gc_content_TCAC1 <- sum(gc_content_TCAC1) * 100

gc_content_ZIC3 <- letterFrequency(DNAString(smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "ZIC3",]$PCR_product_sequence), letters = c("G", "C"), as.prob = TRUE)
gc_content_ZIC3 <- sum(gc_content_ZIC3) * 100

# add GC % column
smash_data_with_primer_summary$GC_procent <- "NO"
smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "ARHGAP6",]$GC_procent <- gc_content_ARHGAP6
smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "CNKSR2",]$GC_procent <- gc_content_CNKSR2
smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "RP2",]$GC_procent <- gc_content_RP2
smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "TCAC1",]$GC_procent <- gc_content_TCAC1
smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "ZIC3",]$GC_procent <- gc_content_ZIC3

smash_data_with_primer_summary$GC_position <- smash_data_with_primer_summary$PCR_product_start + (smash_data_with_primer_summary$PCR_product_length/2)

smash_data_with_primer_summary$GC_procent <- round(as.numeric(smash_data_with_primer_summary$GC_procent), digits = 1)

# adjust ZIC3 GC label position
smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "ZIC3",]$GC_position <- smash_data_with_primer_summary[smash_data_with_primer_summary$closest_gene == "ZIC3",]$GC_position -100


trip_starts <- trip_data[trip_data$repeat_start %in% c("11427816", "21374595", "46836337", "137545369", "137568340"),]$repeat_start
pdf(file = "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/01_plots/Fig_1A_trinucleotide_repeats.pdf", height = 7, width = 7)
for (i in trip_starts){
  temp <- trip_data[trip_data$repeat_start %in% i,]
  temp2 <- smash_data_with_primer_summary[smash_data_with_primer_summary$repeat_start %in% i,]
  zoom.region <- toGRanges(data.frame(temp$contig, temp$repeat_start_minus_300, temp$repeat_end_plus_300))
  
  rectangel_area_forward_primer <- toGRanges(data.frame(temp2$contig,  temp2$PCR_product_start, temp2$PCR_product_start+20))
  rectangel_area_reverse_primer <- toGRanges(data.frame(temp2$contig,  temp2$PCR_product_end-20, temp2$PCR_product_end))
  
  chr <- "chrX"
  start <- temp$repeat_start_minus_300
  end <- temp$repeat_end_plus_300
  
  ccgg_meth.ymax <- ceiling(max(abs(range(ccgg_methylation_gr$methylation_value))))
  ccgg_meth.ymin <- 0
  
  cpg_meth.ymax <- ceiling(max(abs(range(cg_methylation_gr$methylation_value))))
  cpg_meth.ymin <- 0
  
  kp <- plotKaryotype(chromosomes = "chrX", zoom = zoom.region, plot.type = 1)
  at <- autotrack(current.track = 1, total.tracks = 1)
  kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "lightgrey")
  kpPoints(kp, data=ccgg_methylation_gr, y=ccgg_methylation_gr$methylation_value, r0=at$r0, r1=at$r1, col = "red")
  kpAxis(kp, ymax=ccgg_meth.ymax, ymin=ccgg_meth.ymin, r0=at$r0, r1=at$r1)
  kpPoints(kp, data=cg_methylation_gr, y=cg_methylation_gr$methylation_value, r0=at$r0, r1=at$r1, col = "black")
  title(main=paste0(temp$closest_gene))
  kpPlotGenes(kp, data=genes, r0=0.2, r1=0.8)
  kpPlotMarkers(kp, chr=temp$contig, x=temp$repeat_start, labels=temp$repeat_type_and_count, y = 0.1)
  kpPlotMarkers(kp, chr=temp$contig, x=temp$repeat_end, labels=temp$repeat_type_and_count, y = 0.1)
  kpAddBaseNumbers(kp, 
                   tick.dist=100,           
                   minor.tick.dist=50,      
                   digits=10,               
                   cex=0.8,
                   add.units = TRUE)        
  kpRect(kp, chr = seqnames(rectangel_area_forward_primer), x0 = start(rectangel_area_forward_primer), x1 = end(rectangel_area_forward_primer),
         y0 = 0, y1 = 0.1, col = "red", border = "black", lwd = 1)
  kpPlotMarkers(kp, chr=seqnames(rectangel_area_forward_primer), x=(start(rectangel_area_forward_primer)+end(rectangel_area_forward_primer))/2, labels="forward_primer", y = 0.1)
  
  kpRect(kp, chr = seqnames(rectangel_area_reverse_primer), x1 = start(rectangel_area_reverse_primer), x0 = end(rectangel_area_reverse_primer),
         y0 = 0, y1 = 0.1, col = "red", border = "black", lwd = 1)
  kpPlotMarkers(kp, chr=seqnames(rectangel_area_reverse_primer), x=(start(rectangel_area_reverse_primer)+end(rectangel_area_reverse_primer))/2, labels="reverse_primer", y = 0.1)
  
  kpPlotMarkers(kp, chr=temp2$contig, x=temp2$GC_position, labels=paste0(temp2$GC_procent, " CG (%)"), y = 0.5)
  
}

dev.off()



quad_starts <- quad_data[quad_data$repeat_start %in% c("11427816", "21374595", "46836337", "137545369", "137568340"),]$repeat_start
pdf(file = "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/01_plots/Fig_1A_tetranucleotide_repeats.pdf", height = 7, width = 7)
for (i in quad_starts){
  temp <- quad_data[quad_data$repeat_start %in% i,]
  temp2 <- smash_data_with_primer_summary[smash_data_with_primer_summary$repeat_start %in% i,]
  zoom.region <- toGRanges(data.frame(temp$contig, temp$repeat_start_minus_300, temp$repeat_end_plus_300))
  
  rectangel_area_forward_primer <- toGRanges(data.frame(temp2$contig,  temp2$PCR_product_start, temp2$PCR_product_start+20))
  rectangel_area_reverse_primer <- toGRanges(data.frame(temp2$contig,  temp2$PCR_product_end-20, temp2$PCR_product_end))
  
  chr <- "chrX"
  start <- temp$repeat_start_minus_300
  end <- temp$repeat_end_plus_300
  
  ccgg_meth.ymax <- ceiling(max(abs(range(ccgg_methylation_gr$methylation_value))))
  ccgg_meth.ymin <- 0
  
  cpg_meth.ymax <- ceiling(max(abs(range(cg_methylation_gr$methylation_value))))
  cpg_meth.ymin <- 0
  
  kp <- plotKaryotype(chromosomes = "chrX", zoom = zoom.region, plot.type = 1)
  at <- autotrack(current.track = 1, total.tracks = 1)
  kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "lightgrey")
  kpPoints(kp, data=ccgg_methylation_gr, y=ccgg_methylation_gr$methylation_value, r0=at$r0, r1=at$r1, col = "red")
  kpAxis(kp, ymax=ccgg_meth.ymax, ymin=ccgg_meth.ymin, r0=at$r0, r1=at$r1)
  kpPoints(kp, data=cg_methylation_gr, y=cg_methylation_gr$methylation_value, r0=at$r0, r1=at$r1, col = "black")
  title(main=paste0(temp$closest_gene))
  kpPlotGenes(kp, data=genes, r0=0.2, r1=0.8)
  kpPlotMarkers(kp, chr=temp$contig, x=temp$repeat_start, labels=temp$repeat_type_and_count, y = 0.1)
  kpPlotMarkers(kp, chr=temp$contig, x=temp$repeat_end, labels=temp$repeat_type_and_count, y = 0.1)
  kpAddBaseNumbers(kp, 
                   tick.dist=100,           
                   minor.tick.dist=50,      
                   digits=10,               
                   cex=0.8,
                   add.units = TRUE)        
  kpRect(kp, chr = seqnames(rectangel_area_forward_primer), x0 = start(rectangel_area_forward_primer), x1 = end(rectangel_area_forward_primer), y0 = 0, y1 = 0.1, col = "red", border = "black", lwd = 1)
  kpPlotMarkers(kp, chr=seqnames(rectangel_area_forward_primer), x=(start(rectangel_area_forward_primer)+end(rectangel_area_forward_primer))/2, labels="forward_primer", y = 0.1)
  
  kpRect(kp, chr = seqnames(rectangel_area_reverse_primer), x1 = start(rectangel_area_reverse_primer), x0 = end(rectangel_area_reverse_primer), y0 = 0, y1 = 0.1, col = "red", border = "black", lwd = 1)
  kpPlotMarkers(kp, chr=seqnames(rectangel_area_reverse_primer), x=(start(rectangel_area_reverse_primer)+end(rectangel_area_reverse_primer))/2, labels="reverse_primer", y = 0.1)
  
  kpPlotMarkers(kp, chr=temp2$contig, x=temp2$GC_position, labels=paste0(temp2$GC_procent, " CG (%)"), y = 0.5)
  
}
dev.off()