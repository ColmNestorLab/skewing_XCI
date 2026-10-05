library(NanoMethViz)
library(tidyverse)
library(doParallel)
library(foreach)
library(data.table)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)


# ASE tables dir
file_dir_ASE_screen <- dir("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/ASEReadCounter", full.names = T) 

# read in and add file name as column
df_ase_screen <- do.call(rbind, lapply(file_dir_ASE_screen, function(x) cbind(read.csv(x, header = T, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))

# remove dir from file name column
df_ase_screen$name <- gsub(df_ase_screen$name, pattern = "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/ASEReadCounter/", replacement = "")

# WES data dir
file_dir_WES_screen <- dir("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/readcounts/chrX", full.names = T) 

# read in WES data and add file name as column
df_wes_screen <- do.call(rbind, lapply(file_dir_WES_screen, function(x) cbind(read.csv(x, header = F, sep = "\t"), name=strsplit(x,'_ReadCount')[[1]][1])))

# remove dir from file name column
df_wes_screen$name <- gsub(df_wes_screen$name, pattern = "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/readcounts/chrX/chrX.", replacement = "")


# make sure all numerical columns are numerical.
df_ase_screen$position <- as.numeric(df_ase_screen$position)
df_wes_screen$V2 <- as.numeric(df_wes_screen$V2)

# merge WES and ASE dfs
df_merged <- merge(df_ase_screen, df_wes_screen, by.x = c("name", "position", "contig"), by.y = c("name", "V2", "V1"), all.x = T)

# change column order and rename columns
colnames(df_merged) <- c("participant","position", "contig","variantID","refAllele","altAllele",
                         "refCount", "altCount","totalCount","lowMAPQDepth", "lowBaseQDepth" ,"rawDepth" ,"otherBases" ,"improperPairs","refAllele_DNA", "altAllele_DNA",
                         "refCount_WES","altCount_WES")

# make into data table
df_merged <- as.data.table(df_merged)


# add gene names.
chrom_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX")

## load gene annotations.
gtf_raw <- import("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/anno/gencode.v49.annotation.gtf.gz")

# exclude all but standard contigs.
df_merged <- df_merged[df_merged$contig %in% chrom_order]

## add gene annotations to merged WES/ASE df.
gr.ase <- with(df_merged, GRanges(seqnames=contig, IRanges(position,width = 1)))
ol.ase <- findOverlaps(gr.ase, gtf_raw[gtf_raw$type == c("exon", "UTR")])
nm.ase <- tapply(gtf_raw[gtf_raw$type == c("exon", "UTR")]$gene_name[subjectHits(ol.ase)],queryHits(ol.ase),function(x) paste0(unique(x),collapse=";") )
df_merged[, gene := NA]
df_merged$gene[as.numeric(names(nm.ase))] <- nm.ase

# make sure all numerical columns are numerical.
df_merged$altCount <- as.numeric(df_merged$altCount)
df_merged$refCount <- as.numeric(df_merged$refCount)
df_merged$totalCount <- as.numeric(df_merged$totalCount)

df_merged$altCount_WES <- as.numeric(df_merged$altCount_WES)
df_merged$refCount_WES <- as.numeric(df_merged$refCount_WES)
df_merged$totalCount_WES <- as.numeric(df_merged$altCount_WES+df_merged$refCount_WES)


# Add a raw read count filter
df_merged$read_count_filter <- ifelse((df_merged$totalCount_WES >= 10 & df_merged$totalCount >= 8), yes = TRUE, no = FALSE)

# Add raw alt/ref read count filter
df_merged$alt_ref_filter <- ifelse((df_merged$altCount_WES > 9 | df_merged$refCount_WES > 9), yes = TRUE, no = FALSE) 

# Add a $minor_allele_filter tag, where both and minor allele reads need to be equal to or more than 10 % of total reads, either in the WES or RNA-seq data. Essentially we make sure we observe both alleles somewhere in the data.
df_merged$minor_allele_filter <- ifelse((df_merged$refCount_WES/df_merged$totalCount_WES >= 0.1 & df_merged$altCount_WES/df_merged$totalCount_WES >= 0.1) | (df_merged$refCount/df_merged$totalCount >= 0.1 & df_merged$altCount/df_merged$totalCount >= 0.1) , yes = TRUE, no = FALSE)

# Add pass all filter tag
df_merged$pass_all_filters <- ifelse(df_merged$minor_allele_filter == TRUE & df_merged$read_count_filter == TRUE & df_merged$alt_ref_filter == TRUE, yes = TRUE, no = FALSE)

# remove not passing
df_merged <- df_merged[df_merged$pass_all_filters == T]

# calculate effectSize
df_merged$effectSize <- abs(0.5 - df_merged$refCount / (df_merged$refCount + df_merged$altCount))

# store full data
df_merged_full <- df_merged

# add keep filter to only keep hSNP per gene with highest coverage
df_merged_full[, keep := order(totalCount,decreasing = T) == 1, by=c("participant","gene")]

# export
write.table(df_merged_full, "Supplementary_tables/ASE_skewing_table_raw.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

ASE_df <- fread("Supplementary_tables/ASE_skewing_table_raw.tsv")

# prepare ASE df
# remove known escapees. Cover few genes so sensitive to outliers
gois <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/anno/gylemo_Elife_table_D.txt")
gois <- gois[gois$new_category == "inactive_across_tissues",]
gois <- unique(gois$gene)

# filter ASE data frame, keep only inactive genes that passes filtering. Remove any variant with more than one context.
ASE_df_filt <- ASE_df[ASE_df$gene %in% gois & ASE_df$pass_all_filters == T & ASE_df$lowBaseQDepth < 10 & ASE_df$otherBases == 0 & nchar(ASE_df$refAllele_DNA) == 1 & nchar(ASE_df$altAllele_DNA) == 1,]
ASE_df_filt <- ASE_df_filt[nchar(ASE_df_filt$refAllele_DNA) <= 1, ]
ASE_df_filt <- ASE_df_filt[nchar(ASE_df_filt$altAllele_DNA) <= 1, ]

# add scale (making it equivalent to values from 50 to 100 instead of 0 to 0.5).
ASE_df_filt$scaled <- 50 + 100 * ASE_df_filt$effectSize

# Export supplementary table S5
write.table(ASE_df_filt, "Supplementary_tables/supplementary_table5.tsv", quote = F, row.names = F, sep = "\t", col.names = T)

