library(data.table)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggbeeswarm)
library(cowplot)
library(ggalluvial)
library(ggrepel)
library(GenomicRanges)
library(rtracklayer)

source("plot_parameters_TRiX.R")

############ GTEx ############
suppl_table_1_from_eLife_paper <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/anno/Supplementary_Table1.tsv")#[!gene %in% c("PUDP", "ZRSR2", "PRKX", "MED14", "DDX3X", "PJA1", "KDM6A", "PIN4", "KDM5C", "CD99L2", "MSL3", "SEPTIN6", "ATRX", "TXLNG", "ARSD", "EIF1AX", "ZFX", "CA5B", "ENSG00000285756", "RBMX"),]

# Define directory containing TSV files
dir_path_ASE <- "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/GTEX_femmes/ASEReadCounter_output/"

# Get list of TSV files
files_ASE <- list.files(dir_path_ASE, pattern = "\\.table$", full.names = TRUE)

# Read files and add filename as a column
GTEx_femmes_ASE <- rbindlist(lapply(files_ASE, function(file) {
  fread(file)[, filename := basename(file)]
}))

# fix names
GTEx_femmes_ASE$filename <- gsub(GTEx_femmes_ASE$filename, pattern = "\\.table", replacement = "")

# add participant
GTEx_femmes_ASE[, participant := sub("^(([^-]*-){2}).*", "\\1", filename)]
GTEx_femmes_ASE[, participant := sub("-$", "", participant)]  # Remove trailing "-"


# read in WES
# Define directory containing TSV files
dir_path_WES <- "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/GTEX_femmes/WES_readcounts/"

# Get list of TSV files
files_WES <- list.files(dir_path_WES, pattern = "\\.tsv$", full.names = TRUE)

# Read files and add filename as a column
GTEx_femmes_WES <- rbindlist(lapply(files_WES, function(file) {
  fread(file)[, participant := basename(file)]
}))

# fix names
GTEx_femmes_WES$participant <- gsub(GTEx_femmes_WES$participant, pattern = "\\.WES.ReadCount.tsv", replacement = "")


# Merge
GTEx_femmes <- merge(GTEx_femmes_ASE, by.x = c("participant","contig","position", "refAllele", "altAllele"),
                     GTEx_femmes_WES, by.y = c("participant","V1","V2","V3","V4"),
                     all.x = T)


# add metadata
meta_rnaseq <- fread("anno/metadata_RNAseq_full_for_tissues.tsv")[,c("participant", "entity:sample_id", "tissue_id")]
meta_rnaseq$comment <- ""

# add old meta data
extra_anno <- fread("anno/GTEx_Analysis_v8_Annotations_SampleAttributesDS_with_shortnames.txt")[,c("SAMPID_SHORT", "SAMPID", "SMTSD", "SMPTHNTS")]
extra_anno <- extra_anno[!extra_anno$SAMPID %in% meta_rnaseq$`entity:sample_id`]

# Change colnames
colnames(extra_anno) <- c("participant","entity:sample_id", "tissue_id", "comment")

# read in partcipant annotations
participant_anno <- fread("anno/GTEx_participant.tsv")

# rbind metadata
meta_full <- rbind(meta_rnaseq, extra_anno)

# add age to metadata
meta_full <- merge(meta_full,
                   unique(participant_anno[,c("entity:participant_id", "age")]),
                   by.x = "participant",
                   by.y = "entity:participant_id",
                   all.x = T)

# rename tissues to be uniform
meta_full[meta_full$tissue_id == "Adipose - Subcutaneous",]$tissue_id <- "adipose-subc"
meta_full[meta_full$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-subc"
meta_full[meta_full$tissue_id == "Adipose_Visceral_Omentum",]$tissue_id <- "adipose-visc"
meta_full[meta_full$tissue_id == "Adipose - Visceral (Omentum)",]$tissue_id <- "adipose-visc"
meta_full[meta_full$tissue_id == "Adrenal Gland",]$tissue_id <- "adrenal gland"
meta_full[meta_full$tissue_id == "Artery_Aorta",]$tissue_id <- "artery-aort"
meta_full[meta_full$tissue_id == "Artery - Aorta",]$tissue_id <- "artery-aort"
meta_full[meta_full$tissue_id == "Artery_Coronary",]$tissue_id <- "artery-coro"
meta_full[meta_full$tissue_id == "Artery - Coronary",]$tissue_id <- "artery-coro"
meta_full[meta_full$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-tibi"
meta_full[meta_full$tissue_id == "Artery - Tibial",]$tissue_id <- "artery-tibi"
meta_full[meta_full$tissue_id == "Brain_Amygdala",]$tissue_id <- "brain-amyg"
meta_full[meta_full$tissue_id == "Brain - Amygdala",]$tissue_id <- "brain-amyg"
meta_full[meta_full$tissue_id == "Brain_Anterior_cingulate_cortex_BA24",]$tissue_id <- "brain-ante"
meta_full[meta_full$tissue_id == "Brain - Anterior cingulate cortex (BA24)",]$tissue_id <- "brain-ante"
meta_full[meta_full$tissue_id == "Brain - Caudate (basal ganglia)",]$tissue_id <- "brain-caud"
meta_full[meta_full$tissue_id == "Brain_Caudate_basal_ganglia",]$tissue_id <- "brain-caud"
meta_full[meta_full$tissue_id == "Brain_Cerebellar_Hemisphere",]$tissue_id <- "brain-cehe"
meta_full[meta_full$tissue_id == "Brain - Cerebellar Hemisphere",]$tissue_id <- "brain-cehe"
meta_full[meta_full$tissue_id == "Brain - Cerebellum",]$tissue_id <- "brain-cere"
meta_full[meta_full$tissue_id == "Brain_Cerebellum",]$tissue_id <- "brain-cere"
meta_full[meta_full$tissue_id == "Brain - Cortex",]$tissue_id <- "brain-cort"
meta_full[meta_full$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-cort"
meta_full[meta_full$tissue_id == "Brain_Frontal_Cortex_BA9",]$tissue_id <- "brain-frco"
meta_full[meta_full$tissue_id == "Brain - Frontal Cortex (BA9)",]$tissue_id <- "brain-frco"
meta_full[meta_full$tissue_id == "Brain - Hippocampus",]$tissue_id <- "brain-hipp"
meta_full[meta_full$tissue_id == "Brain_Hippocampus",]$tissue_id <- "brain-hipp"
meta_full[meta_full$tissue_id == "Brain - Hypothalamus",]$tissue_id <- "brain-hypo"
meta_full[meta_full$tissue_id == "Brain_Hypothalamus",]$tissue_id <- "brain-hypo"
meta_full[meta_full$tissue_id == "Brain - Nucleus accumbens (basal ganglia)",]$tissue_id <- "brain-nucl"
meta_full[meta_full$tissue_id == "Brain_Nucleus_accumbens_basal_ganglia",]$tissue_id <- "brain-nucl"
meta_full[meta_full$tissue_id == "Brain - Putamen (basal ganglia)",]$tissue_id <- "brain-puta"
meta_full[meta_full$tissue_id == "Brain_Putamen_basal_ganglia",]$tissue_id <- "brain-puta"
meta_full[meta_full$tissue_id == "Brain - Spinal cord (cervical c-1)",]$tissue_id <- "brain-spin"
meta_full[meta_full$tissue_id == "Brain_Spinal_cord_cervical_c-1",]$tissue_id <- "brain-spin"
meta_full[meta_full$tissue_id == "Brain - Substantia nigra",]$tissue_id <- "brain-subs"
meta_full[meta_full$tissue_id == "Brain_Substantia_nigra",]$tissue_id <- "brain-subs"
meta_full[meta_full$tissue_id == "Breast_Mammary_Tissue",]$tissue_id <- "breast"
meta_full[meta_full$tissue_id == "Breast - Mammary Tissue",]$tissue_id <- "breast"
meta_full[meta_full$tissue_id == "Minor_Salivary_Gland",]$tissue_id <- "salivary gland"
meta_full[meta_full$tissue_id == "Cervix_Ectocervix",]$tissue_id <- "cervix-ecto"
meta_full[meta_full$tissue_id == "Cervix_Endocervix",]$tissue_id <- "cervix-endo"
meta_full[meta_full$tissue_id == "Colon - Sigmoid",]$tissue_id <- "colon-sigm"
meta_full[meta_full$tissue_id == "Colon_Sigmoid",]$tissue_id <- "colon-sigm"
meta_full[meta_full$tissue_id == "Colon - Transverse",]$tissue_id <- "colon-tran"
meta_full[meta_full$tissue_id == "Colon_Transverse",]$tissue_id <- "colon-tran"
meta_full[meta_full$tissue_id == "Colon - Transverse",]$tissue_id <- "colon-tran"
meta_full[meta_full$tissue_id == "Esophagus_Gastroesophageal_Junction",]$tissue_id <- "esophagus-gaju"
meta_full[meta_full$tissue_id == "Esophagus - Gastroesophageal Junction",]$tissue_id <- "esophagus-gaju"
meta_full[meta_full$tissue_id == "Esophagus_Mucosa",]$tissue_id <- "esophagus-muco"
meta_full[meta_full$tissue_id == "Esophagus - Mucosa",]$tissue_id <- "esophagus-muco"
meta_full[meta_full$tissue_id == "Esophagus_Muscularis",]$tissue_id <- "esophagus-musc"
meta_full[meta_full$tissue_id == "Esophagus - Muscularis",]$tissue_id <- "esophagus-musc"
meta_full[meta_full$tissue_id == "Cells_Cultured_fibroblasts",]$tissue_id <- "fibroblasts"
meta_full[meta_full$tissue_id == "Heart_Atrial_Appendage",]$tissue_id <- "heart-atri"
meta_full[meta_full$tissue_id == "Heart - Atrial Appendage",]$tissue_id <- "heart-atri"
meta_full[meta_full$tissue_id == "Minor Salivary Gland",]$tissue_id <- "salivary gland"
meta_full[meta_full$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-vent"
meta_full[meta_full$tissue_id == "Heart - Left Ventricle",]$tissue_id <- "heart-vent"
meta_full[meta_full$tissue_id == "Kidney_Cortex",]$tissue_id <- "kidney-cort"
meta_full[meta_full$tissue_id == "Kidney - Cortex",]$tissue_id <- "kidney-cort"
meta_full[meta_full$tissue_id == "Kidney_Medulla",]$tissue_id <- "kidney-medu"
meta_full[meta_full$tissue_id == "Cells_EBV-transformed_lymphocytes",]$tissue_id <- "lymphocytes"
meta_full[meta_full$tissue_id == "Cells - EBV-transformed lymphocytes",]$tissue_id <- "lymphocytes"
meta_full[meta_full$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
meta_full[meta_full$tissue_id == "Muscle - Skeletal",]$tissue_id <- "muscle"
meta_full[meta_full$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
meta_full[meta_full$tissue_id == "Nerve - Tibial",]$tissue_id <- "nerve"
meta_full[meta_full$tissue_id == "Skin - Not Sun Exposed (Suprapubic)",]$tissue_id <- "skin-supr"
meta_full[meta_full$tissue_id == "Skin_Not_Sun_Exposed_Suprapubic",]$tissue_id <- "skin-supr"
meta_full[meta_full$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-lleg"
meta_full[meta_full$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-lleg"
meta_full[meta_full$tissue_id == "Skin - Sun Exposed (Lower leg)",]$tissue_id <- "skin-lleg"
meta_full[meta_full$tissue_id == "Small_Intestine_Terminal_Ileum",]$tissue_id <- "small intestine"
meta_full[meta_full$tissue_id == "Small Intestine - Terminal Ileum",]$tissue_id <- "small intestine"
meta_full[meta_full$tissue_id == "Adrenal_Gland",]$tissue_id <- "adrenal gland"
meta_full[meta_full$tissue_id == "Bladder",]$tissue_id <- "bladder"
meta_full[meta_full$tissue_id == "Liver",]$tissue_id <- "liver"
meta_full[meta_full$tissue_id == "Lung",]$tissue_id <- "lung"
meta_full[meta_full$tissue_id == "Pancreas",]$tissue_id <- "pancreas"
meta_full[meta_full$tissue_id == "Pituitary",]$tissue_id <- "pituitary"
meta_full[meta_full$tissue_id == "Spleen",]$tissue_id <- "spleen"
meta_full[meta_full$tissue_id == "Stomach",]$tissue_id <- "stomach"
meta_full[meta_full$tissue_id == "Thyroid",]$tissue_id <- "thyroid"
meta_full[meta_full$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"
meta_full[meta_full$tissue_id == "Whole Blood",]$tissue_id <- "whole blood"
meta_full[meta_full$tissue_id == "Uterus",]$tissue_id <- "uterus"
meta_full[meta_full$tissue_id == "Vagina",]$tissue_id <- "vagina"
meta_full[meta_full$tissue_id == "Ovary",]$tissue_id <- "ovary"
meta_full[meta_full$tissue_id == "Fallopian_Tube",]$tissue_id <- "fallopian tube"

GTEx_femmes_meta <- merge(GTEx_femmes, meta_full, by.x=c("filename", "participant")  , by.y= c("entity:sample_id", "participant"), all.x =T)

GTEx_femmes_meta <- GTEx_femmes_meta[!GTEx_femmes_meta$tissue_id %in% c("Cells - Cultured fibroblasts", "lymphocytes", "fibroblasts"),]

count_samples_per_sample <- unique(GTEx_femmes_meta[,c("participant", "tissue_id")]) %>% dplyr::group_by(participant) %>% dplyr::count()

# require at least 3 tissues
count_samples_per_sample <- count_samples_per_sample[count_samples_per_sample$n > 2,]

GTEx_femmes_meta <- GTEx_femmes_meta[GTEx_femmes_meta$participant %in% count_samples_per_sample$participant,]

length(unique(GTEx_femmes_meta$participant))
nrow(unique(GTEx_femmes_meta[,c("participant", "tissue_id")]))


# Ensure femmes is a data.table
setDT(GTEx_femmes_meta)

# Load GENCODE GTF
gtf <- import("anno/gencode.v47.basic.annotation.gtf")

# Filter for exons and UTRs
gtf_filtered <- gtf[gtf$type %in% c("exon", "UTR")]

# Convert GTF to GRanges
gtf_gr <- GRanges(
  seqnames = seqnames(gtf_filtered),
  ranges = IRanges(start(gtf_filtered), end(gtf_filtered)),
  gene_name = mcols(gtf_filtered)$gene_name
)

# Convert GTEx_femmes_meta to GRanges (assuming 'position' is a single base position)
GTEx_femmes_gr <- GRanges(
  seqnames = GTEx_femmes_meta$contig,
  ranges = IRanges(start = GTEx_femmes_meta$position, end = GTEx_femmes_meta$position)
)

# Find overlaps
hits <- findOverlaps(GTEx_femmes_gr, gtf_gr)

# Add gene names to GTEx_femmes_meta data.table
GTEx_femmes_meta[queryHits(hits), gene_name := mcols(gtf_gr)$gene_name[subjectHits(hits)]]

# View annotated GTEx_femmes_meta data
head(GTEx_femmes_meta)


# remove unused columns
GTEx_femmes_meta <- GTEx_femmes_meta[,c("filename", "participant", "contig", "position", "refAllele", "altAllele", "variantID", "refCount", "altCount", "totalCount","V5",    "V6",    "V7",    "tissue_id", "comment",   "age", "gene_name")]

colnames(GTEx_femmes_meta) <- c("filename", "participant", "contig", "position", "refAllele", "altAllele", "variantID", "refCount", "altCount", "totalCount","refVariantDepth",    "altVariantDepth",    "totalVariantDepth",    "tissue_id", "comment",   "age", "gene_name")


# filter
GTEx_femmes_meta$read_count_filter <- ifelse((GTEx_femmes_meta$totalVariantDepth >= 20 & GTEx_femmes_meta$totalCount >= 7), yes = TRUE, no = FALSE)

GTEx_femmes_meta$alt_ref_filter <- ifelse((GTEx_femmes_meta$altVariantDepth >= 10 & GTEx_femmes_meta$refVariantDepth >= 10), yes = TRUE, no = FALSE) 

GTEx_femmes_meta$minor_allele_filter <- ifelse((GTEx_femmes_meta$refVariantDepth/GTEx_femmes_meta$totalVariantDepth >= 0.1 & GTEx_femmes_meta$altVariantDepth/GTEx_femmes_meta$totalVariantDepth >= 0.1) | (GTEx_femmes_meta$refCount/GTEx_femmes_meta$totalCount >= 0.1 & GTEx_femmes_meta$altCount/GTEx_femmes_meta$totalCount >= 0.1) , yes = TRUE, no = FALSE)

GTEx_femmes_meta$pass_all_filters <- ifelse(
  GTEx_femmes_meta$read_count_filter == TRUE & 
    GTEx_femmes_meta$alt_ref_filter == TRUE, 
  yes = TRUE, no = FALSE)


# calculate AE
GTEx_femmes_meta$effectSize <- abs(0.5 - GTEx_femmes_meta$refCount / (GTEx_femmes_meta$refCount + GTEx_femmes_meta$altCount))

# add KEEP tag to the hetSNP with the highest read count per tissue and individual
GTEx_femmes_meta[pass_all_filters == TRUE, keep := order(totalCount,decreasing = T) == 1, by=c("tissue_id","gene_name", "participant")]

# Remove genes for screening for sXCI females.
GTEx_femmes_not_filtered <- GTEx_femmes_meta[GTEx_femmes_meta$contig == "chrX",]
GTEx_femmes_filt <- GTEx_femmes_meta[GTEx_femmes_meta$gene %in% suppl_table_1_from_eLife_paper$gene & keep == T,]

# Check no. of samples and participants
length(unique(GTEx_femmes_filt$participant))
nrow(unique(GTEx_femmes_filt[,c("participant", "tissue_id")]))

# get tissues with 10 or more samples.
tissues_to_keep <- unique(GTEx_femmes_filt[,c("tissue_id", "participant")]) %>% dplyr::group_by(tissue_id) %>% dplyr::count() %>% subset(n>=10)

# Exclude tissues with very few samples
GTEx_femmes_filt <- GTEx_femmes_filt[GTEx_femmes_filt$tissue_id %in% tissues_to_keep$tissue_id,]

# Rename tissues
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "adipose-subc",]$tissue_id <- "adipose-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "adipose-visc",]$tissue_id <- "adipose-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "adrenal gland",]$tissue_id <- "adrenal gland"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "artery-aort",]$tissue_id <- "artery-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "artery-coro",]$tissue_id <- "artery-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "artery-tibi",]$tissue_id <- "artery-3"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-amyg",]$tissue_id <- "brain-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-ante",]$tissue_id <- "brain-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-caud",]$tissue_id <- "brain-3"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-cehe",]$tissue_id <- "brain-4"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-cere",]$tissue_id <- "brain-5"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-cort",]$tissue_id <- "brain-6"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-frco",]$tissue_id <- "brain-7"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-hipp",]$tissue_id <- "brain-8"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-hypo",]$tissue_id <- "brain-9"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-nucl",]$tissue_id <- "brain-10"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-puta",]$tissue_id <- "brain-11"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-spin",]$tissue_id <- "brain-12"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "brain-subs",]$tissue_id <- "brain-13"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "breast",]$tissue_id <- "breast"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "colon-sigm",]$tissue_id <- "colon-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "colon-tran",]$tissue_id <- "colon-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "esophagus-gaju",]$tissue_id <- "esophagus-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "esophagus-muco",]$tissue_id <- "esophagus-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "esophagus-musc",]$tissue_id <- "esophagus-3"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "heart-atri",]$tissue_id <- "heart-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "heart-vent",]$tissue_id <- "heart-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "kidney-cort",]$tissue_id <- "kidney"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "kidney cortex",]$tissue_id <- "kidney"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "liver",]$tissue_id <- "liver"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "lung",]$tissue_id <- "lung"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "muscle",]$tissue_id <- "muscle"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "nerve",]$tissue_id <- "nerve"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "ovary",]$tissue_id <- "ovary"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "pancreas",]$tissue_id <- "pancreas"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "pituitary",]$tissue_id <- "pituitary"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "salivary gland",]$tissue_id <- "salivary gland"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "skin-supr",]$tissue_id <- "skin-1"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "skin-lleg",]$tissue_id <- "skin-2"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "small intestine",]$tissue_id <- "small intestine"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "spleen",]$tissue_id <- "spleen"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "stomach",]$tissue_id <- "stomach"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "thyroid",]$tissue_id <- "thyroid"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "uterus",]$tissue_id <- "uterus"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "vagina",]$tissue_id <- "vagina"
GTEx_femmes_filt[GTEx_femmes_filt$tissue_id == "whole blood",]$tissue_id <- "whole blood"

# calculate metrics for plotting.
suppl_table_2_metrics <- GTEx_femmes_filt %>% dplyr::group_by(participant, tissue_id) %>% rstatix::get_summary_stats(effectSize, type = "common")
suppl_table_2_metrics$tissue_id <- factor(suppl_table_2_metrics$tissue_id)

suppl_table_2_metrics_order <- suppl_table_2_metrics %>% dplyr::group_by(participant) %>% rstatix::get_summary_stats(median, type = "common")

suppl_table_2_metrics_order2 <- suppl_table_2_metrics %>% dplyr::group_by(tissue_id) %>% rstatix::get_summary_stats(median, type = "common")

##### Figure 3b, show degree of skewing for TRiXi and GTEx whole blood #####
##### read in data #####
df_miniTRiX <- fread("data/minitrix_results.tsv")
df_TRiX <- fread("data/trix_results.tsv")

df_miniTRiX$method <- "mini"
df_TRiX$method <- "full"

melt_df_miniTRiX <- melt(data.table(df_miniTRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing"))
melt_df_miniTRiX$value <- as.numeric(melt_df_miniTRiX$value)

melt_df_TRiX <- melt(data.table(df_TRiX), id.vars = c("sample", "method"), measure.vars = c("mean_skewing","ARHGAP6_skewing", "RP2_skewing", "CNKSR2_skewing", "TCAC1_skewing", "ZIC3_skewing"))
melt_df_TRiX$value <- as.numeric(melt_df_TRiX$value)

# rbind
TRiXi_data <- rbind(melt_df_miniTRiX,
                    melt_df_TRiX)

# calculate descriptive statistics.
stats_TRiXi_data <- TRiXi_data %>% dplyr::group_by(sample, method) %>% rstatix::get_summary_stats(value, type = "common")

# calculate descriptive statistics of descriptive statistics, to plot lines mean, and 1sd and 2sd.
stats_stats_TRiXi_data <- stats_TRiXi_data %>% rstatix::get_summary_stats(mean, type = "common")

# Classify XCI pattern of individuals.
stats_TRiXi_data$'XCI pattern' <- "mosaic"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+1*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "skewed"
stats_TRiXi_data[stats_TRiXi_data$mean > (stats_stats_TRiXi_data$mean+2*stats_stats_TRiXi_data$sd),]$'XCI pattern' <- "extremely skewed"
stats_TRiXi_data$'XCI pattern' <- factor(stats_TRiXi_data$'XCI pattern', levels = c("extremely skewed", "skewed", "mosaic"))
stats_TRiXi_data_short <- stats_TRiXi_data[,c("sample","variable","mean")]

stats_suppl_table_2_metrics_order_short <- suppl_table_2_metrics[suppl_table_2_metrics$tissue_id == "whole blood",c("participant","variable","median")]

colnames(stats_TRiXi_data_short) <- c("sample", "variable",  "skewing")
colnames(stats_suppl_table_2_metrics_order_short) <- c("sample", "variable",  "skewing")

# not scaled df
smash_dat <- rbind(stats_TRiXi_data_short, stats_suppl_table_2_metrics_order_short)

smash_dat$variable <- gsub(smash_dat$variable, pattern = "value", replacement = "TRiXi - whole blood")
smash_dat$variable <- gsub(smash_dat$variable, pattern = "effectSize", replacement = "GTEx - whole blood")

scale_range <- function(x, new_min = 50, new_max = 100, old_min = 0, old_max = 0.5) {
  (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
}

test_scale <- smash_dat[smash_dat$variable %in% c("GTEx - whole blood"),]
test_scale$skewing <- scale_range(test_scale$skewing)

fin_scaled <- rbind(test_scale, smash_dat[smash_dat$variable %in% c("TRiXi - whole blood"),])
fin_scaled$skewed <- ifelse(fin_scaled$skewing >= 80, yes = "skewed", no = "not_skewed")
fin_scaled_count <- fin_scaled %>% dplyr::group_by(variable, skewed) %>% dplyr::count()

fin_scaled_count <- fin_scaled_count %>%
  group_by(variable) %>%
  mutate(pct = n / sum(n))

# plot density
pdf("01_plots/Figure_3C.pdf", height = 5)
plot_grid(rel_widths = c(0.9,0,1), ncol = 2,
          ggdensity(fin_scaled, x="skewing", col = "variable", add = "median") + 
            geom_vline(xintercept = 80, lty=2) + 
            annotate(geom="text", x=c(85,85), y=c(0.015,0.005), 
                     label= round(fin_scaled_count[fin_scaled_count$skewed == "skewed",]$pct, 3)*100)
)
dev.off()
