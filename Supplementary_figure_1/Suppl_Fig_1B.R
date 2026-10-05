library(rtracklayer)
library(GenomicRanges)

# directory containing .bw files
bw_dir <- "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/Kaplan_2023_female_blood_cell_samples"

# list all bigwig files
bw_files <- list.files(bw_dir, pattern="\\.bigwig$", full.names=TRUE)

region <- GRanges(
  seqnames = "chrX",
  ranges = IRanges(
    start = c(11427654, 21374408, 46836203, 137545307, 137568039),
    end   = c(11427890, 21374708, 46836692, 137545689, 137568606)
  )
)

bw_list <- lapply(bw_files, function(f) {
  import(f, format="BigWig", which=region)
})
names(bw_list) <- basename(bw_files)


df_list <- lapply(names(bw_list), function(nm) {
  df <- as.data.frame(bw_list[[nm]])
  df$file <- nm
  df
})

df <- do.call(rbind, df_list)

df$file <- gsub(df$file, pattern = ".hg38.bigwig", replacement = "")

df$loc <- "LOL"
df[df$start >= 11427654 & df$start <= 11427890,]$loc <- "TRiXi1"
df[df$start >= 21374408 & df$start <= 21374708,]$loc <- "TRiXi2"
df[df$start >= 46836203 & df$start <= 46836692,]$loc <- "TRiXi3"
df[df$start >= 137545307 & df$start <= 137545689,]$loc <- "TRiXi4"
df[df$start >= 137568039 & df$start <= 137568606,]$loc <- "TRiXi5"
df <- df[df$loc != "LOL",]


pattern <- "^GSM[0-9]+_.*[-_]Z[0-9A-Z]+$"

has_full <- grepl(pattern, df$file)

df$GSM <- NA_character_
df$ZID <- NA_character_
df$tissue <- NA_character_

# GSM
df$GSM[has_full] <- sub("^(GSM[0-9]+)_.*$", "\\1", df$file[has_full])

# ZID (works for both -Z and _Z)
df$ZID[has_full] <- sub("^.*[-_](Z[0-9A-Z]+)$", "\\1", df$file[has_full])

# tissue (core fix here)
df$tissue[has_full] <- sub("^GSM[0-9]+_(.*)[-_]Z[0-9A-Z]+$", "\\1", df$file[has_full])

# normalize
df$tissue <- gsub("-", "_", df$tissue)

unique(df$tissue)

# fix WBC
df_WBC <- df[is.na(df$tissue),]

df_WBC <- df_WBC %>% tidyr::separate(file, into = c("GSM"), sep = "-", remove = F) 
df_WBC <- df_WBC %>% tidyr::separate(GSM, into = c("GSM"), sep = "_") 
df_WBC$tissue <- "WBC"

df_WBC$ZID <- df_WBC$GSM

df <- df[!is.na(df$tissue),]

df <- rbind(df, df_WBC)

library(tidyverse)

file <- "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/Kaplan_2023_female_blood_cell_samples/GSE186458_series_matrix.txt"

lines <- readLines(file)

meta_lines <- lines[grepl("^!Sample", lines)]

meta_split <- strsplit(meta_lines, "\t")
max_len <- max(lengths(meta_split))

meta_df <- do.call(rbind, lapply(meta_split, function(x) {
  length(x) <- max_len
  x
})) %>% as.data.frame(stringsAsFactors = FALSE)

# ---- CLEAN + TRANSPOSE (NO dplyr magic) ----
mat <- as.matrix(meta_df)

rownames(mat) <- make.unique(mat[,1])   # use first column as rownames
mat <- mat[,-1]                         # remove first column

meta_t <- as.data.frame(t(mat), stringsAsFactors = FALSE)

# ---- ADD SAMPLE ID ----
meta_t$sample <- meta_t$`!Sample_geo_accession`

# ---- CLEAN CHARACTERISTICS ----
char_cols <- grep("^!Sample_characteristics", colnames(meta_t), value = TRUE)

meta_char <- meta_t %>%
  select(sample, all_of(char_cols)) %>%
  pivot_longer(-sample, values_to = "val") %>%
  filter(!is.na(val), val != "") %>%
  separate(val, into = c("key", "value"), sep = ": ", fill = "right", extra = "merge") %>%
  mutate(key = stringr::str_trim(key),
         value = stringr::str_trim(value)) %>%
  distinct() %>%
  pivot_wider(names_from = key, values_from = value)

# ---- FINAL TABLE ----
meta_clean <- meta_t %>%
  select(-all_of(char_cols)) %>%
  left_join(meta_char, by = "sample")

# ---- CLEAN NAMES ----
colnames(meta_clean) <- gsub("^!Sample_", "", colnames(meta_clean))

meta_clean[] <- lapply(meta_clean, function(x) {
  if (is.character(x)) gsub('"', '', x) else x
})

meta_clean_short <- unique(meta_clean[,c("sample", "title", "geo_accession", "source_name_ch1")])

meta_clean_short$ZID <- sub("^.*[-_](Z[0-9A-Z]+)$", "\\1", meta_clean_short$title)

# make sure we only have female samples
meta_WBC <- fread("data/Kaplan_2023_female_blood_cell_samples/female_WBC_meta.txt")
meta_other_tissues <- fread("data/Kaplan_2023_female_blood_cell_samples/samples.tsv")

meta_WBC <- meta_WBC %>% tidyr::separate(V1, into = "ZID", sep = "_")
meta_other_tissues <- meta_other_tissues %>% tidyr::separate(alias, into = "ZID", sep = "_")
meta_other_tissues_femme <- meta_other_tissues[meta_other_tissues$biological_sex == "female",]

samples_keep <- unique(c(meta_WBC$ZID, meta_other_tissues_femme$ZID))

# merge df and meta
df_meta <- merge(df[df$ZID %in% samples_keep,], 
                 meta_clean_short, 
                 by = "ZID",
                 all.x = T)

# add CG or CCGG tag
df_meta$type <- "CG"
df_meta[df_meta$start > 11427770 & df_meta$start < 11427772, ]$type  <- "CCGG"
df_meta[(df_meta$start > 21374428 & df_meta$start < 21374430) | (df_meta$start > 21374576 & df_meta$start < 21374578) | (df_meta$start > 21374669 & df_meta$start < 21374671) , ]$type  <- "CCGG"
df_meta[df_meta$start > 46836650 & df_meta$start < 46836652, ]$type  <- "CCGG"
df_meta[(df_meta$start > 137545622 & df_meta$start < 137545624) | (df_meta$start > 137545633 & df_meta$start < 137545635) , ]$type  <- "CCGG"
df_meta[(df_meta$start > 137568262 & df_meta$start < 137568264) | (df_meta$start > 137568507 & df_meta$start < 137568509) | (df_meta$start > 137568525 & df_meta$start < 137568527) , ]$type  <- "CCGG"

# get summaries
df_meta_stats <- df_meta[df_meta$type == "CCGG",] %>% dplyr::group_by(start, loc, type, tissue) %>% rstatix::get_summary_stats(score, type = "common")

# blood cells
blood_cells <- c("WBC","Blood_B", "Blood_B_Mem", "Blood_Granulocytes", "Blood_Monocytes", "Blood_NK", "Blood_T_CD3", "Blood_T_CD4", "Blood_T_CD8","Blood_T_CenMem_CD4","Blood_T_Eff_CD8","Blood_T_EffMem_CD4", "Blood_T_EffMem_CD8","Blood_T_Naive_CD4","Blood_T_Naive_CD8","Bone_marrow_Erythrocyte_progenitors" )

# it the sample blood cells or not?
df_meta_stats$blood_cells_or_no <- ifelse(df_meta_stats$tissue %in% blood_cells, yes = "blood cells", no = "no")

# keep only CCGG sites
df_meta_stats <- df_meta_stats[df_meta_stats$type == "CCGG",]

# save table
write.table(df_meta_stats, "source_data/Suppl_Fig_1B_source_data.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# read in source data
Suppl_Fig_1B_source_data <- fread("source_data/Suppl_Fig_1B_source_data.tsv")

pdf("01_plots/Suppl_fig_1B.pdf", width = 12, height = 6)

ggplot(Suppl_Fig_1B_source_data,
       aes(x=factor(start),
           y=mean,
           col=type,
           ymin=mean-sd,
           ymax=mean+sd)) + 
  stat_summary(fun.data = mean_sdl)+
  geom_quasirandom(alpha=.25)+
  geom_hline(yintercept = 0.5, lty = 2, col = "black") + 
  scale_color_manual(values = c("CG"="black", "CCGG"="red")) +
  ggplot2::facet_grid(blood_cells_or_no~loc, space = "free", scale = "free")

dev.off()

