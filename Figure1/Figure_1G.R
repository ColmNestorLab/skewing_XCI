library(tidyverse)
library(NanoMethViz)
library(data.table)
library(GenomicRanges)
library(tidyr)
library(dbscan)
library(scanstatistics)
library(stringr)
library(dplyr)
library(cowplot)
library(doParallel)
library(foreach)
library(Rsamtools)
library(GenomicAlignments)

#### Manual version below (pre-processing) ####

# Set BAM file path
bam_file <- "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/overlaps.tagged.phased.SNV.only.bam"

# # Define what to extract
param <- ScanBamParam(
  what = c("qname", "rname", "pos", "strand", "seq", "cigar")
)

# Read the BAM
reads_list <- scanBam(bam_file, param = param)
df <- as.data.frame(reads_list[[1]])

# Compute read length from sequence
df$read_length <- width(df$seq)

# Compute alignment end position from CIGAR
# Use GenomicAlignments for accurate end calculation
gal <- readGAlignments(bam_file)
end_pos <- end(gal)

# add stop for making GRanges object
df$end <- end_pos

# Convert to GRanges object
gr <- GRanges(seqnames = df$rname, ranges = IRanges(start = df$pos, end = df$end))

# Interval to check for overlap
query <- data.frame(chr = c("chrX","chrX","chrX","chrX","chrX","chrX"), start = c(46836331-1,137568340-1,137545369-1,11427816-1,21374595-1, 67545318-1), stop = c(46836331+1,137568340+1,137545369+1,11427816+1,21374595+1, 67545318+1))
query_gr <- GRanges(seqnames = query$chr, ranges = IRanges(start = query$start, end = query$stop))

# Find overlaps
overlap_idx <- findOverlaps(gr, query_gr)

# Extract overlapping rows
hitz <- df[queryHits(overlap_idx), ]

# Read list of read names
read_names <- unique(hitz$qname)

# Export
write.table(read_names, "read_names.txt", quote = F, row.names = F)

# Using samtools I extract all reads using the read_names.txt file:
# samtools view -@ 64 -N read_names.txt -b merged.allez.chrX.bam > hitz.merged.allez.chrX.bam

# I then get the actual sequences from the BAM files and convert it to a tsv file.
# samtools view hitz.merged.allez.chrX.bam | cut -f1,2,3,4,5,6,7,8,9,10 > hitz.merged.allez.chrX.incl.AR.tsv
# Read it in.
hitz_tsv <- fread("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/hitz.merged.allez.chrX.tsv")

grep_df <- data.frame(repeat_name = c("ARHGAP6", "CNKSR2", "RP2", "TCAC1", "ZIC3", "AR"),
                      repeat_start = c(11427816, 21374595, 46836331, 137545369, 137568340, 67545318),
                      patterns = toupper(c("GCGAAGGGggagggggaa", "GGAGCGGAGCGGCGGAGg", "cgagaccctgtGAAAGAA", "TACTAGAGCAGTACTTAT", "AGCCCCCTGGATCGATCG", "GCCAGTTTGCTGCTGCTG")),
                      repeat_type = c("GAG",
                                      "CAG",
                                      "AAAG",
                                      "TCAC",
                                      "ATCT",
                                      "CAG")
)

# Create a list to store results
split_dfs <- list()

# Loop through each pattern and filter matching rows
for (i in seq_along(grep_df$patterns)) {
  pattern <- grep_df$patterns[i]
  repeat_name <- grep_df$repeat_name[i]
  split_dfs[[repeat_name]] <- hitz[grep(pattern, hitz$seq, ignore.case = TRUE), ]
}

# Split dfs
ARHGAP6_df <- split_dfs$ARHGAP6
CNKSR2_df <- split_dfs$CNKSR2
RP2_df <- split_dfs$RP2
TCAC1_df <- split_dfs$TCAC1
ZIC3_df <- split_dfs$ZIC3
AR_df <- split_dfs$AR

# repeat types
ARHGAP6_pattern <- "GAG"
CNKSR2_pattern <- "CAG"
RP2_pattern <- "AAAG"
TCAC1_pattern <- "TCAC"
ZIC3_pattern <- "ATCT"
AR_pattern <- "CAG"

# Function to find the longest uninterrupted repeat and its sequence
longest_repeat_info <- function(seq, pattern) {
  matches <- str_extract_all(seq, paste0("(", pattern, ")+"))[[1]]  # Extract all contiguous occurrences
  if (length(matches) == 0) return(c(0, ""))  # Return 0 length and empty string if no match
  longest_match <- matches[which.max(nchar(matches))]  # Get the longest uninterrupted match
  return(c(nchar(longest_match), longest_match))  # Return length and sequence
}

total_repeat_bases <- function(seq, pattern) {
  matches <- str_extract_all(seq, paste0("(", pattern, ")+"))[[1]]  # Extract contiguous occurrences
  if (length(matches) == 0) return(0)  # Return 0 if no match
  longest_match <- matches[which.max(nchar(matches))]  # Get the longest match
  return(nchar(longest_match))  # Return the length of the longest match
}


# Apply function 
ARHGAP6_df <- ARHGAP6_df %>%
  rowwise() %>%
  mutate(info = list(longest_repeat_info(seq, "GAG"))) %>%
  mutate(total_repeat_bases = as.numeric(info[[1]]), 
         longest_repeat_sequence = info[[2]]) %>%
  select(-info)  # Remove intermediate list column

ARHGAP6_df <- ARHGAP6_df %>%
  mutate(total_repeat_bases = sapply(seq, total_repeat_bases, pattern = ARHGAP6_pattern))

ARHGAP6_df$repeats <- ARHGAP6_df$total_repeat_bases/3


# Apply function 
CNKSR2_df <- CNKSR2_df %>%
  rowwise() %>%
  mutate(info = list(longest_repeat_info(seq, "CAG"))) %>%
  mutate(total_repeat_bases = as.numeric(info[[1]]), 
         longest_repeat_sequence = info[[2]]) %>%
  select(-info)  # Remove intermediate list column

CNKSR2_df <- CNKSR2_df %>%
  mutate(total_repeat_bases = sapply(seq, total_repeat_bases, pattern = CNKSR2_pattern))

CNKSR2_df$repeats <- CNKSR2_df$total_repeat_bases/3


# Apply function 
RP2_df <- RP2_df %>%
  rowwise() %>%
  mutate(info = list(longest_repeat_info(seq, "AAAG"))) %>%
  mutate(total_repeat_bases = as.numeric(info[[1]]), 
         longest_repeat_sequence = info[[2]]) %>%
  select(-info)  # Remove intermediate list column

RP2_df <- RP2_df %>%
  mutate(total_repeat_bases = sapply(seq, total_repeat_bases, pattern = RP2_pattern))

RP2_df$repeats <- RP2_df$total_repeat_bases/4


# Apply function 
TCAC1_df <- TCAC1_df %>%
  rowwise() %>%
  mutate(info = list(longest_repeat_info(seq, "TCAC"))) %>%
  mutate(total_repeat_bases = as.numeric(info[[1]]), 
         longest_repeat_sequence = info[[2]]) %>%
  select(-info)  # Remove intermediate list column

TCAC1_df <- TCAC1_df %>%
  mutate(total_repeat_bases = sapply(seq, total_repeat_bases, pattern = TCAC1_pattern))

TCAC1_df$repeats <- TCAC1_df$total_repeat_bases/4


# Apply function 
ZIC3_df <- ZIC3_df %>%
  rowwise() %>%
  mutate(info = list(longest_repeat_info(seq, "ATCT"))) %>%
  mutate(total_repeat_bases = as.numeric(info[[1]]), 
         longest_repeat_sequence = info[[2]]) %>%
  select(-info)  # Remove intermediate list column

ZIC3_df <- ZIC3_df %>%
  mutate(total_repeat_bases = sapply(seq, total_repeat_bases, pattern = ZIC3_pattern))

ZIC3_df$repeats <- ZIC3_df$total_repeat_bases/4


# Apply function 
AR_df <- AR_df %>%
  rowwise() %>%
  mutate(info = list(longest_repeat_info(seq, "CAG"))) %>%
  mutate(total_repeat_bases = as.numeric(info[[1]]), 
         longest_repeat_sequence = info[[2]]) %>%
  select(-info)  # Remove intermediate list column

AR_df <- AR_df %>%
  mutate(total_repeat_bases = sapply(seq, total_repeat_bases, pattern = AR_pattern))

AR_df$repeats <- AR_df$total_repeat_bases/3

# Check repeat length
plot_grid(
  ggplot(ARHGAP6_df, aes(x=repeats)) + geom_histogram(binwidth = 1),
  ggplot(CNKSR2_df, aes(x=repeats)) + geom_histogram(binwidth = 1),
  ggplot(RP2_df, aes(x=repeats)) + geom_histogram(binwidth = 1),
  ggplot(TCAC1_df, aes(x=repeats)) + geom_histogram(binwidth = 1),
  ggplot(ZIC3_df, aes(x=repeats)) + geom_histogram(binwidth = 1),
  ggplot(AR_df, aes(x=repeats)) + geom_histogram(binwidth = 1)
)

ARHGAP6_df_14 <- ARHGAP6_df[ARHGAP6_df$repeats == 14,]
ARHGAP6_df_17_18_19 <- ARHGAP6_df[ARHGAP6_df$repeats %in% c(17,18,19),]

CNKSR2_df_11 <- CNKSR2_df[CNKSR2_df$repeats == 11,]
CNKSR2_df_12 <- CNKSR2_df[CNKSR2_df$repeats == 12,]

RP2_df_12 <- RP2_df[RP2_df$repeats == 12,]
RP2_df_14 <- RP2_df[RP2_df$repeats == 14,]

TCAC1_df_10 <- TCAC1_df[TCAC1_df$repeats == 10,]
TCAC1_df_12 <- TCAC1_df[TCAC1_df$repeats == 12,]

ZIC3_df_11 <- ZIC3_df[ZIC3_df$repeats == 11,]
ZIC3_df_12 <- ZIC3_df[ZIC3_df$repeats == 12,]

AR_df_20 <- AR_df[AR_df$repeats == 20,]
AR_df_23 <- AR_df[AR_df$repeats == 23,]


write.table(ARHGAP6_df_14$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ARHGAP6_df_14_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(ARHGAP6_df_17_18_19$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ARHGAP6_df_17_18_19_read_names.txt", quote = F, row.names = F, col.names = F)

write.table(CNKSR2_df_11$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/CNKSR2_df_11_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(CNKSR2_df_12$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/CNKSR2_df_12_read_names.txt", quote = F, row.names = F, col.names = F)

write.table(RP2_df_12$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/RP2_df_12_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(RP2_df_14$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/RP2_df_14_read_names.txt", quote = F, row.names = F, col.names = F)

write.table(TCAC1_df_10$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/TCAC1_df_10_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(TCAC1_df_12$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/TCAC1_df_12_read_names.txt", quote = F, row.names = F, col.names = F)

write.table(ZIC3_df_11$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ZIC3_df_11_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(ZIC3_df_12$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ZIC3_df_12_read_names.txt", quote = F, row.names = F, col.names = F)

write.table(AR_df_20$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/AR_df_20_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(AR_df_23$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/AR_df_23_read_names.txt", quote = F, row.names = F, col.names = F)


## after extracting the reads from the BAM files using an awk script (see below), I manually checked all loci that the reads are correctly phased across the repeats. Manually moving reads to the other haplogroup where the initial assignment above was incorrect.


# ARHGAP6
# delete two reads, as they cannot be assigned.
write.table(ARHGAP6_df_14[ARHGAP6_df_14$qname != "84991687-7de0-4ecf-9ac2-e82c66e34e14",]$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/ARHGAP6_df_14_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(ARHGAP6_df_17_18_19[ARHGAP6_df_17_18_19$qname != "352b5515-8efb-4796-83f8-151f2e74ce4b",], "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/ARHGAP6_df_17_18_19_read_names.txt", quote = F, row.names = F, col.names = F)

# ARHGAP6_df_14
#84991687-7de0-4ecf-9ac2-e82c66e34e14 del

# ARHGAP6_df_17_18_19
#352b5515-8efb-4796-83f8-151f2e74ce4b del

# CNKSR2
# move two reads.
write.table(CNKSR2_df_11[!CNKSR2_df_11$qname %in% c("83390a26-7f07-4e7c-8f01-a02f50b4b2dc", "ee9e42da-e984-4633-b7ae-889e08e40936"),]$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/CNKSR2_df_11_read_names.txt", quote = F, row.names = F, col.names = F)
write.table(c(CNKSR2_df_12$qname, "83390a26-7f07-4e7c-8f01-a02f50b4b2dc", "ee9e42da-e984-4633-b7ae-889e08e40936"), "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/CNKSR2_df_12_read_names.txt", quote = F, row.names = F, col.names = F)

# CNKSR2_df_11 ---> CNKSR2_df_12
#83390a26-7f07-4e7c-8f01-a02f50b4b2dc move
#ee9e42da-e984-4633-b7ae-889e08e40936 move

# RP2
# deleting a total of 3 reads that could not be assigned perfectly.
write.table(RP2_df_12[!RP2_df_12$qname %in% c("386ee39c-afb4-4590-9b48-890eec05f3f7",
                                           "aba8c605-04f2-4b9d-85ce-3f072fd24504",
                                           "27213315-9e75-4327-97dc-ed2f53fe1881"),]$qname, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/RP2_df_12_read_names.txt", quote = F, row.names = F, col.names = F)
#write.table(RP2_df_14[!RP2_df_14$V1 %in% c("ccde924b-dad4-4f55-a424-cb7d5486ced4"),]$V1, "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/RP2_df_14_read_names.txt", quote = F, row.names = F, col.names = F)

#RP2_12
#"386ee39c-afb4-4590-9b48-890eec05f3f7" del
#"aba8c605-04f2-4b9d-85ce-3f072fd24504" del
#"27213315-9e75-4327-97dc-ed2f53fe1881" del

#RP2_14
#unchanged

#TCAC1
# Removed no reads.

# ZIC3
# looks perfect, no changes.

# AR
# looks perfect, no changes.

# I then use the following script to produce one BAM file per repeat and phase.
#samtools view -@ 64 -N ${repeat}_df_${repeat_count}_read_names.txt -b hitz.merged.allez.chrX.bam > ${repeat}_df_${repeat_count}_hitz.merged.allez.chrX.bam

#### Plot split reads using split bam files #####

# get exons
exon_tibble <- get_exons_hg38()

# Read in modbams
mbr_ARHGAP6 <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ARHGAP6_df_14_hitz.merged.allez.chrX.bam",
      "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ARHGAP6_df_17_18_19_hitz.merged.allez.chrX.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
)

mbr_CNKSR2 <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/CNKSR2_df_11_hitz.merged.allez.chrX.bam",
      "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/CNKSR2_df_12_hitz.merged.allez.chrX.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
)

mbr_RP2 <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/RP2_df_12_hitz.merged.allez.chrX.bam",
      "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/RP2_df_14_hitz.merged.allez.chrX.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
)

mbr_TCAC1 <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/TCAC1_df_12_hitz.merged.allez.chrX.bam",
      "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/TCAC1_df_10_hitz.merged.allez.chrX.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
)

mbr_ZIC3 <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ZIC3_df_11_hitz.merged.allez.chrX.bam",
      "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/ZIC3_df_12_hitz.merged.allez.chrX.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
  
)

mbr_AR <- ModBamResult(
  methy = ModBamFiles(
    samples = c("haplotype1", "haplotype2"),
    c("/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/AR_df_20_hitz.merged.allez.chrX.bam",
      "/home/bjogy93/Desktop/TREX_proj/nanopore_data/analysis/split_bams/manual/AR_df_23_hitz.merged.allez.chrX.bam")
  ),
  samples = data.frame(
    sample = c("haplotype1", "haplotype2"),
    group = c(1, 2)
  ),
  exons = exon_tibble
  
)


ARHGAP6_PCR_prod_plot <- plot_region(mbr_ARHGAP6, chr = "chrX", start = 11427654, end = 11427890, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 11427816) + geom_vline(xintercept = c(11427770, 11427772), col = "red")
CNKSR2_PCR_prod_plot <- plot_region(mbr_CNKSR2, chr = "chrX", start = 21374408, end = 21374708, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 21374595) + geom_vline(xintercept = c(21374428,21374430, 21374576,21374578, 21374669, 21374671), col = "red")
RP2_PCR_prod_plot <- plot_region(mbr_RP2, chr = "chrX", start = 46836203, end = 46836692, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 46836331) + geom_vline(xintercept = c(46836650, 46836652), col = "red")
TCAC1_PCR_prod_plot <- plot_region(mbr_TCAC1, chr = "chrX", start = 137545307, end = 137545689, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 137545369) + geom_vline(xintercept = c(137545622, 137545624, 137545633, 137545635), col = "red")
ZIC3_PCR_prod_plot <- plot_region(mbr_ZIC3, chr = "chrX", start = 137568039, end = 137568606, avg_method = "median", gene_anno = F) + geom_vline(xintercept = 137568340) + geom_vline(xintercept = c(137568262,137568264, 137568507,137568509, 137568525,137568527), col = "red")

# Below calculate XCI skew using the CCGG sites.
ARHGAP6_meth <- query_methy(mbr_ARHGAP6, c("chrX"), c(11427770), c(11427772))
CNKSR2_meth <- query_methy(mbr_CNKSR2, c("chrX","chrX","chrX"), c(21374428,21374576,21374669), c(21374430,21374578,21374671))
RP2_meth <- query_methy(mbr_RP2, c("chrX"), c(46836650), c(46836652))
TCAC1_meth <- query_methy(mbr_TCAC1, c("chrX","chrX"), c(137545622,137545633), c(137545624,137545635))
ZIC3_meth <- query_methy(mbr_ZIC3, c("chrX","chrX","chrX"), c(137568262,137568507,137568525), c(137568264,137568509,137568527))


ARHGAP6_meth_mean <- ARHGAP6_meth
CNKSR2_meth_mean <- CNKSR2_meth %>% dplyr::group_by(read_name, sample) %>% dplyr::summarise(meanz = median(mod_prob)) %>% dplyr::ungroup()
RP2_meth_mean <- RP2_meth
TCAC1_meth_mean <- TCAC1_meth %>% dplyr::group_by(read_name, sample) %>% dplyr::summarise(meanz = median(mod_prob)) %>% dplyr::ungroup()
ZIC3_meth_mean <- ZIC3_meth %>% dplyr::group_by(read_name, sample) %>% dplyr::summarise(meanz = median(mod_prob)) %>% dplyr::ungroup()

ARHGAP6_meth_mean$meth <- "undetermined"
ARHGAP6_meth_mean[ARHGAP6_meth_mean$mod_prob < 0.1,]$meth <- "Xa"
ARHGAP6_meth_mean[ARHGAP6_meth_mean$mod_prob > 0.9,]$meth <- "Xi"

CNKSR2_meth_mean$meth <- "undetermined"
CNKSR2_meth_mean[CNKSR2_meth_mean$meanz < 0.1,]$meth <- "Xa"
CNKSR2_meth_mean[CNKSR2_meth_mean$meanz > 0.9,]$meth <- "Xi"

RP2_meth_mean$meth <- "undetermined"
RP2_meth_mean[RP2_meth_mean$mod_prob < 0.1,]$meth <- "Xa"
RP2_meth_mean[RP2_meth_mean$mod_prob > 0.9,]$meth <- "Xi"

TCAC1_meth_mean$meth <- "undetermined"
TCAC1_meth_mean[TCAC1_meth_mean$meanz < 0.1,]$meth <- "Xa"
TCAC1_meth_mean[TCAC1_meth_mean$meanz > 0.9,]$meth <- "Xi"

ZIC3_meth_mean$meth <- "undetermined"
ZIC3_meth_mean[ZIC3_meth_mean$meanz < 0.1,]$meth <- "Xa"
ZIC3_meth_mean[ZIC3_meth_mean$meanz > 0.9,]$meth <- "Xi"


ARHGAP6_meth_filt_count <- ARHGAP6_meth_mean %>%
  dplyr::count(sample, meth) %>%
  tidyr::complete(sample, meth, fill = list(n = 0))

CNKSR2_meth_filt_count <- CNKSR2_meth_mean %>%
  dplyr::count(sample, meth) %>%
  tidyr::complete(sample, meth, fill = list(n = 0))

RP2_meth_filt_count <- RP2_meth_mean %>%
  dplyr::count(sample, meth) %>%
  tidyr::complete(sample, meth, fill = list(n = 0))

TCAC1_meth_filt_count <- TCAC1_meth_mean %>%
  dplyr::count(sample, meth) %>%
  tidyr::complete(sample, meth, fill = list(n = 0))

ZIC3_meth_filt_count <- ZIC3_meth_mean %>%
  dplyr::count(sample, meth) %>%
  tidyr::complete(sample, meth, fill = list(n = 0))

ARHGAP6_wide <- pivot_wider(ARHGAP6_meth_filt_count[ARHGAP6_meth_filt_count$meth != "undetermined",], names_from = c(sample,meth), values_from = n)
CNKSR2_wide <- pivot_wider(CNKSR2_meth_filt_count[CNKSR2_meth_filt_count$meth != "undetermined",], names_from = c(sample,meth), values_from = n)
RP2_wide <- pivot_wider(RP2_meth_filt_count[RP2_meth_filt_count$meth != "undetermined",], names_from = c(sample,meth), values_from = n)
TCAC1_wide <- pivot_wider(TCAC1_meth_filt_count[TCAC1_meth_filt_count$meth != "undetermined",], names_from = c(sample,meth), values_from = n)
ZIC3_wide <- pivot_wider(ZIC3_meth_filt_count[ZIC3_meth_filt_count$meth != "undetermined",], names_from = c(sample,meth), values_from = n)
                                                                                                                                          
ARHGAP6_meth <- query_methy(mbr_ARHGAP6, c("chrX"), c(11427770), c(11427770+1))
CNKSR2_meth <- query_methy(mbr_CNKSR2, c("chrX","chrX","chrX"), c(21374429,21374576,21374669), c(21374429+1,21374576+1,21374669+1))
RP2_meth <- query_methy(mbr_RP2, c("chrX"), c(46836650), c(46836650+1))
TCAC1_meth <- query_methy(mbr_TCAC1, c("chrX","chrX"), c(137545622,137545633), c(137545622+1,137545633+1))
ZIC3_meth <- query_methy(mbr_ZIC3, c("chrX","chrX", "chrX"), c(137568262,137568507,137568525), c(137568262+1,137568507+1,137568525+1))
AR_meth <- query_methy(mbr_AR, c("chrX","chrX"), c(67545256,67545296), c(67545256+1,67545296+1))

ARHGAP6_wide_skew <- ARHGAP6_wide %>% mutate(H1_Xa_skew = (haplotype1_Xi + haplotype2_Xa) / (haplotype1_Xa + haplotype1_Xi + haplotype2_Xa + haplotype2_Xi))
CNKSR2_wide_skew <- CNKSR2_wide %>% mutate(H1_Xa_skew = (haplotype1_Xi + haplotype2_Xa) / (haplotype1_Xa + haplotype1_Xi + haplotype2_Xa + haplotype2_Xi))
RP2_wide_skew <- RP2_wide %>% mutate(H1_Xa_skew = (haplotype1_Xi + haplotype2_Xa) / (haplotype1_Xa + haplotype1_Xi + haplotype2_Xa + haplotype2_Xi))
TCAC1_wide_skew <- TCAC1_wide %>% mutate(H1_Xa_skew = (haplotype1_Xi + haplotype2_Xa) / (haplotype1_Xa + haplotype1_Xi + haplotype2_Xa + haplotype2_Xi))
ZIC3_wide_skew <- ZIC3_wide %>% mutate(H1_Xa_skew = (haplotype1_Xi + haplotype2_Xa) / (haplotype1_Xa + haplotype1_Xi + haplotype2_Xa + haplotype2_Xi))

round(base::mean(
  c(as.numeric(ARHGAP6_wide_skew$H1_Xa_skew),
    as.numeric(CNKSR2_wide_skew$H1_Xa_skew),
    as.numeric(RP2_wide_skew$H1_Xa_skew),
    as.numeric(TCAC1_wide_skew$H1_Xa_skew),
    as.numeric(ZIC3_wide_skew$H1_Xa_skew)
  )
),3
)

ARHGAP6_PCR_prod_plot_fin <- ARHGAP6_PCR_prod_plot + labs(caption= paste0(as.numeric(ARHGAP6_wide_skew$H1_Xa_skew),", ","CCGG is 8th from right"))
CNKSR2_PCR_prod_plot_fin <- CNKSR2_PCR_prod_plot + labs(caption= paste0(as.numeric(CNKSR2_wide_skew$H1_Xa_skew),", ","CCGG is 5rd and 15th from right, 3rd from left"))
RP2_PCR_prod_plot_fin <- RP2_PCR_prod_plot + labs(caption= paste0(as.numeric(RP2_wide_skew$H1_Xa_skew),", ","CCGG is 2nd from right"))
TCAC1_PCR_prod_plot_fin <- TCAC1_PCR_prod_plot + labs(caption= paste0(as.numeric(TCAC1_wide_skew$H1_Xa_skew),", ","CCGG is 2nd and 3rd from right"))
ZIC3_PCR_prod_plot_fin <- ZIC3_PCR_prod_plot + labs(caption= paste0(as.numeric(ZIC3_wide_skew$H1_Xa_skew),", ","CCGG is 7th, 9th and 17th from right"))

genes_PCR_prod_plot_manual <- 
  plot_grid(nrow=1,
            ARHGAP6_PCR_prod_plot_fin,
            CNKSR2_PCR_prod_plot_fin,
            RP2_PCR_prod_plot_fin,
            TCAC1_PCR_prod_plot_fin,
            ZIC3_PCR_prod_plot_fin)


pdf("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/01_plots/Figure_1G.pdf", height = 6, width = 24)
genes_PCR_prod_plot_manual
dev.off()



