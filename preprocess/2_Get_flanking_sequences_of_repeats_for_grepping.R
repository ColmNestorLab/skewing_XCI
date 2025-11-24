# Read in libraries.
library(data.table)
library(Biostrings)
library(BSgenome)
library(dplyr)
library(tidyr)
library(purrr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Before this script, the bash script 'bash_preprocess_repeats.sh' is run. It outputs all the files required for this script.

# Read in methylation data. Remove positions with a read depth lower than 4.
df_meth <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/50_meth_summarized.sorted.samples.merged.bed", col.names = c("contig", "CG_start", "CG_end", "read_depth", "methylation_value"))[read_depth > 4]
df_unmeth <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/not_50_meth_summarized.sorted.samples.merged.bed", col.names = c("contig", "CG_start", "CG_end", "read_depth", "methylation_value"))[read_depth > 4]

# Merge the methylation data.
df_smash_meth <- rbind(df_meth, df_unmeth)

# read in CG and CCGG positions. Keep only 50% methylated CCGG positions.
df_CG_pos <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/CGpos.txt", col.names = c("CG_position"))
df_CCGG_pos <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/CCGGpos.txt", col.names = c("CCGG_position"))[CCGG_position %in% df_meth$CG_start]

# Read in triplet, quadruplet and quintuple repeat positions
triplet_repeat_position <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/triplet_repeats_position.txt", col.names = c("repeat_start", "repeat_end"))
quad_repeat_position <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_position.txt", col.names = c("repeat_start", "repeat_end"))
quin_repeat_position <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quintuple_repeats_position.txt", col.names = c("repeat_start", "repeat_end"))

# Add/remove 300 to end/start of the CpG
triplet_repeat_position$start_pos_minus_300 <- triplet_repeat_position$repeat_start - 300
triplet_repeat_position$end_pos_plus_300 <- triplet_repeat_position$repeat_end + 300
quad_repeat_position$start_pos_minus_300 <- quad_repeat_position$repeat_start - 300
quad_repeat_position$end_pos_plus_300 <- quad_repeat_position$repeat_end + 300
quin_repeat_position$start_pos_minus_300 <- quin_repeat_position$repeat_start - 300
quin_repeat_position$end_pos_plus_300 <- quin_repeat_position$repeat_end + 300


check_overlap <- function(df1, df2, position_col, range_start_col, range_end_col) {
  # Initialize a column to store results in df1
  df1$overlap <- FALSE
  df1$range_start <- NA
  
  # Loop through each position in df1
  for (i in 1:nrow(df1)) {
    position <- df1[[position_col]][i]
    
    # Find the index where the position is within the ranges in df2
    overlap_indices <- which(position >= df2[[range_start_col]] & position <= df2[[range_end_col]])
    
    # If there is an overlap, update df1
    if (length(overlap_indices) > 0) {
      df1$overlap[i] <- TRUE
      df1$range_start[i] <- df2[[range_start_col]][overlap_indices[1]]
    }
  }
  
  return(df1)
}

# See if a CCGG position is found within 300bp of the repeat.
result_trip <- check_overlap(df_CCGG_pos, triplet_repeat_position, "CCGG_position", "start_pos_minus_300", "end_pos_plus_300")
result_quad <- check_overlap(df_CCGG_pos, quad_repeat_position, "CCGG_position", "start_pos_minus_300", "end_pos_plus_300")
result_quin <- check_overlap(df_CCGG_pos, quin_repeat_position, "CCGG_position", "start_pos_minus_300", "end_pos_plus_300")

# Keep only overlaps
hits_result_trip <- result_trip[result_trip$overlap == T]
hits_result_quad <- result_quad[result_quad$overlap == T]
hits_result_quin <- result_quin[result_quin$overlap == T]

# merge with repeat position
hits_result_trip <- merge(hits_result_trip, triplet_repeat_position, by.x = c("range_start"), by.y = c("start_pos_minus_300"))
hits_result_quad <- merge(hits_result_quad, quad_repeat_position, by.x = c("range_start"), by.y = c("start_pos_minus_300"))
hits_result_quin <- merge(hits_result_quin, quin_repeat_position, by.x = c("range_start"), by.y = c("start_pos_minus_300"))

hits_result_trip$contig <- "chrX"
hits_result_quad$contig <- "chrX"
hits_result_quin$contig <- "chrX"

hits_result_trip <- hits_result_trip[,c("contig", "repeat_start", "repeat_end", "CCGG_position", "range_start", "end_pos_plus_300")]
hits_result_quad <- hits_result_quad[,c("contig", "repeat_start", "repeat_end", "CCGG_position", "range_start", "end_pos_plus_300")]
hits_result_quin <- hits_result_quin[,c("contig", "repeat_start", "repeat_end", "CCGG_position", "range_start", "end_pos_plus_300")]

colnames(hits_result_trip) <- c("contig", "repeat_start", "repeat_end", "CCGG_position", "repeat_start_minus_300", "repeat_end_plus_300")
colnames(hits_result_quad) <- c("contig", "repeat_start", "repeat_end", "CCGG_position", "repeat_start_minus_300", "repeat_end_plus_300")
colnames(hits_result_quin) <- c("contig", "repeat_start", "repeat_end", "CCGG_position", "repeat_start_minus_300", "repeat_end_plus_300")

hits_result_trip_unique_positions <- unique(hits_result_trip[,-c("CCGG_position")])
hits_result_quad_unique_positions <- unique(hits_result_quad[,-c("CCGG_position")])
hits_result_quin_unique_positions <- unique(hits_result_quin[,-c("CCGG_position")])

# keep CCGG position
#hits_result_trip_unique_positions <- unique(hits_result_trip)
#hits_result_quad_unique_positions <- unique(hits_result_quad)
#hits_result_quin_unique_positions <- unique(hits_result_quin)

# Function to check for overlap and return multiple matches
check_range_overlap <- function(df1, df2, start1_col, stop1_col, start2_col, stop2_col, extra_cols) {
  df1 <- df1 %>%
    rowwise() %>%
    mutate(
      overlaps = list(which(df2[[start2_col]] <= get(stop1_col) & df2[[stop2_col]] >= get(start1_col)))
    ) %>%
    ungroup() %>%
    filter(lengths(overlaps) > 0) %>%  # Keep only rows with overlaps
    unnest_longer(overlaps) %>%
    mutate(
      overlap_start = df2[[start2_col]][overlaps],
      overlap_read_depth = df2[[extra_cols[1]]][overlaps],
      overlap_methylation_value = df2[[extra_cols[2]]][overlaps]
    ) %>%
    dplyr::select(-overlaps)
  
  return(df1)
}

# Apply the function to check overlap and add the overlap_start, read_depth, and methylation_value columns to hits_result_trip/quad/quin_unique_positions
hits_result_trip_unique_positions_overlap_GCs <- check_range_overlap(hits_result_trip_unique_positions, df_smash_meth, "repeat_start_minus_300", "repeat_end_plus_300", "CG_start", "CG_end", c("read_depth", "methylation_value"))

hits_result_quad_unique_positions_overlap_GCs <- check_range_overlap(hits_result_quad_unique_positions, df_smash_meth, "repeat_start_minus_300", "repeat_end_plus_300", "CG_start", "CG_end", c("read_depth", "methylation_value"))

hits_result_quin_unique_positions_overlap_GCs <- check_range_overlap(hits_result_quin_unique_positions, df_smash_meth, "repeat_start_minus_300", "repeat_end_plus_300", "CG_start", "CG_end", c("read_depth", "methylation_value"))


# Add column on whether or not the CpG is 50% methylated.
hits_result_trip_unique_positions_overlap_GCs$methylated_50 <- ifelse(hits_result_trip_unique_positions_overlap_GCs$overlap_methylation_value > 0.34 & hits_result_trip_unique_positions_overlap_GCs$overlap_methylation_value < 0.66, yes = "meth_50", no = "meth_not_50")

hits_result_quad_unique_positions_overlap_GCs$methylated_50 <- ifelse(hits_result_quad_unique_positions_overlap_GCs$overlap_methylation_value > 0.34 & hits_result_quad_unique_positions_overlap_GCs$overlap_methylation_value < 0.66, yes = "meth_50", no = "meth_not_50")

# count the number of methylated/unmethylated CpGs within 300bp per repeat that has a methylated CCGG within 300bp of either end of the repeat.
hits_result_trip_unique_positions_overlap_GCs_summarized <- hits_result_trip_unique_positions_overlap_GCs[,c("contig", "repeat_start", "repeat_end", "repeat_start_minus_300", "repeat_end_plus_300", "methylated_50")] %>% dplyr::group_by(repeat_start_minus_300, methylated_50) %>% dplyr::mutate(counts = n()) %>% unique()

hits_result_quad_unique_positions_overlap_GCs_summarized <- hits_result_quad_unique_positions_overlap_GCs[,c("contig", "repeat_start", "repeat_end", "repeat_start_minus_300", "repeat_end_plus_300", "methylated_50")] %>% dplyr::group_by(repeat_start_minus_300, methylated_50) %>% dplyr::mutate(counts = n()) %>% unique()




# Define the genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Function to get genomic sequence with extended range based on start position
get_genomic_sequence <- function(chromosome, start, extend_start = 20, extend_stop = 100) {
  # Adjust start and stop positions
  adj_start <- start - extend_start
  adj_stop <- start + extend_stop
  
  # Ensure the chromosome is in the correct format (e.g., "chr1")
  chr <- chromosome
  
  # Retrieve the sequence
  seq <- getSeq(genome, names = chr, start = adj_start, end = adj_stop)
  
  return(as.character(seq))
}

# Function to add sequence column to data frame
add_sequence_column <- function(df, chromosome_col, start_col) {
  df$sequence <- mapply(function(start) {
    get_genomic_sequence(df[[chromosome_col]][1], start)
  }, df[[start_col]])
  
  return(df)
}

# add sequence column. I will use the 10 bases 5' to the repeat to grep the repeat sequences from the GTEx fastq files.
hits_result_trip_unique_positions_overlap_GCs_to_grep <- add_sequence_column(unique(hits_result_trip_unique_positions_overlap_GCs[,c("contig", "repeat_start"),]), "contig", "repeat_start")
hits_result_quad_unique_positions_overlap_GCs_to_grep <- add_sequence_column(unique(hits_result_quad_unique_positions_overlap_GCs[,c("contig", "repeat_start"),]), "contig", "repeat_start")
hits_result_quin_unique_positions_overlap_GCs_to_grep <- add_sequence_column(unique(hits_result_quin_unique_positions_overlap_GCs[,c("contig", "repeat_start"),]), "contig", "repeat_start")

# Export tables for grep'ing.
write.table(hits_result_trip_unique_positions_overlap_GCs_to_grep, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/trip_repeats_sequences_to_grep.tsv", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(hits_result_quad_unique_positions_overlap_GCs_to_grep, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_sequences_to_grep.tsv", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(hits_result_quin_unique_positions_overlap_GCs_to_grep, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quin_repeats_sequences_to_grep.tsv", quote = F, row.names = F, col.names = F, sep = "\t")

#### Export tables for plotting ####
trip_repeats_meta <- unique(hits_result_trip_unique_positions_overlap_GCs[,c("contig","repeat_start","repeat_end","repeat_start_minus_300","repeat_end_plus_300")])
quad_repeats_meta <- unique(hits_result_quad_unique_positions_overlap_GCs[,c("contig","repeat_start","repeat_end","repeat_start_minus_300","repeat_end_plus_300")])

write.table(trip_repeats_meta, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/trip_repeats_meta.tsv", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(quad_repeats_meta, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_meta.tsv", quote = F, row.names = F, col.names = F, sep = "\t")
