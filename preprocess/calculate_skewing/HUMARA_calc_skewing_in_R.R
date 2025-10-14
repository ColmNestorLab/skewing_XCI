library(data.table)
library(tidyverse)

dataz <- fread("data/data_rerun_by_Bjorn/HUMARA_data_blue_channel.csv")
dataz <- dataz %>% tidyr::separate(2, sep = "_", into = c("digest", "sample", "extended_sample"), remove = F)

dataz$gene <- "LOL"
dataz[which(dataz$Size > 225 & dataz$Size < 320),]$gene <- "AR"

data_filt <- dataz[dataz$gene != "LOL",]

metadataz <- unique(fread("data/data_rerun_by_Bjorn/HUMARA_metadata.csv", header = T)[,c("filename", "date_run")])
colnames(metadataz) <- c("Sample File Name", "date_run")

# add metadata
data_filt_with_meta <- merge(data_filt, metadataz, by = "Sample File Name")

# merge sample and date_run column of each sample.
data_filt_with_meta$sample_date <- paste0(data_filt_with_meta$sample, "_", data_filt_with_meta$date_run)

# get control samples
data_filt_with_meta_controls <- data_filt_with_meta[data_filt_with_meta$sample == "S",]

# find the latest run of each sample.
sample_dates_to_keep <- unique(data_filt_with_meta[,c("date_run", "sample")]) %>%
  group_by(sample) %>%
  slice_max(date_run) %>%
  ungroup()

# merge sample and date_run column of each sample I want to keep.
sample_dates_to_keep$sample_date <- paste0(sample_dates_to_keep$sample, "_", sample_dates_to_keep$date_run)

# Keep only latest run of each sample.
data_filt_with_meta_keep <- data_filt_with_meta[data_filt_with_meta$sample_date %in% sample_dates_to_keep$sample_date & data_filt_with_meta$sample != "S",]

# double-check there is only one entry per sample.
#hue <- unique(data_filt_with_meta_keep[,c("date_run", "sample")]) %>% dplyr::group_by(sample) %>% dplyr::count()

# add control samples again.
data_filt_with_meta_keep_fin <- rbind(data_filt_with_meta_keep,data_filt_with_meta_controls )

# Find the two peaks for each peak region, for each treatment.
two_highest_peaks_per_undigested_sample <- data_filt_with_meta_keep_fin[data_filt_with_meta_keep_fin$digest == "U",] %>% 
  dplyr::group_by(sample, gene, date_run) %>%
  dplyr::slice_max(order_by = Area, n = 2) %>%
  dplyr::ungroup()

# get data frame to match undigested values to.
digested_samples <- data_filt_with_meta_keep_fin[data_filt_with_meta_keep_fin$digest == "D",]

# Function to find the closest value
find_closest <- function(size, comparison_sizes) {
  comparison_sizes[which.min(abs(comparison_sizes - size))]
}

# Match closest Sizes
result <- two_highest_peaks_per_undigested_sample %>%
  rowwise() %>%
  mutate(closest_Size = {
    # Extract current sample and gene
    current_sample <- sample
    current_gene <- gene
    
    # Subset digested_samples for the same sample and gene
    subset_digested_samples <- digested_samples[digested_samples$Area > 300,] %>% filter(sample == current_sample & gene == current_gene)
    
    # Find the closest value from df2, handle empty subsets
    if (nrow(subset_digested_samples) > 0) {
      find_closest(Size, subset_digested_samples$Size)
    } else {
      NA  # Return NA if no matching sample and gene in df2
    }
  }) %>%
  ungroup()

# match magic
result$match_magic <- paste0(result$sample, "_", result$gene, "_", result$closest_Size)
digested_samples$match_magic <- paste0(digested_samples$sample, "_", digested_samples$gene, "_", digested_samples$Size)

# It works up until here. Just continue selecting and appending the additional data from the digested samples based on the Size and gene column in the new 'result' data frame. From there I should easily enough be able to wrangle the data to fit the below scripts.

# get full digested data
digested_data <- digested_samples[digested_samples$match_magic %in% result$match_magic,]

# rbind relevant columns
smash <- rbind(result[,c("sample","gene","digest", "Size", "Area", "Height", "match_magic", "date_run")], 
               digested_data[,c("sample","gene","digest", "Size", "Area", "Height", "match_magic", "date_run")])

# Here I add a column for each gene, digest and sample which specifies which of the two values are the highest and the lowest values. This makes making the data frame wide easier in the next step.
df <- smash %>%
  group_by(gene, digest, sample, date_run ) %>%
  mutate(high_low = ifelse(Size == max(Size), "HIGH", "LOW")) %>%
  ungroup()

# convert to data frame
df <- as.data.frame(df)

# Here I find measurements where the peak of the area goes out of bounds in the Thermo Fisher tool (i.e. values above 32000).
stored_over <- df[df$Height > 32000,]

# make the data frame wide format. Exclude control sample from here on out.
df_wide <- df[!df$sample %in% c("S", "C", "PC"),] %>%
  tidyr::unite("gene_status", digest, gene, high_low, sep = "_") %>%
  pivot_wider(
    names_from = gene_status, 
    values_from = 'Area', 
    id_cols = c(sample)
  )


# calculate skewing (from PMC9681199)
#Allele Ratio Mock Digestion (Rm)=allele 1 peak height / allele 2 peak height
#Allele Ratio HpaII Digestion (Rh)=allele 1 peak height / allele 2 peak height
#Normalized Ratio (Rn)=Rh/Rm
#XCI percentage = [Rn/(Rn +1)] * 100

df_wide_calc <- 
  df_wide %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(
    AR_ratio_mock = sum(as.numeric(U_AR_HIGH) / as.numeric(U_AR_LOW), na.rm = T),
    AR_ratio_digest = sum(as.numeric(D_AR_HIGH) / as.numeric(D_AR_LOW), na.rm = T)
    )

# get column names in which I want to change anything above 2.5 to NA (and thus excluding homozygous repeats).
ratio_mock_column_names <- c("AR_ratio_mock")

# Store rows with extremely low or high ratios, this because I want to know the repeat length regardless of whether the individual is homozygous for the repeat or not.
rows_with_low_high <- as.data.frame(df_wide_calc[apply(df_wide_calc[ratio_mock_column_names], 1, function(x) any(x>2.5|x<0.4)), c("sample", ratio_mock_column_names)])

# melt for easier handling
rows_with_low_high <- melt(setDT(df_wide_calc[,c("sample", ratio_mock_column_names)]), id.vars = "sample")
rows_with_low_high <- rows_with_low_high[rows_with_low_high$value < 0.4 | rows_with_low_high$value > 2.5,]

# decide which one to keep, i.e. which repeat lenght is correct.
rows_with_low_high$high_or_low_homozygous <- ifelse(rows_with_low_high$value > 2.5, yes = "keep_high", no = "keep_low")

# split variable column
rows_with_low_high <- rows_with_low_high %>% tidyr::separate(variable, into = c("gene"), sep = "_")

# Replace values > 2.5 with NA in matching columns
df_wide_calc[ratio_mock_column_names] <- lapply(df_wide_calc[ratio_mock_column_names], function(col) {
  col[col > 2.5 | col < 0.4] <- NA
  return(col)
})

df_wide_calc_final <- 
  df_wide_calc %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(
    AR_ratio_normalized = sum(AR_ratio_digest/AR_ratio_mock),
    AR_skewing =  round(abs(AR_ratio_normalized/(AR_ratio_normalized + 1)) * 100,1),
    AR_skewing =  ifelse(AR_skewing > 50, yes = AR_skewing, no = 100-AR_skewing)
  ) %>% dplyr::select(-AR_ratio_normalized)


##### Get repeat lengths ####
# Here I convert the size to repeat count, using a male sample with known repeat lengths.
sizing_standard <- fread("anno/sizing_standard_HUMARA.tsv")

# control sample peaks.
df_control <- df[df$sample == "S" & df$high_low == "HIGH",]

# append the repeat length to the full data frame.
repeat_df <- merge(df[df$sample != "S" & df$digest == "U",], df_control[,c("gene", "Size", "date_run")], by = c("gene", "date_run"), suffixes = c("_sample","_repeat"), all = T)

# add the size_standard df
repeat_df_meta <- merge(repeat_df, sizing_standard, by = "gene")

# calculate the difference in bases compared to the known repeat length sample.
repeat_df_meta$base_difference <- repeat_df_meta$Size_sample-repeat_df_meta$Size_repeat

# Get how many repeats that difference equates to. Round it.
repeat_df_meta$repeat_difference <- round(repeat_df_meta$base_difference/repeat_df_meta$repeat_type, 2)

# Calculate the no. of repeats for each allele. If the value is lower than 0 (i.e. negative), then absolute it and subtract. If positive just add to repeat count.
repeat_df_meta$sample_repeat_count <- ifelse(repeat_df_meta$repeat_difference < 0, 
                                             repeat_df_meta$repeats - abs(repeat_df_meta$repeat_difference), 
                                             repeat_df_meta$repeats + repeat_df_meta$repeat_difference)

# round it
repeat_df_meta$sample_repeat_count_rounded <- round(repeat_df_meta$sample_repeat_count, 0)

# remove excess columns
repeat_df_fin <- unique(repeat_df_meta[!is.na(repeat_df_meta$sample),c("sample", "gene", "sample_repeat_count_rounded", "high_low")])

# add which repeat length is correct for the homozygous samples.
repeat_df_fin_smash <- merge(repeat_df_fin, rows_with_low_high, by = c("sample", "gene"), all=T)

# change the repeat length to show which of the lengths that are correct for the homozygous samples.
try(repeat_df_fin_smash[which(repeat_df_fin_smash$high_or_low_homozygous == "keep_high" & repeat_df_fin_smash$high_low == "LOW"),]$sample_repeat_count_rounded <- "homozygous", silent = T)
try(repeat_df_fin_smash[which(repeat_df_fin_smash$high_or_low_homozygous == "keep_low" & repeat_df_fin_smash$high_low == "HIGH"),]$sample_repeat_count_rounded <- "homozygous", silent = T)

# make wide.
wide_dt <- dcast(setDT(repeat_df_fin_smash), sample ~ gene + high_low, value.var = "sample_repeat_count_rounded")

# change colnames
names(wide_dt) <- gsub("HIGH", "allele1", names(wide_dt))
names(wide_dt) <- gsub("LOW", "allele2", names(wide_dt))

# finalize table
final_table <- merge(df_wide_calc_final, wide_dt, by = "sample")

# As a final step, change skewing into "over" if the value is above 32000, as stored in 'stored_over'. If its labeled as homozygous already then keep that classification. Also change values stored in stored_below to 'UNDER' as these peaks tend to be unreliable.
heh <- melt(setDT(final_table), id.vars = "sample")
heh_test <- heh #[!is.na(heh$value) & heh$value != "homozygous",]

#stored_over_test <- stored_over[,c("sample", "gene")]
#stored_over_test$gene <- paste0(stored_over_test$gene, "_skewing")
#stored_over_test$cat <- "OVER"

#stored_below_test <- stored_below[,c("sample", "gene")]
#stored_below_test$gene <- paste0(stored_below_test$gene, "_skewing")
#stored_below_test$cat <- "UNDER"  
#stored_over_under <- rbind(stored_over_test, stored_below_test)

#stored_over_under <- rbind(stored_over_test)

heh_test2 <- heh

# Change all skewing values with the "OVER" tag to NA.
try(heh_test2[which(heh_test2$cat == "OVER" & !is.na(heh_test2$value)),]$value <- "OVER", silent = T)
#try(heh_test2[which(heh_test2$cat == "UNDER" & !is.na(heh_test2$value)),]$value <- "UNDER", silent = T)

# make wide yet again
final_table_wide_fixed <- dcast(unique(heh_test2[,c("sample", "variable", "value")]), sample ~ variable, value.var = "value")

# finalize
final_table_wide_fixed <- final_table_wide_fixed[final_table_wide_fixed$sample != "S",]

# final finalization, change column order
final_table_wide_fixed <- dplyr::select(final_table_wide_fixed, "sample",
                                        "AR_skewing","AR_allele1","AR_allele2")

# Export
write.table(final_table_wide_fixed, "data/HUMARA_results.tsv", sep = "\t", quote = F, row.names = F)
