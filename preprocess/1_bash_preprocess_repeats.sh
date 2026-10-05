# Downloaded methylation data:
#C003VO55	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/C003VO/CD8-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/C003VO55.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
#S001U353	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/S001U3/regulatory_T_cell/Bisulfite-Seq/CNAG/S001U353.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
#C003VO55	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/C003VO/CD8-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/C003VO55.CPG_methylation_calls.bs_cov.GRCh38.20160531.bw
#S001U353	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/S001U3/regulatory_T_cell/Bisulfite-Seq/CNAG/S001U353.CPG_methylation_calls.bs_cov.GRCh38.20160531.bw

# Used:
#S009W451	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/S009W4/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/S009W451.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
#S009W451	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/S009W4/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/S009W451.CPG_methylation_calls.bs_cov.GRCh38.20160531.bw
#P586	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/329.2/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/P586.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
#P586	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/329.2/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/P586.CPG_methylation_calls.bs_cov.GRCh38.20160531.bw
#S007DD51	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/S007DD/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/S007DD51.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
#S007DD51	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/S007DD/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/S007DD51.CPG_methylation_calls.bs_cov.GRCh38.20160531.bw
#P582	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/145.2/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/P582.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
#P582	http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/145.2/CD4-positive_alpha-beta_T_cell/Bisulfite-Seq/CNAG/P582.CPG_methylation_calls.bs_cov.GRCh38.20160531.bw

#############################################################################################
###### Identify CpGs on the X-chromosome that are 50% methylated (actually 34% to 66%) ######
#############################################################################################

# data locations
methylation_data_location=/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data
reference_fasta_location=/home/bjogy93/Desktop/TREX_proj/TREX_bash/reference_genomes
repeat_output_location=/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output

# softwares
run_bigwigtowig=/home/bjogy93/Desktop/TREX_proj/TREX_bash/software/bigWigToWig
run_wig2bed=/home/bjogy93/Desktop/TREX_proj/TREX_bash/software/bedops_linux_x86_64-v2.4.41/bin/wig2bed

# Pre-process the two files per sample, one with methylation values and one with read coverage (both from blueprint WGBS).
#cat /run/user/1000/gvfs/smb-share:server=nestornas.local,share=piranha_share1/TREX/Dora/TREX/new_pipeline/WGBS_bigwigs.txt | while read sample_id	file_name	type
#do

# Convert bigwig to wig and then to bed.
#$run_bigwigtowig -chrom=chrX $methylation_data_location/$file_name $methylation_data_location/$sample_id.$type.wig
#$run_wig2bed < $methylation_data_location/$sample_id.$type.wig > $methylation_data_location/$sample_id.$type.bed

# remove intermediate files.
#rm $methylation_data_location/$sample_id.$type.wig

#done


# Combine methylation score and coverage from the two bed files.
#cut -f1 /run/user/1000/gvfs/smb-share:server=nestornas.local,share=piranha_share1/TREX/Dora/TREX/new_pipeline/WGBS_bigwigs.txt | uniq | while read	sample_id
#do

# Sort the files by the first four columns, following sorting of the CpG position.
#join -j 2 -o 1.1,1.2,1.3,1.5,2.5 <(sort -k2 $methylation_data_location/$sample_id.cov.bed) <(sort -k2 $methylation_data_location/$sample_id.call.bed) > $methylation_data_location/$sample_id.merged.bed

#done

#bedops -u $methylation_data_location/C003VO55.merged.bed $methylation_data_location/P582.merged.bed $methylation_data_location/P586.merged.bed $methylation_data_location/S001U353.merged.bed $methylation_data_location/S007DD51.merged.bed $methylation_data_location/S009W451.merged.bed | sed 's/ /\t/g' >  $methylation_data_location/samples.merged.bed

#sort -k2 $methylation_data_location/samples.merged.bed > $methylation_data_location/sorted.samples.merged.bed

# merge the files. This was done in R because it is so much easier.
#library(data.table)
#library(tidyverse)
#df <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/sorted.samples.merged.bed")
#df_stats <- df %>% group_by(V1, V2, V3) %>% summarise(sum_coverage = sum(V4), mean_methylation = mean(V5))
#colnames(df_stats) <- c("contig", "CpG_start" , "CpG_stop", "sum_coverage", "mean_methylation")
#write.table(df_stats, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/summarized.sorted.samples.merged.bed", quote = F, row.names = F, col.names = F, sep = "\t")
#df_stats_50 <- df_stats[df_stats$mean_methylation > 0.34 & df_stats$mean_methylation < 0.66 & df_stats$sum_coverage > 4,]
#write.table(df_stats_50, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/50_meth_summarized.sorted.samples.merged.bed", quote = F, row.names = F, col.names = F, sep = "\t")
#df_stats_not_50 <- df_stats[!(df_stats$mean_methylation > 0.34 & df_stats$mean_methylation < 0.66) & df_stats$sum_coverage > 4,]
#write.table(df_stats_not_50, "/home/bjogy93/Desktop/TREX_proj/TREX_bash/blueprint_methylation_data/not_50_meth_summarized.sorted.samples.merged.bed", quote = F, row.names = F, col.names = F, sep = "\t")


############################################################################################
############ Identify triplet, quadruple, quintuple repeats on the X-chromosome ############
############################################################################################


# Identify triplet repeats

# Function to generate triplet repeat combinations with at least one C or G and at most two C or G
generate_triplet_repeats() {
    nucleotides=("A" "C" "G" "T")
    for i in "${nucleotides[@]}"; do
        for j in "${nucleotides[@]}"; do
            for k in "${nucleotides[@]}"; do
                triplet="$i$j$k"
                # Count the number of Cs and Gs
                count_cg=$(echo "$triplet" | grep -o '[CG]' | wc -l)
                # Include triplets with at least one and at most two Cs or Gs
                if [ "$count_cg" -ge 1 ] && [ "$count_cg" -le 2 ]; then
                    echo "$triplet"
                fi
            done
        done
    done
}

# Call the function and save the output to a file
output_file=$repeat_output_location/triplet_repeats.txt
generate_triplet_repeats > "$output_file"


# Get the location of all the repeats from the above function (requiring at least 6 repeats).
while read name
do
  zgrep -b -o -P -h -i "($name){6,}" $reference_fasta_location/chrX.fasta >> 2.txt
done  < $repeat_output_location/triplet_repeats.txt

# Sort
sort -k1n,1 2.txt | awk -F":" '{print  $1 "\t" length($2)+$1}'  | awk '{ if ($1 > old+1) print; old = $0; }' > $repeat_output_location/triplet_repeats_position.txt

# remove intermediate
rm 2.txt



# Identify quadruple repeats
base=( A C G T )

for a in "${base[@]}" 
do
for b in "${base[@]}" 
do
for c in "${base[@]}" 
do
for d in "${base[@]}" 
do
echo "$a" "$b" "$c" "$d"
done
done
done
done |  awk '{print $0 "\t" gsub(/C/, "") "\t" gsub(/G/, "")}' | 
	awk  '$7 = $5+$6' | 
	awk '{if ($7<3) print $1 $2 "\t" $3 $4 }' | 
	awk '{if ($1!=$2) print $1 $2}' > $repeat_output_location/quad_repeats.txt

# Get the location of all the repeats from the above function (requiring at least 6 repeats).
while read name
do
  zgrep -b -o -P -h -i "($name){6,}" $reference_fasta_location/chrX.fasta >> 2.txt
done  < $repeat_output_location/quad_repeats.txt

# Sort
sort -k1n,1 2.txt | awk -F":" '{print  $1 "\t" length($2)+$1}'  | awk '{ if ($1 > old+1) print; old = $0; }' > $repeat_output_location/quad_repeats_position.txt

# remove intermediate
rm 2.txt


# Identify quintuple repeats
base=( A C G T )

for a in "${base[@]}" 
do
for b in "${base[@]}" 
do
for c in "${base[@]}" 
do
for d in "${base[@]}" 
do
for e in "${base[@]}" 
do
echo "$a" "$b" "$c" "$d" "$e"
done
done
done
done
done |  awk '{print $0 "\t" gsub(/C/, "") "\t" gsub(/G/, "")}' | 
	awk  '$8 = $6+$7' | 
	awk '{if ($7<3) print $1 $2 $3 $4 $5}' > $repeat_output_location/quintuple_repeats.txt

# Get the location of all the repeats from the above function (requiring at least 6 repeats).
while read name
do
  zgrep -b -o -P -h -i "($name){6,}" $reference_fasta_location/chrX.fasta >> 2.txt
done  < $repeat_output_location/quintuple_repeats.txt

# Sort
sort -k1n,1 2.txt | awk -F":" '{print  $1 "\t" length($2)+$1}'  | awk '{ if ($1 > old+1) print; old = $0; }' > $repeat_output_location/quintuple_repeats_position.txt

# remove intermediate
rm 2.txt


# Create list of the position of all CG sites on chrX
grep -i -b -o 'CG' $reference_fasta_location/chrX.fasta | grep -oE '[0-9]+' > $repeat_output_location/CGpos.txt

# Create list of the position of all CCGG sites on chrX
grep -i -b -o 'CCGG' $reference_fasta_location/chrX.fasta | grep -oE '[0-9]+' > $repeat_output_location/CCGGpos.txt


# From here, move over to R to do lots of things.
# Script called do_stuff.R.
# After that, back to bash to do grep'ing of the sequences.




# Candidate triplet repeat sequences from the WGS .txt fastq files #
# Input TSV file containing DNA sequences
input_tsv="/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/trip_repeats_sequences_to_grep.tsv"

# Directory with FASTQ files converted to text
fastq_dir="/run/user/1000/gvfs/smb-share:server=nestornas.local,share=piranha_share1/TREX/GTEx_50_females_whole_genome_seq/fastq"

# Output directory for grep results
output_dir="/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/trip"
#mkdir -p "$output_dir"

# Read the DNA sequences and patterns from the TSV file
mapfile -t sequences_and_patterns < <(awk -F'\t' '{print $3"\t"$5}' "$input_tsv")

# Export variables for parallel execution
export fastq_dir
export output_dir

# Function to process each sequence and FASTQ file combination
process_file() {
  local line="$1"
  local file="$2"

  # Split the sequence and pattern
  local sequence=$(echo "$line" | cut -f1)
  local pattern=$(echo "$line" | cut -f2)

  local flanking_sequence="${sequence:0:10}"
  local file_name
  file_name=$(basename "$file" .fastq)
  
  # Replace any special characters in the pattern to make it file-name safe (optional step)
  safe_pattern=$(echo "$pattern" | sed 's/[^a-zA-Z0-9]/_/g')

  # Modify the output file name to include both sequence and pattern
  local grep_output_file="${output_dir}/${file_name}_${sequence}_${safe_pattern}_grep.txt"

  # Combine both grep patterns into a single grep command with extended regex (-E)
  grep -E "${flanking_sequence}.*${pattern}" "$file" > "$grep_output_file"
  
  echo "Grep results saved to $grep_output_file"
}

export -f process_file

# Run the process_file function in parallel, ensuring each sequence is matched with its corresponding pattern
parallel -j 100 process_file ::: "${sequences_and_patterns[@]}" ::: "$fastq_dir"/*.fastq

echo "Script completed."




# Candidate quadruplet repeat sequences from the WGS .txt fastq files #
# Input TSV file containing DNA sequences
input_tsv="/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_sequences_to_grep.tsv"

# Directory with FASTQ files converted to text
fastq_dir="/run/user/1000/gvfs/smb-share:server=nestornas.local,share=piranha_share1/TREX/GTEx_50_females_whole_genome_seq/fastq"

# Output directory for grep results
output_dir="/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/quad"
#mkdir -p "$output_dir"

# Read the DNA sequences and patterns from the TSV file
mapfile -t sequences_and_patterns < <(awk -F'\t' '{print $3"\t"$5}' "$input_tsv")

# Export variables for parallel execution
export fastq_dir
export output_dir

# Function to process each sequence and FASTQ file combination
process_file() {
  local line="$1"
  local file="$2"

  # Split the sequence and pattern
  local sequence=$(echo "$line" | cut -f1)
  local pattern=$(echo "$line" | cut -f2)

  local flanking_sequence="${sequence:0:10}"
  local file_name
  file_name=$(basename "$file" .fastq)
  
  # Replace any special characters in the pattern to make it file-name safe (optional step)
  safe_pattern=$(echo "$pattern" | sed 's/[^a-zA-Z0-9]/_/g')

  # Modify the output file name to include both sequence and pattern
  local grep_output_file="${output_dir}/${file_name}_${sequence}_${safe_pattern}_grep.txt"

  # Combine both grep patterns into a single grep command with extended regex (-E)
  grep -E "${flanking_sequence}.*${pattern}" "$file" > "$grep_output_file"
  
  echo "Grep results saved to $grep_output_file"
}

export -f process_file

# Run the process_file function in parallel, ensuring each sequence is matched with its corresponding pattern
parallel -j 100 process_file ::: "${sequences_and_patterns[@]}" ::: "$fastq_dir"/*.fastq

echo "Script completed."




# Candidate quintuplet repeat sequences from the WGS .txt fastq files #
# Input TSV file containing DNA sequences
input_tsv="/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quin_repeats_sequences_to_grep.tsv"

# Directory with FASTQ files converted to text
fastq_dir="/run/user/1000/gvfs/smb-share:server=nestornas.local,share=piranha_share1/TREX/GTEx_50_females_whole_genome_seq/fastq"

# Output directory for grep results
output_dir="/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/quin"
#mkdir -p "$output_dir"

# Read the DNA sequences and patterns from the TSV file
mapfile -t sequences_and_patterns < <(awk -F'\t' '{print $3"\t"$5}' "$input_tsv")

# Export variables for parallel execution
export fastq_dir
export output_dir

# Function to process each sequence and FASTQ file combination
process_file() {
  local line="$1"
  local file="$2"

  # Split the sequence and pattern
  local sequence=$(echo "$line" | cut -f1)
  local pattern=$(echo "$line" | cut -f2)

  local flanking_sequence="${sequence:0:10}"
  local file_name
  file_name=$(basename "$file" .fastq)
  
  # Replace any special characters in the pattern to make it file-name safe (optional step)
  safe_pattern=$(echo "$pattern" | sed 's/[^a-zA-Z0-9]/_/g')

  # Modify the output file name to include both sequence and pattern
  local grep_output_file="${output_dir}/${file_name}_${sequence}_${safe_pattern}_grep.txt"

  # Combine both grep patterns into a single grep command with extended regex (-E)
  grep -E "${flanking_sequence}.*${pattern}" "$file" > "$grep_output_file"
  
  echo "Grep results saved to $grep_output_file"
}

export -f process_file

# Run the process_file function in parallel, ensuring each sequence is matched with its corresponding pattern
parallel -j 100 process_file ::: "${sequences_and_patterns[@]}" ::: "$fastq_dir"/*.fastq

echo "Script completed."






# Candidate autosomal repeat sequences from the WGS .txt fastq files #
# Input TSV file containing DNA sequences
input_tsv="/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/autosomal_repeats_sequences_to_grep.tsv"

# Directory with FASTQ files converted to text
fastq_dir="/run/user/1000/gvfs/smb-share:server=nestornas.local,share=piranha_share1/TREX/GTEx_50_females_whole_genome_seq/fastq"

# Output directory for grep results
output_dir="/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/autosomal"
#mkdir -p "$output_dir"

# Read the DNA sequences and patterns from the TSV file
mapfile -t sequences_and_patterns < <(awk -F'\t' '{print $3"\t"$5}' "$input_tsv")

# Export variables for parallel execution
export fastq_dir
export output_dir

# Function to process each sequence and FASTQ file combination
process_file() {
  local line="$1"
  local file="$2"

  # Split the sequence and pattern
  local sequence=$(echo "$line" | cut -f1)
  local pattern=$(echo "$line" | cut -f2)

  local flanking_sequence="${sequence:0:10}"
  local file_name
  file_name=$(basename "$file" .fastq)
  
  # Replace any special characters in the pattern to make it file-name safe (optional step)
  safe_pattern=$(echo "$pattern" | sed 's/[^a-zA-Z0-9]/_/g')

  # Modify the output file name to include both sequence and pattern
  local grep_output_file="${output_dir}/${file_name}_${sequence}_${safe_pattern}_grep.txt"

  # Combine both grep patterns into a single grep command with extended regex (-E)
  grep -E "${flanking_sequence}.*${pattern}" "$file" > "$grep_output_file"
  
  echo "Grep results saved to $grep_output_file"
}

export -f process_file

# Run the process_file function in parallel, ensuring each sequence is matched with its corresponding pattern
parallel -j 100 process_file ::: "${sequences_and_patterns[@]}" ::: "$fastq_dir"/*.fastq

echo "Script completed."

