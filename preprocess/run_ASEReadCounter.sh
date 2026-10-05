cat sample.txt | while read sample
do

echo BAM /home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/direct_RNAseq/bams/$sample.bam
echo variant /home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/vcfs/$sample'_'fin.bcftools.vcf
echo output /home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/ASEReadCounter/$sample.ASEReadCounter.table

### Run ASEReadCounter ###
/home/bjogy93/Desktop/software/gatk-4.6.2.0/gatk ASEReadCounter \
--reference /home/bjogy93/Desktop/indexes/gatk/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
--input /home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/direct_RNAseq/bams/$sample.bam \
--variant /home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/vcfs/$sample'_'fin.bcftools.vcf \
--output /home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/revision_new_experiment/ASEReadCounter/$sample.ASEReadCounter.table \
--output-format RTABLE \
--min-base-quality 20 \
-L chrX

done
