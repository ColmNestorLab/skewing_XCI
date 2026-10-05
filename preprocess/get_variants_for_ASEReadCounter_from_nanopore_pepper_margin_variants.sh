#### SelectVariants ####
/home/bjogy93/Desktop/software/gatk-4.6.2.0/gatk SelectVariants \
-R /home/bjogy93/Desktop/indexes/gatk/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-V PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
-O filtered_PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
--exclude-filtered \
--create-output-variant-index \
--exclude-non-variants \
--exclude-filtered \
-restrict-alleles-to BIALLELIC


# Make vcfs to use with ASEReadCounter
# Keep only heterozygous variants
bcftools filter -i'GT="het"' filtered_PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz > fin.bcftools.vcf

# Extract WES readcounts for #
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %AD]\t[ %DP]\n' fin.bcftools.vcf | sed 's/,/\t/g' > ReadCount.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %AD]\t[ %DP]\n' filtered_PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz | sed 's/,/\t/g' > ReadCount_unfiltered.tsv

