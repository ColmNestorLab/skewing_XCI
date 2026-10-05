library(QDNAseq)
library(QDNAseq.hg38) # installed with devtools and their github.
library(Biobase)

bam_directory <- "low_pass_bams/fastqs"
bam_files <- list.files(path = bam_directory, pattern = "*.bam$", full.names = TRUE)

##### 500k bins #####
bins <- getBinAnnotations(binSize = 500, genome = "hg38")  # Adjust the bin size if needed
readCounts <- binReadCounts(bins, bamfiles = bam_files)

binsChromosomes <- fData(readCounts)$chromosome
table(binsChromosomes)  # This will show the counts of bins for each chromosome

readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE, chromosomes = "Y")
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)

copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

pdf(file = "01_plots/qDNAseq_output_500k_bins_segment.pdf")
plot(copyNumbersSegmented)
dev.off()

pdf(file = "01_plots/qDNAseq_output_500k_bins.pdf")
plot(readCountsFiltered)
dev.off()