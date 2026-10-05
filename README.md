# TRiXi: A Multi-Target Tandem Repeat-Based Method for Accurate Detection of X-Inactivation Skewing in Humans

*Dóra Goldmann<sup>1,5</sup>, Björn Gylemo<sup>1,5</sup>, Maike Bensberg<sup>1,5</sup>, Ingela Johansson<sup>1</sup>, Svenja Fornahl<sup>1</sup>, Júlia Goldmann<sup>1</sup>, Lisa Haglund<sup>1</sup>, Magdalena Vitkova<sup>1</sup>, Janos Kondri<sup>1</sup>, Sandra Hellberg<sup>1</sup>, Engla Haglund<sup>1</sup>, Johanna Ungerstedt<sup>2,3,4</sup>, Charlotta Dabrosin<sup>3</sup>, Shadi Jafari<sup>1,6</sup>, Johnny Ludvigsson<sup>1,6</sup> & Colm E. Nestor<sup>1,6</sup>*

<sup>1</sup> Division of Children's and Women's Health, Department of Biomedical and Clinical Sciences, Faculty of Medicine and Health Sciences, Linköping University, Linköping, Sweden

<sup>2</sup> Hematology Clinic, Linköping University Hospital, Linköping, Sweden

<sup>3</sup> Division of Surgery, Orthopedics and Oncology (KOO), Department of Biomedical and Clinical Sciences, Faculty of Medicine and Health Sciences, Linköping University, Linköping, Sweden

<sup>4</sup> Science for Life Laboratory, Linköping University, Linköping, Sweden

<sup>5</sup> joint first authorship

<sup>6</sup> joint senior authorship

Corresponding author: colm.nestor@liu.se

##

Here you can find scripts used to perform analyses and generate figures in the above publication (link to be added). Author of bioinformatic analysis scripts: Björn Gylemo.

## Abstract
Preferential usage of either X-chromosome, skewed X-inactivation (sXCI), can profoundly alter the penetrance and expressivity of X-linked traits in females, often with life-threatening consequences. Despite its clinical relevance, the true frequency, origins, and stability of sXCI in human tissues remain poorly understood, in part because many existing methods cannot determine XCI status in a substantial proportion of females, limiting their applicability in large-scale population studies. 

Here, we present TRiXi (Tandem Repeat-based Identification of X-Inactivation), a simple and highly sensitive technique capable of detecting sXCI from as little as 10 ng of archived human DNA. Applying TRiXi to over 1,000 neonatal cord blood samples from healthy infants, we demonstrate that pronounced sXCI at birth is common in humans. Use of longitudinal samples from the same donors revealed that sXCI patterns in skewed neonates were stably maintained into adolescence, whereas in non-skewed females, patterns fluctuated considerably over time. 

We further demonstrate that TRiXi accurately assesses skewing by comparison to orthogonal approaches, including allele-specific expression, allele-specific DNA methylation and the HUMARA assay. Using low-pass short-read sequencing, whole-genome long-read sequencing, and targeted exon sequencing, we exclude most known genetic drivers of sXCI in extremely skewed infants, suggesting contributions from stochastic processes or currently unknown genetic modifiers. Further, by inferring X-inactivation patterns in 4,571 adult tissues from 258 females in the GTEx cohort, we demonstrate that skewed XCI is not restricted to blood, revealing widespread sXCI across non-hematopoietic tissues in females. Our study represents the most comprehensive analysis of neonatal sXCI to date and reveals sXCI as a widespread and stable epigenetic trait. This can profoundly modify the penetrance of X-linked conditions with immediate implications for genetic diagnostics and our understanding of human disease variability. 

TRiXi provides a convenient and accessible approach to inform diagnosis of X-linked traits by assaying sXCI. Studying individuals with sXCI may yield critical insights into (i) the initiation of XCI and its underlying genomic elements, (ii) previously unrecognized lethal X-linked traits and (iii) improved polygenic risk prediction through incorporation of X-chromosome variation.
##

## Instructions
*Code has been tested on Ubuntu 22.04.5 LTS and R version 4.6.1 (2026-06-24)*

*Hardware used AMD Ryzen Threadripper PRO 5995WX 64-Cores, 256 GB of RAM. However, should work on any system with enough RAM to support the in-memory operations in R.*

*To run the code, download the scripts and install the required R packages. See package_versions.tsv for packages and package version or use the renv.lock file.*

*Data used in this paper includes: WGS data from GTEx to identify suitable X-linked repeats, tables exported from GeneMapper5 (from 3500 Genetic Analyzer runs), RNA-seq and WES from GTEx, long-read WGS (ONT), direct RNA-seq data, TRiXi results.*  
