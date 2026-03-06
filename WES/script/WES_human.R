
# WHOLE EXOME SEQUENCING ANALYSIS
# Author: Najneen Rejwana
# Dataset: SRR10341617
# Source: https://www.ebi.ac.uk/ena/browser/view/SRR10341617
# Platform: Illumina HiSeq 2500 (Paired-end)
# VCF: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266735
#-------------------------------------------------------------------------------------

# 1. SET WORKING DIRECTORY-----------------------------------------------------------


setwd("F:/github/GENOMICS/whole_exome_seq/script")
getwd()

R.version.string

# 2. INSTALL REQUIRED PACKAGES----------------------------------------------------------

# Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
install.packages(c(
  "fastqcr",
  "ggplot2",
  "tidyverse",
  "vcfR",
  "progressr"
))

# Bioconductor packages
BiocManager::install(c(
  "ShortRead",
  "Rsubread",
  "VariantAnnotation",
  "GenomicRanges",
  "BSgenome.Hsapiens.UCSC.hg38",
  "edgeR",
  "DESeq2",
  "biomaRt"
))

# 3. LOAD LIBRARIES
###############################

library(ShortRead)
library(fastqcr)
library(Rsubread)
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(vcfR)
library(tidyverse)
library(progressr)

# Enable progress bar
handlers(global = TRUE)
handlers("txtprogressbar")

# 4. LOAD VCF FILE-----------------------------------------------------------------------
vcf_file <- "F:/github/GENOMICS/whole_exome_seq/data/GSM8254208_R1.intersect.ann.vcf"
vcf <- read.vcfR(vcf_file)


# 5. EXPLORE VCF STRUCTURE---------------------------------------------------------------
summary(vcf)
str(vcf)
print(vcf)
head(vcf)

# Metadata
vcf@meta

# 6. DEPTH OF COVERAGE ANALYSIS-----------------------------------------------------------
class(coverage)
#[1] "standardGeneric"
#attr(,"package")
#[1] "methods"

#That means coverage is not your calculated vector.

# Extract depth values
dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Calculate mean coverage per variant
coverage_values <- rowMeans(dp, na.rm = TRUE)

#plot-------
library(ggplot2)

coverage_df <- data.frame(coverage = coverage_values)

ggplot(coverage_df, aes(x = coverage)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  labs(
    title = "Depth of Coverage Distribution",
    x = "Coverage Depth",
    y = "Frequency"
  )





# 7. CONVERT VCF TO DATAFRAME-------------------------------------------------------------------

vcf_df <- as.data.frame(vcf@fix)

# Extract annotation field
info_df <- vcfR::extract_info_tidy(
  vcf,
  info_fields = c("ANN")
)





# 8. PARSE SNPEFF ANNOTATIONS------------------------------------

mutations <- info_df %>%
  separate(
    ANN,
    into = c(
      "Allele",
      "Consequence",
      "Impact",
      "Gene",
      "Feature",
      "Biotype",
      "Rank",
      "HGVS_C",
      "HGVS_P",
      "cDNA_pos",
      "CDS_pos",
      "Protein_pos",
      "Distance",
      "Errors"
    ),
    sep = "\\|",
    extra = "drop"
  )

head(mutations)

# 9. SAVE MUTATION TABLE------------------------

write.csv(
  mutations,
  "mutations_output.csv",
  row.names = FALSE
)

file.info("mutations_output.csv")


# 10. VISUALIZATION-----------------------------------------------------------------

# Mutation consequences-----------------------------------------------------
ggplot(mutations, aes(x = Consequence)) +
  geom_bar(fill = "purple") +
  labs(
    title = "Mutation Consequences Distribution",
    x = "Consequence",
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mutation impact--------------------------------------------------------------
ggplot(mutations, aes(x = Impact)) +
  geom_bar(fill = "lightgreen") +
  labs(
    title = "Mutation Impact Distribution",
    x = "Impact",
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gene mutation frequency-----------------------------------------------------------
top_genes <- mutations %>%
  count(Gene, sort = TRUE) %>%
  head(20)

ggplot(top_genes, aes(x = reorder(Gene, n), y = n)) +
  geom_bar(stat = "identity", fill = "salmon") +
  coord_flip() +
  labs(
    title = "Top 20 Mutated Genes",
    x = "Gene",
    y = "Mutation Count"
  )

# Mutation types per gene (top 20 genes)---------------------------------------------------------
library(dplyr)
library(ggplot2)

# Get top 20 mutated genes
top_genes <- mutations %>%
  count(Gene, sort = TRUE) %>%
  slice_head(n = 20)

# Filter dataset for those genes
mutations_top <- mutations %>%
  filter(Gene %in% top_genes$Gene)

# Plot
ggplot(mutations_top, aes(x = reorder(Gene, Gene, length), fill = Consequence)) +
  geom_bar() +
  coord_flip() +
  labs(
    title = "Mutation Types in Top 20 Mutated Genes",
    x = "Gene",
    y = "Mutation Count"
  ) +
  theme_minimal()

# Mutation Spectrum-----------------------------------------------------

ref_alt <- paste(vcf@fix[,4], ">", vcf@fix[,5])

ggplot(data.frame(mutation = ref_alt), aes(x = mutation)) +
  geom_bar(fill = "darkred") +
  labs(
    title = "Mutation Spectrum",
    x = "Base Substitution",
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 45))


# Variant Distribution Across Chromosomes--------------------------------------------------
chrom_df <- as.data.frame(vcf@fix)

ggplot(chrom_df, aes(x = CHROM)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Variant Distribution Across Chromosomes",
    x = "Chromosome",
    y = "Variant Count"
  )