# Genomics Repository

This repository contains genomics-related datasets and analysis workflows.

#Currently, it includes resources for **Whole Exome Sequencing (WES)** analysis.


## Whole Exome Sequencing (WES) Pipeline – Brain Cancer Dataset
Dataset: SRR10341617 (Illumina HiSeq 2500, paired-end)
VCF: GSE266735

This pipeline performs variant discovery and analysis from whole exome sequencing reads using R and standard bioinformatics tools

# Workflow Steps

1) Set Working Directory & Install Packages
Prepare environment and install R/Bioconductor packages.

2) Load Libraries
Load all required packages for WES analysis.

3) Load VCF File
Read the annotated VCF file using vcfR.

4) Explore VCF Structure
Inspect metadata, variant fields, and genotype information.

5) Depth of Coverage Analysis

Extract DP (read depth) from VCF.
Compute mean coverage per variant.
Visualize coverage distribution.

6) Convert VCF to DataFrame
Extract variant information to a tabular format for downstream processing.

7) Parse ANNOVAR / SnpEff Annotations
Split the INFO ANN column into columns: Consequence, Impact, Gene, Protein change, etc.

8) Save Mutation Table
Export cleaned mutation information to CSV.

9) Visualization
- Generate the following plots:
- Mutation consequences distribution
- Mutation impact distribution
- Top 20 mutated genes
- Mutation types per gene (top 20)
- Mutation spectrum (base substitutions)
- Variant distribution across chromosomes
