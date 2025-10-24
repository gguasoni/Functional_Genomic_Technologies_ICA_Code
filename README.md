## Task
Students were randomly assigned an RNA-seq study from Gene Expression Omnibus (GEO) to analyse. Students were initially presented with fastq files and had to conduct quality control checks to determine if the files needed to be trimmed, and if so, were trimmed by the course instructor using Trimmomatic and then subsequently asligned using STAR. Students were then tasked to perform differential gene expression analysis and functional enirchment analysis using their assigned data files using R/RStudio and produce an 8-page report of their analysis. 

## Script Information
- fgt_assignment.R
  - This script processes RNA-Seq data for quality control, differential expression analysis, functional enrichment, and visualization. It integrates various bioinformatics packages (DESeq2, edgeR, limma, etc.) to perform quality checks, normalization, differential expression (DE) analysis, and gene annotation. It also provides visualization tools such as boxplots, PCA plots, heatmaps, and functional enrichment analysis to facilitate downstream biological discovery.
