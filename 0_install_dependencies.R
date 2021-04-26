# Install R packages

# Install Bioconductor package manager and the core Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
library("BiocManager")
BiocManager::install()

# WGCNA - https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
BiocManager::install("WGCNA")

# topGO - https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
BiocManager::install("topGO")

# genefilter - http://bioconductor.org/packages/release/bioc/html/genefilter.html
BiocManager::install("genefilter")

# BioMart database - https://bioconductor.org/packages/release/bioc/html/biomaRt.html
BiocManager::install("biomaRt")

# graphviz library for plotting - https://www.bioconductor.org/packages/release/bioc/html/Rgraphviz.html
BiocManager::install("Rgraphviz")

# Differential expression analysis - http://bioconductor.org/packages/release/bioc/html/edgeR.html
BiocManager::install("edgeR")

# Test two sets of gene lists - https://www.bioconductor.org/packages/release/bioc/html/GeneOverlap.html
BiocManager::install("GeneOverlap")

BiocManager::install("venn")

BiocManager::install("reshape2")

BiocManager::install("ggplot2")
