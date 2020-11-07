# Exploring Gene Expression Data Using Weighted Gene Co-expression Network Analysis Across Multiple Cancer Types

## CS5228 Knowledge Discovery and Data Mining
A.V. AKILA RAVIHANSA PERERA – A0212216X

TRAN KHANH HUNG – A0212253W


## Method

- WGCNA

Weighted gene correlation network analysis (WGCNA) is a powerful method that uses a topological overlap module approach for constructing co-expression networks based on gene expression data. This method involves reconstructing gene co-expression modules and summarizing modules using module eigengenes (ME) and intramodular hub genes.

- Gene Enrichment and Pathway Analysis

Biologically interesting modules were identified using Fisher's exact test. The overlapping and union sets of genes from theses interesting gene module pairs were 
subjected to Gene Set Enrichment Analysis using topGO package).

## Dataset

Three gene expression datasets for three cancer types (GBM, OV, BRCA) were selected from TCGA (The Cancer Genome Atlas).

 - Glioblastoma Multiforme (GBM) gene expression by RNAseq
 https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HiSeqV2_PANCAN.gz
 
 - Ovarian Serous Cystadenocarcinoma (OV) gene expression by RNAseq
 https://tcga.xenahubs.net/download/TCGA.OV.sampleMap/HiSeqV2_PANCAN.gz
 
 - Breast Invasive Carcinoma (BRCA) gene expression by RNAseq
https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz


## Analysis

###  Data Exploration 

Scale Free Topology Model       |  Mean Connectivity      | Selected Soft Threshold
:-------------------------:|:-------------------------:|:-------------------------:
![](results/2_SFTM_Fit_GBM.png)  |  ![](results/2_Mean_Connectivity_GBM.png) | 9
![](results/2_SFTM_Fit_OV.png)  |  ![](results/2_Mean_Connectivity_OV.png) | 6
![](results/2_SFTM_Fit_BRCA.png)  |  ![](results/2_Mean_Connectivity_BRCA.png) | 12


### Clustering Tree

![GBM](results/2_Clustering_Tree_GBM.png)
<br/><br/>

![GBM](results/2_Clustering_Tree_OV.png)
<br/><br/>

![GBM](results/2_Clustering_Tree_BRCA.png)
<br/><br/>

### Gene Expression Network

Network Heatmap       |      Eigengene Adjacency Heatmap      |      Eigengene Dendrogram
:-------------------------:|:-------------------------:|:-------------------------:
![](results/3_Network_heatmap_GBM.png) |  ![](results/3_Eigengene_heatmap_GBM.png)  |  ![](results/3_Eigengene_dendrogram_GBM.png)
![](results/3_Network_heatmap_OV.png) |  ![](results/3_Eigengene_heatmap_OV.png)  |  ![](results/3_Eigengene_dendrogram_OV.png)
![](results/3_Network_heatmap_BRCA.png) |  ![](results/3_Eigengene_heatmap_BRCA.png) |  ![](results/3_Eigengene_dendrogram_BRCA.png)


### Overlapping Genes in Gene Modules

Item                     | BRCA and GBM               |  GBM and OV               | OV and BRCA
:-----------------------:|:-------------------------:|:-------------------------:|:-------------------------:
|||
--- Intersection --- |||
Name | [grey_cyan_genes.txt](results/4_BRCA_GBM_Intersection_P0.999999996738616_grey_cyan_genes.txt) | [lightcyan_brown_genes.txt](results/4_GBM_OV_Intersection_P0.999999978350642_lightcyan_brown_genes.txt) | [grey_black_genes.txt](results/4_OV_BRCA_Intersection_P0.999999999999291_grey_black_genes.txt) 
Gene module pair with highest p-value | 0.999999996738616 | 0.999999978350642 | 0.999999999999291
Gene count in top module pair | 14 | 11 | 35
topGO plot | [grey_cyan_topGO.pdf](results/4_BRCA_GBM_Intersection_grey_cyan_topGOPlot_fullnames.pdf) | [lightcyan_brown_topGO.pdf](results/4_GBM_OV_Intersection_lightcyan_brown_topGOPlot_fullnames.pdf) | [grey_black_topGO.pdf](results/4_OV_BRCA_Intersection_grey_black_topGOPlot_fullnames.pdf)
topGO analysis | [grey_cyan_topGO.csv](results/4_BRCA_GBM_Intersection_grey_cyan_summary_topGO_analysis.csv) | [lightcyan_brown_topGO.csv](results/4_GBM_OV_Intersection_lightcyan_brown_summary_topGO_analysis.csv) | [grey_black_topGO.csv](results/4_OV_BRCA_Intersection_grey_black_summary_topGO_analysis.csv)
|||
|||
--- Union --- |||
Name | [grey_cyan_txt](results/4_BRCA_GBM_Union_P0.999999996738616_grey_cyan_genes.txt) | [lightcyan_brown_genes.txt](results/4_GBM_OV_Union_P0.999999978350642_lightcyan_brown_genes.txt) | [grey_black_genes.txt](results/4_OV_BRCA_Union_P0.999999999999291_grey_black_genes.txt)
Gene module pair with highest p-value | 0.999999996738616 | 0.999999978350642 | 0.999999999999291
Gene count in top module pair | 6066 | 1949 | 3262
topGO plot | [grey_cyan_topGP.pdf](results/4_BRCA_GBM_Union_grey_cyan_topGOPlot_fullnames.pdf) | [lightcyan_brown_topGO.pdf](results/4_GBM_OV_Union_lightcyan_brown_topGOPlot_fullnames.pdf) | [grey_black_topGO.pdf](results/4_OV_BRCA_Union_grey_black_topGOPlot_fullnames.pdf)
topGO analysis | [grey_cyan_topGO.csv](results/4_BRCA_GBM_Union_grey_cyan_summary_topGO_analysis.csv) | [lightcyan_brown_topGO.csv](results/4_GBM_OV_Union_lightcyan_brown_summary_topGO_analysis.csv) | [grey_black_topGO.csv](results/4_OV_BRCA_Union_grey_black_summary_topGO_analysis.csv)


## Project Structure and Run Instructions

 - Download and extract datasets to `./data` directory
 - Install dependencies (R packages)
    - Run `0_install_dependencies.R`
    
 - Gene filtering
    - Run `1_wgcna_cluster.R`
    
 - Build gene expression network and identify gene modules
    - Run `2_module_detection.R`
    
 - Generate gene network plots 
    - Run `3_network_visualization.R`

 - Gene Enrichment Analysis
    - Run `4_gene_enrichment.R`
    
