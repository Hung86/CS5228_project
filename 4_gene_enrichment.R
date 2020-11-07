library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(igraph)
library(edgeR)
library(WGCNA)
library(GeneOverlap)
source("utility_functions.R")



lnames = load(file = "./Routput/2-networkConstruction.RData");
lnames = load(file = "./Routput/1-dataInput.RData");


dim(gbmExpr)
dim(ovExpr)
dim(brcaExpr)

####################Operation for all datasets################################################
geneModules_ToFile("./results/4_GBM", gbmMEs)
geneModules_ToFile("./results/4_OV", ovMEs)
geneModules_ToFile("./results/4_BRCA", brcaMEs)


####High P-value
bg_genes = names(gbmExpr)

gbm_ov_gene_overlap = findGeneOverlap(gbmExpr, ovExpr, gbmModuleColors, ovModuleColors, gbmMEs, ovMEs)
geneOverlap_ToFile("./results/4_GBM_OV", gbm_ov_gene_overlap)
gbm_ov_topGOEnrichment = topGO_GeneEnrichment(bg_genes, gbm_ov_gene_overlap)
topGO_GeneEnrichment_ToFile("./results/4_GBM_OV", gbm_ov_topGOEnrichment)

ov_brca_gene_overlap = findGeneOverlap(ovExpr, brcaExpr, ovModuleColors, brcaModuleColors, ovMEs, brcaMEs)
geneOverlap_ToFile("./results/4_OV_BRCA", ov_brca_gene_overlap)
ov_brca_topGOEnrichment = topGO_GeneEnrichment(bg_genes, ov_brca_gene_overlap)
topGO_GeneEnrichment_ToFile("./results/4_OV_BRCA", ov_brca_topGOEnrichment)

brca_gbm_gene_overlap = findGeneOverlap(brcaExpr, gbmExpr, brcaModuleColors, gbmModuleColors, brcaMEs, gbmMEs)
geneOverlap_ToFile("./results/4_BRCA_GBM", brca_gbm_gene_overlap)
brca_gbm_topGOEnrichment = topGO_GeneEnrichment(bg_genes, brca_gbm_gene_overlap)
topGO_GeneEnrichment_ToFile("./results/4_BRCA_GBM", brca_gbm_topGOEnrichment)


####Low P-value
geneOverlap_ToFile("./results/4_GBM_OV_lowP", gbm_ov_gene_overlap, FALSE)
gbm_ov_topGOEnrichment_lowP = topGO_GeneEnrichment(bg_genes, gbm_ov_gene_overlap, FALSE)
topGO_GeneEnrichment_ToFile("./results/4_GBM_OV_lowP", gbm_ov_topGOEnrichment_lowP)

geneOverlap_ToFile("./results/4_OV_BRCA_lowP", ov_brca_gene_overlap, FALSE)
ov_brca_topGOEnrichment_lowP = topGO_GeneEnrichment(bg_genes, ov_brca_gene_overlap, FALSE)
topGO_GeneEnrichment_ToFile("./results/4_OV_BRCA_lowP", ov_brca_topGOEnrichment_lowP)

geneOverlap_ToFile("./results/4_BRCA_GBM_lowP", brca_gbm_gene_overlap, FALSE)
brca_gbm_topGOEnrichment_lowP = topGO_GeneEnrichment(bg_genes, brca_gbm_gene_overlap, FALSE)
topGO_GeneEnrichment_ToFile("./results/4_BRCA_GBM_lowP", brca_gbm_topGOEnrichment_lowP)










