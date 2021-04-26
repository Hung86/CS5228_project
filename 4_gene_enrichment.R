source("utility_functions.R")

#### load data from previous steps ####
lnames = load(file = "./Routput/1_dataInput.RData");
lnames = load(file = "./Routput/2_networkConstruction.RData");

dim(gbmExpr)
dim(ovExpr)
dim(brcaExpr)


#### find gene modules ####
geneModules_ToFile("./results/4_ME_BRCA", brcaMEs)
geneModules_ToFile("./results/4_ME_GBM", gbmMEs)
geneModules_ToFile("./results/4_ME_OV", ovMEs)


#### find overlap in gene module pairs across cancer types ####
brca_gbm_gene_overlap = findGeneOverlap(brcaExpr, gbmExpr, brcaModuleColors, gbmModuleColors, brcaMEs, gbmMEs)
gbm_ov_gene_overlap = findGeneOverlap(gbmExpr, ovExpr, gbmModuleColors, ovModuleColors, gbmMEs, ovMEs)
ov_brca_gene_overlap = findGeneOverlap(ovExpr, brcaExpr, ovModuleColors, brcaModuleColors, ovMEs, brcaMEs)


# Save gene overlap data
save(ov_brca_gene_overlap, gbm_ov_gene_overlap, brca_gbm_gene_overlap,
     file = "./Routput/4_geneOverlap.RData")


#### set background gene set ####
bg_genes = names(gbmExpr)


#### find significant gene module pairs (lowest p-value)
geneOverlap_ToFile("./results/4_GBM_OV_lowP", gbm_ov_gene_overlap, FALSE)
geneOverlap_VennDiagram("./results/4","GBM","OV", gbm_ov_gene_overlap,FALSE)
gbm_ov_topGOEnrichment_lowP = topGO_GeneEnrichment(bg_genes, gbm_ov_gene_overlap, FALSE)
topGO_GeneEnrichment_ToFile("./results/4_GBM_OV_lowP", gbm_ov_topGOEnrichment_lowP)

geneOverlap_ToFile("./results/4_OV_BRCA_lowP", ov_brca_gene_overlap, FALSE)
geneOverlap_VennDiagram("./results/4","OV","BRCA", ov_brca_gene_overlap, FALSE)

ov_brca_topGOEnrichment_lowP = topGO_GeneEnrichment(bg_genes, ov_brca_gene_overlap, FALSE)
topGO_GeneEnrichment_ToFile("./results/4_OV_BRCA_lowP", ov_brca_topGOEnrichment_lowP)

geneOverlap_ToFile("./results/4_BRCA_GBM_lowP", brca_gbm_gene_overlap, FALSE)
geneOverlap_VennDiagram("./results/4","BRCA","GBM", brca_gbm_gene_overlap, FALSE)
brca_gbm_topGOEnrichment_lowP = topGO_GeneEnrichment(bg_genes, brca_gbm_gene_overlap, FALSE)
topGO_GeneEnrichment_ToFile("./results/4_BRCA_GBM_lowP", brca_gbm_topGOEnrichment_lowP)
