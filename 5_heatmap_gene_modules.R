library(ggplot2)
library(scales)
library(reshape2)

#### load data from previous steps ####
lnames = load(file = "./Routput/4_geneOverlap.RData");

dim(brca_gbm_gene_overlap)
dim(gbm_ov_gene_overlap)
dim(ov_brca_gene_overlap)

########### Generate heatmaps for each pairwise analysis of gene modules ###########
generateGeneModuleHeatmap <- function(geneOverlapObj, title, xl, yl, prefix)
{
  positionGeneList = geneOverlapObj$positionList
  pval_matrix = geneOverlapObj$modulePairMatrix
  v = positionGeneList[order(sapply(positionGeneList, function(x) x[1], simplify = TRUE), decreasing = FALSE)][1]
  idx_i = v[[1]][3]
  idx_j = v[[1]][4]
  melted_mat <- melt(pval_matrix)

  png(file = paste0(prefix, "_heatmap_gene_module_pairs.png"), width = 1200, height = 1000);
  par(cex.main = 1.5)

  p <- ggplot(data = melted_mat, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(colour = "gray80") +
    theme_bw(10) +
    ggtitle(title) +
    xlab(xl) +
    ylab(yl) +
    scale_fill_gradient2(low = muted("darkblue"), mid = "white", high = muted("red"),
                         midpoint = 0, space = "Lab", na.value = "grey10", guide = "colourbar",
                         limits = c(-0.049, 0.049)) +
    geom_rect(size = 10, colour = "orange", mapping = aes_string(x = "Var1", y = "Var2"),
              xmin = idx_i, xmax = idx_i, ymin = idx_j, ymax = idx_j) +
    theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  dev.off()
}


generateGeneModuleHeatmap(brca_gbm_gene_overlap,
                          "Overlap of Gene Modules - BRCA and GBM",
                          "BRCA gene modules", "GBM gene modules", "./results/5_BRCA_GBM")

generateGeneModuleHeatmap(gbm_ov_gene_overlap,
                          "Overlap of Gene Modules - GBM and OV",
                          "GBM gene modules", "OV gene modules", "./results/5_GBM_OV")

generateGeneModuleHeatmap(ov_brca_gene_overlap,
                          "Overlap of Gene Modules - OV and BRCA",
                          "OV gene modules", "BRCA gene modules", "./results/5_OV_BRCA")
