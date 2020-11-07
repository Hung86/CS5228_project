# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);

# Load the WGCNA package
library(WGCNA)

# Use gplots for heatmap color panel (fix color mapping)
library(gplots)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()


# Load the expression and trait data saved in the first part
lnames = load(file = "./Routput/1-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "./Routput/2-networkConstruction.RData");
lnames
nGenes = ncol(gbmExpr)
nSamples = nrow(gbmExpr)


#################### GBM ##########################

png(file = "./results/3_Network_heatmap_GBM.png",width=10,height=10,units="in",res=1200);
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = gbmDissTOM^7;
diag(plotDiss) = NA;
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, gbmGeneTree, gbmModuleColors, main = "Network Heatmap Plot - GBM", col=myheatcol)
dev.off()


###Visualizing the network of eigengenes######
MEs = moduleEigengenes(gbmExpr, gbmModuleColors)$eigengenes
n_MEs = ncol(MEs)
sample_MEs = nrow(MEs)

# Plot the dendrogram
png(file = "./results/3_Eigengene_dendrogram_GBM.png",width=10,height=10,units="in",res=1200);
par(cex = 1.0)
plotEigengeneNetworks(MEs , "Eigengene Dendrogram - GBM", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix
png(file = "./results/3_Eigengene_heatmap_GBM.png",width=10,height=10,units="in",res=1200);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene Adjacency Heatmap - GBM", marDendro = c(0,4,2,0), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

##################### OV #######################

png(file = "./results/3_Network_heatmap_OV.png",width=10,height=10,units="in",res=1200);
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = ovDissTOM^7;
diag(plotDiss) = NA;
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, ovGeneTree, ovModuleColors, main = "Network Heatmap Plot - OV", col=myheatcol)
dev.off()

###Visualizing the network of eigengenes######
MEs = moduleEigengenes(ovExpr, ovModuleColors)$eigengenes
n_MEs = ncol(MEs)
sample_MEs = nrow(MEs)

# Plot the dendrogram
png(file = "./results/3_Eigengene_dendrogram_OV.png",width=10,height=10,units="in",res=1200);
par(cex = 1.0)
plotEigengeneNetworks(MEs , "Eigengene Dendrogram - OV", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix
png(file = "./results/3_Eigengene_heatmap_OV.png",width=10,height=10,units="in",res=1200);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene Adjacency Heatmap - OV", marDendro = c(0,4,2,0), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


##################### BRCA #######################

png(file = "./results/3_Network_heatmap_BRCA.png",width=10,height=10,units="in",res=1200);
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = brcaDissTOM^7;
diag(plotDiss) = NA;
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, brcaGeneTree, brcaModuleColors, main = "Network Heatmap Plot - BRCA", col=myheatcol)
dev.off()

###Visualizing the network of eigengenes######
MEs = moduleEigengenes(brcaExpr, brcaModuleColors)$eigengenes
n_MEs = ncol(MEs)
sample_MEs = nrow(MEs)

# Plot the dendrogram
png(file = "./results/3_Eigengene_dendrogram_BRCA.png",width=10,height=10,units="in",res=1200);
par(cex = 1.0)
plotEigengeneNetworks(MEs , "Eigengene Dendrogram - BRCA", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix
png(file = "./results/3_Eigengene_heatmap_BRCA.png",width=10,height=10,units="in",res=1200);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene Adjacency Heatmap - BRCA", marDendro = c(0,4,2,0), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

