# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
# enableWGCNAThreads()

allowWGCNAThreads()

# Load the data saved in the first part
lnames = load(file = "./Routput/1-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
dim(gbmExpr)
dim(ovExpr)
dim(brcaExpr)

#################gbmExpr######################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(gbmExpr, powerVector = powers, verbose = 5)
png(file = "./results/2_SFTM_Fit_GBM.png", width = 600, height = 700);
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit (signed R^2)",type="n",
     main = paste("Scale Independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(file = "./results/2_Mean_Connectivity_GBM.png", width = 600, height = 700);
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#################ovExpr######################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(ovExpr, powerVector = powers, verbose = 5)
png(file = "./results/2_SFTM_Fit_OV.png", width = 600, height = 700);
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit (signed R^2)",type="n",
     main = paste("Scale Independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(file = "./results/2_Mean_Connectivity_OV.png", width = 600, height = 700);
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#################brcaExpr######################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(brcaExpr, powerVector = powers, verbose = 5)
png(file = "./results/2_SFTM_Fit_BRCA.png", width = 600, height = 700);
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit (signed R^2)",type="n",
     main = paste("Scale Independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(file = "./results/2_Mean_Connectivity_BRCA.png", width = 600, height = 700);
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



############# GBM ##################################################################################################
##Co-expression similarity and adjacency
softPower = 9;
adjacency = adjacency(gbmExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
gbmDissTOM = 1-TOM


# Call the hierarchical clustering function
gbmGeneTree = hclust(as.dist(gbmDissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = gbmGeneTree, distM = gbmDissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#plotDendroAndColors(gbmGeneTree, dynamicColors, "Dynamic Tree Cut",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Gene Dendrogram and Module Colors - GBM")

########################################
MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(gbmExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
gbmModuleColors = merge$colors;
# Eigengenes of the new merged modules:
gbmMEs = merge$newMEs;

png(file = "./results/2_Clustering_Tree_GBM.png",width=10,height=10,units="in",res=1200);
plotDendroAndColors(gbmGeneTree, cbind(dynamicColors, gbmModuleColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    main = "Gene Dendrogram and Module Colors - GBM",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

############# OV ##################################################################################################
##Co-expression similarity and adjacency
softPower = 6;
adjacency = adjacency(ovExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
ovDissTOM = 1-TOM

# Call the hierarchical clustering function
ovGeneTree = hclust(as.dist(ovDissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = ovGeneTree, distM = ovDissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#####################################3
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#plotDendroAndColors(ovGeneTree, dynamicColors, "Dynamic Tree Cut",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Gene dendrogram and module colors")

########################################
MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(ovExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
ovModuleColors = merge$colors;
# Eigengenes of the new merged modules:
ovMEs = merge$newMEs;


png(file = "./results/2_Clustering_Tree_OV.png",width=10,height=10,units="in",res=1200);
plotDendroAndColors(ovGeneTree, cbind(dynamicColors, ovModuleColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    main = "Gene Dendrogram and Module Colors - OV",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


############# BRCA ##################################################################################################
##Co-expression similarity and adjacency
softPower = 12;
adjacency = adjacency(brcaExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
brcaDissTOM = 1-TOM

# Call the hierarchical clustering function
brcaGeneTree = hclust(as.dist(brcaDissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = brcaGeneTree, distM = brcaDissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(brcaExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
brcaModuleColors = merge$colors;
# Eigengenes of the new merged modules:
brcaMEs = merge$newMEs;

png(file = "./results/2_Clustering_Tree_BRCA.png",width=10,height=10,units="in",res=1200);
plotDendroAndColors(brcaGeneTree, cbind(dynamicColors, brcaModuleColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    main = "Gene Dendrogram and Module Colors - BRCA",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


###########################################################################################################

# png(file = "./results/2_Clustering_Tree_ALL.png",width=10,height=10,units="in",res=1200);
# plotDendroAndColors(cbind(gbmModuleColors, ovModuleColors, brcaModuleColors),
#                     c("GBM", "OV", "BRCA"),
#                     main = "Gene Dendrogram and Module Colors - All",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
gbmModuleLabels = match(gbmModuleColors, colorOrder)-1;
ovModuleLabels = match(ovModuleColors, colorOrder)-1;
brcaModuleLabels = match(brcaModuleColors, colorOrder)-1;

# Save module colors and labels for use in subsequent parts
save(gbmMEs, gbmModuleLabels, gbmModuleColors, gbmGeneTree, gbmDissTOM,
     ovMEs, ovModuleLabels, ovModuleColors, ovGeneTree, ovDissTOM,
     brcaMEs, brcaModuleLabels, brcaModuleColors, brcaGeneTree, brcaDissTOM,
     file = "./Routput/2-networkConstruction.RData")
