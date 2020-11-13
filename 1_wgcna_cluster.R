# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
library(genefilter)
source("utility_functions.R")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

#Read in the data set

rawGBM= read.csv("./data/GBM_HiSeqV2_PANCAN", sep = "\t", header = T, row.names=1)
rawOV= read.csv("./data/OV_HiSeqV2_PANCAN", sep = "\t", header = T, row.names=1)
rawBRCA= read.csv("./data/BRCA_HiSeqV2_PANCAN", sep = "\t", header = T, row.names=1)


# Take a quick look at what is in the data set:
dim(rawGBM);
dim(rawOV);
dim(rawBRCA);

filteredGBM = filteringGeneDataset(rawGBM, 0.5)
filteredOV = filteringGeneDataset(rawOV, 0.5)
filteredBRCA = filteringGeneDataset(rawBRCA, 0.5)

dim(filteredGBM)
dim(filteredOV)
dim(filteredBRCA)


gbmExpr = as.data.frame(t(filteredGBM[,]))
dim(gbmExpr)

ovExpr = as.data.frame(t(filteredOV[,]))
dim(ovExpr)

brcaExpr = as.data.frame(t(filteredBRCA[,]))
dim(brcaExpr)

###################################
gbm = goodSamplesGenes(gbmExpr , verbose = 3);
gbm$allOK
###################################
if (!gbm$allOK)
{
	# Optionally, print the gene and sample names that were removed:
	if (sum(!gbm$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(gbmExpr)[!gbm$goodGenes], collapse = ", ")));
	if (sum(!gbm$goodSamples)>0)
	printFlush(paste("Removing samples:", paste(rownames(gbmExpr)[!gbm$goodSamples], collapse = ", ")));
	# Remove the offending genes and samples from the data:
	gbmExpr = gbmExpr [gbm$goodSamples, gbm$goodGenes]
}

###################################
ov = goodSamplesGenes(ovExpr , verbose = 3);
ov$allOK
###################################
if (!ov$allOK)
{
	# Optionally, print the gene and sample names that were removed:
	if (sum(!ov$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(ovExpr)[!ov$goodGenes], collapse = ", ")));
	if (sum(!ov$goodSamples)>0)
	printFlush(paste("Removing samples:", paste(rownames(ovExpr)[!ov$goodSamples], collapse = ", ")));
	# Remove the offending genes and samples from the data:
	ovExpr = ovExpr [ov$goodSamples, ov$goodGenes]
}

###################################
brca = goodSamplesGenes(brcaExpr , verbose = 3);
brca$allOK
###################################
if (!brca$allOK)
{
	# Optionally, print the gene and sample names that were rembrcaed:
	if (sum(!brca$goodGenes)>0)
	printFlush(paste("Rembrcaing genes:", paste(names(brcaExpr)[!brca$goodGenes], collapse = ", ")));
	if (sum(!brca$goodSamples)>0)
	printFlush(paste("Rembrcaing samples:", paste(rownames(brcaExpr)[!brca$goodSamples], collapse = ", ")));
	# Rembrcae the offending genes and samples from the data:
	brcaExpr = brcaExpr [brca$goodSamples, brca$goodGenes]
}


gbmPatientTree = hclust(dist(gbmExpr), method = "average")
png(file = "./results/1_Sample_Clustering_GBM.png", width = 1400, height = 1500);
par(cex = 1.0);
par(mar = c(1,1,1,1))
plot(gbmPatientTree , main = "Sample Clustering - GBM",
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

ovPatientTree = hclust(dist(ovExpr), method = "average")
sizeGrWindow(14,10)
png(file = "./results/1_Sample_Clustering_OV.png", width = 1400, height = 1500);
par(cex = 1.0);
par(mar = c(1,1,1,1))
plot(ovPatientTree , main = "Sample Clustering - OV",
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()


brcaPatientTree = hclust(dist(brcaExpr), method = "average")
sizeGrWindow(14,10)
png(file = "./results/1_Sample_Clustering_BRCA.png", width = 1400, height = 1500);
par(cex = 1.0);
par(mar = c(1,1,1,1))
plot(brcaPatientTree , main = "Sample Clustering - BRCA",
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()


save(gbmExpr, ovExpr, brcaExpr, file = "./Routput/1_dataInput.RData")

collectGarbage()

