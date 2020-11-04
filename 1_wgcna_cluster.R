# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
library(genefilter)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the data set
#tcgaData = read.delim("GBM_HiSeqV2_PANCAN");

rawGBM= read.csv("GBM_HiSeqV2_PANCAN",sep = "\t",header = T,row.names=1)
rawOV= read.csv("OV_HiSeqV2_PANCAN",sep = "\t",header = T,row.names=1)
rawBRCA= read.csv("BRCA_HiSeqV2_PANCAN",sep = "\t",header = T,row.names=1)


# Take a quick look at what is in the data set:
dim(rawGBM);
dim(rawOV);
dim(rawBRCA);


x<- as.matrix(rawGBM)
filteredGBM<- genefilter::varFilter(x, var.func = IQR, 
            var.cutoff = 0.75, filterByQuantile = TRUE)
dim(filteredGBM)

y<- as.matrix(rawOV)
filteredOV<- genefilter::varFilter(y, var.func = IQR, 
            var.cutoff = 0.75, filterByQuantile = TRUE)
dim(filteredOV)

z<- as.matrix(rawBRCA)
filteredBRCA<- genefilter::varFilter(z, var.func = IQR, 
            var.cutoff = 0.75, filterByQuantile = TRUE)
dim(filteredBRCA)
#######filtering##################
#f1 <- pOverA(0.25, 2)
#f2 <- function(x) (IQR(x) > 0.5)
#ff <- filterfun(f1, f2)
#filtered_data <- tcgaData[genefilter(tcgaData, ff), ]
#dim(filtered_data)

#f1 <- pOverA(0.8, 4)
#f2 <- function(x) (IQR(x) > 0.5)
#ff <- filterfun(f1, f2)
#interested_genes<- tcgaData[genefilter(tcgaData, ff), ]
#dim(interested_genes)
###################################
gbmExpr = as.data.frame(t(filteredGBM[,]))
dim(gbmExpr)

ovExpr = as.data.frame(t(filteredOV[,]))
dim(ovExpr)

brcaExpr = as.data.frame(t(filteredBRCA[,]))
dim(brcaExpr)

#datExpr0 = as.data.frame(t(tcgaData[, -c(1)]))
#names(datExpr0) = tcgaData$sample;

#rownames(datExpr0) = tcgaData$Subject
#names(datExpr0) = tcgaData$Subject

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
sizeGrWindow(12,9)
#pdf(file = "Patient_Clustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(gbmPatientTree , main = "GBM patient clustering", 
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)


ovPatientTree = hclust(dist(ovExpr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(ovPatientTree , main = "OV patient clustering", 
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)



brcaPatientTree = hclust(dist(brcaExpr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(brcaPatientTree , main = "BRCA patient clustering", 
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)



save(gbmExpr, ovExpr, brcaExpr, file = "1-dataInput.RData")

collectGarbage()

