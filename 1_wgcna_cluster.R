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

tcgaData = read.csv("GBM_HiSeqV2_PANCAN",sep = "\t",header = T,row.names=1)

# Take a quick look at what is in the data set:
dim(tcgaData);
names(tcgaData);

x<- as.matrix(tcgaData)
filtered_data <- genefilter::varFilter(x, var.func = IQR, 
            var.cutoff = 0.75, filterByQuantile = TRUE)
dim(filtered_data)
dim(tcgaData)
#######filtering##################
#f1 <- pOverA(0.25, 2)
#f2 <- function(x) (IQR(x) > 0.5)
#ff <- filterfun(f1, f2)
#filtered_data <- tcgaData[genefilter(tcgaData, ff), ]
#dim(filtered_data)

f1 <- pOverA(0.8, 4)
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
interested_genes<- tcgaData[genefilter(tcgaData, ff), ]
dim(interested_genes)
###################################
datExpr0 = as.data.frame(t(filtered_data[,]))
dim(datExpr0)
#datExpr0 = as.data.frame(t(tcgaData[, -c(1)]))
#names(datExpr0) = tcgaData$sample;

#rownames(datExpr0) = tcgaData$Subject
#names(datExpr0) = tcgaData$Subject

###################################
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
###################################
if (!gsg$allOK)
{
	# Optionally, print the gene and sample names that were removed:
	if (sum(!gsg$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
	if (sum(!gsg$goodSamples)>0)
	printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
	# Remove the offending genes and samples from the data:
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK





#datExpr0_T$variance = apply(datExpr0_T, 1, var)
#datExpr  = datExpr0_T[datExpr0_T$variance >= quantile(datExpr0_T$variance, c(.90)), ]
#datExpr$variance <- NULL

datExpr  = datExpr0
col = ncol(datExpr)
row = nrow(datExpr)

hierarchicalTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
#pdf(file = "Patient_Clustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(hierarchicalTree , main = "Patient clustering to detect outliers", 
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)



save(datExpr, file = "GBM_1-dataInput.RData")

collectGarbage()

