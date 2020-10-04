# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the data set
tcgaData = read.csv("tcga_wgcna_counts.csv");
# Take a quick look at what is in the data set:
dim(tcgaData);
names(tcgaData);


###################################
datExpr0 = as.data.frame(tcgaData[, -c(1)])
#datExpr0 = as.data.frame(t(tcgaData[, -c(1)]))

rownames(datExpr0) = tcgaData$Subject
#names(datExpr0) = tcgaData$Subject

###################################
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
###################################

hierarchicalTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9)
#pdf(file = "cancerClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(hierarchicalTree , main = "Cancer clustering to detect outliers", 
		sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)



datExpr =datExpr0

save(datExpr , file = "1-dataInput.RData")

collectGarbage()
