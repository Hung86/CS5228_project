####################Filtering Functions################################################

filteringGeneDataset <- function(raw_gene_csv, variance)
{
	removedGenes = startsWith(rownames(raw_gene_csv) , "?")
	keepGenes = raw_gene_csv[!removedGenes,]
	x<- as.matrix(keepGenes)
	remainingGenes<- genefilter::varFilter(x, var.func = IQR, 
            	var.cutoff = variance, filterByQuantile = TRUE)
	return (remainingGenes)
}

####################Gene Enrichment Functions################################################
geneModules_ToFile <- function(prefix, MEs)
{

	moduleLabels = substring(names(MEs), 3)
	print(typeof(moduleLabels ))
	fileName = paste(prefix, "_GeneModules.txt", sep="")
	if (file.exists(fileName)) {
		file.remove(fileName)
	}
	lapply(moduleLabels , write, fileName, append=TRUE, ncolumns=100)
}


findGeneOverlap <- function(dataset1, dataset2, moduleColors1, moduleColors2, MEs1, MEs2)
{

	modNames1 = substring(names(MEs1), 3)
	modNames2 = substring(names(MEs2), 3)
	intersectionGeneList <- list()
	unionListGeneList <- list()
	positionGeneList <- list()
	modulePairList <- list()
	idx = 1
	for (m1 in modNames1 )
	{
		moduleGene1 <- dataset1[moduleColors1==m1]
		geneNames1 = names(moduleGene1)
		for (m2 in modNames2)
		{
			moduleGene2 <- dataset2[moduleColors2==m2]
			geneNames2 = names(moduleGene2)
			go.obj <- newGeneOverlap(geneNames1 , geneNames2)
			go.obj <- testGeneOverlap(go.obj)
			
			if (getContbl(go.obj)[4][1] > 0)
			{	
					
					modulePairList[idx] <- list(c(m1,m2)) 
					positionGeneList[idx] <- list(c(getPval(go.obj), idx))
					intersectionGeneList[idx]  <- list(getIntersection(go.obj))
					unionListGeneList[idx] <- list(getUnion(go.obj))
					idx = idx+1
			}
		}
	}
	return (list(intersectionList=intersectionGeneList,
			 unionList=unionListGeneList,
			 positionList=positionGeneList,
			 modulePairList=modulePairList))
}



geneOverlap_ToFile <- function(prefix, geneOverlapOjb, high_pvalue = TRUE){
	positionGeneList = geneOverlapOjb$positionList
	modulePairList = geneOverlapOjb$modulePairList
	unionListGeneList = geneOverlapOjb$unionList
	intersectionGeneList = geneOverlapOjb$intersectionList
	geneOverlapOrder = NULL
	
	if (high_pvalue) {
		geneOverlapOrder = positionGeneList[order(sapply(positionGeneList, function(x) x[1], simplify=TRUE), decreasing=TRUE)][0:3]
	} else {
		geneOverlapOrder = positionGeneList[order(sapply(positionGeneList, function(x) x[1], simplify=TRUE), decreasing=FALSE)][0:3]
	}
	
	for (v in geneOverlapOrder)
	{
		pValue = v[[1]]
		idx = v[[2]]
		module_labels = modulePairList[[idx]] 
		candidate_list = intersectionGeneList[[idx]]
		fileName = paste(prefix, "_Intersection_P",pValue, "_" , module_labels[[1]], "_", module_labels[[2]], "_genes.txt", sep="")
		if (file.exists(fileName)) {
			file.remove(fileName)
		}
		lapply(candidate_list, write, fileName, append=TRUE, ncolumns=30000)
	}
	
	for (v in geneOverlapOrder)
	{
		pValue = v[[1]]
		idx = v[[2]]
		module_labels = modulePairList[[idx]] 
		candidate_list = unionListGeneList[[idx]]
		fileName = paste(prefix, "_Union_P",pValue, "_" , module_labels[[1]], "_", module_labels[[2]], "_genes.txt", sep="")
		if (file.exists(fileName)) {
			file.remove(fileName)
		}
		lapply(candidate_list, write, fileName, append=TRUE, ncolumns=30000)
	}

}



topGO_GeneEnrichment <- function (bg_genes, geneOverlapOjb, high_pvalue = TRUE)
{
	positionGeneList = geneOverlapOjb$positionList
	modulePairList = geneOverlapOjb$modulePairList
	unionListGeneList = geneOverlapOjb$unionList
	intersectionGeneList = geneOverlapOjb$intersectionList
	
	intersection_fishers = list()
	intersection_enrichments = list()
	intersection_GOdata= list()
	intersection_pairs = list()

	union_fishers = list()
	union_enrichments = list()
	union_GOdata= list()
	union_pairs = list()
	
	geneOverlapOrder = NULL
	if (high_pvalue) {
		geneOverlapOrder = positionGeneList[order(sapply(positionGeneList, function(x) x[1], simplify=TRUE), decreasing=TRUE)][0:3]
	} else {
		geneOverlapOrder = positionGeneList[order(sapply(positionGeneList, function(x) x[1], simplify=TRUE), decreasing=FALSE)][0:3]
	}

	db = useMart('ENSEMBL_MART_ENSEMBL',
			dataset='hsapiens_gene_ensembl', 
			host="uswest.ensembl.org")

	go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), 
						filters='external_gene_name', 
						values=bg_genes, mart=db)
					
	gene_2_GO=unstack(go_ids[,c(1,2)])	
	i = 1
	for (v in geneOverlapOrder)
	{
		idx = v[[2]]
		module_labels = modulePairList[[idx]] 
		candidate_list = intersectionGeneList[[idx]]
		# remove any candidate genes without GO annotation
		keep = candidate_list %in% go_ids[,2]
		keep = which(keep==TRUE)
		candidate_list=candidate_list[keep]
		geneList=factor(as.integer(bg_genes %in% candidate_list))
		names(geneList)= bg_genes
		GOdata=new('topGOdata', ontology='BP',
					allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
		weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher') 
		# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
		allGO=usedGO(GOdata)
		all_res=GenTable(GOdata, weightFisher = weight_fisher_result,
						orderBy='weightFisher', topNodes=length(allGO))
		
		#performing BH correction on our p values
		p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
		 
		# create the file with all the statistics from GO analysis
		all_res_final=cbind(all_res,p.adj)
		all_res_final=all_res_final[order(all_res_final$p.adj),]
		
		intersection_fishers[i] = list(weight_fisher_result)
		intersection_enrichments[i] = list(all_res_final)
		intersection_GOdata[i]= list(GOdata)
		intersection_pairs[i] = list(module_labels)
		i = i + 1
	}

	i = 1
	for (v in geneOverlapOrder)
	{
		idx = v[[2]]
		module_labels = modulePairList[[idx]] 
		candidate_list = unionListGeneList[[idx]]
		# remove any candidate genes without GO annotation
		keep = candidate_list %in% go_ids[,2]
		keep = which(keep==TRUE)
		candidate_list=candidate_list[keep]
		geneList=factor(as.integer(bg_genes %in% candidate_list))
		names(geneList)= bg_genes
		GOdata=new('topGOdata', ontology='BP',
					allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
		weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher') 
		# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
		allGO=usedGO(GOdata)
		all_res=GenTable(GOdata, weightFisher = weight_fisher_result,
						orderBy='weightFisher', topNodes=length(allGO))
		
		#performing BH correction on our p values
		p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
		 
		# create the file with all the statistics from GO analysis
		all_res_final=cbind(all_res,p.adj)
		all_res_final=all_res_final[order(all_res_final$p.adj),]
		
		union_fishers[i] = list(weight_fisher_result)
		union_enrichments[i] = list(all_res_final)
		union_GOdata[i]= list(GOdata)
		union_pairs[i] = list(module_labels)
		i = i + 1
	}
	
	return (list(intersection_fishers = intersection_fishers,
				intersection_enrichments = intersection_enrichments,
				intersection_GOdata = intersection_GOdata,
				intersection_pairs = intersection_pairs,
				union_fishers = union_fishers,
				union_enrichments = union_enrichments,
				union_GOdata = union_GOdata,
				union_pairs = union_pairs))
}

topGO_GeneEnrichment_ToFile <- function(prefix, topGOEnrichment)
{
	intersection_fishers = topGOEnrichment$intersection_fishers
	intersection_enrichments = topGOEnrichment$intersection_enrichments
	intersection_GOdata = topGOEnrichment$intersection_GOdata
	intersection_pairs = topGOEnrichment$intersection_pairs

	union_fishers = topGOEnrichment$union_fishers
	union_enrichments = topGOEnrichment$union_enrichments
	union_GOdata = topGOEnrichment$union_GOdata
	union_pairs = topGOEnrichment$union_pairs
	
	len = length(intersection_enrichments)

	for (i in 1:len)
	{
		fisher = intersection_fishers[[i]]
		res = intersection_enrichments[[i]]
		GOdata = intersection_GOdata[[i]]
		pair = intersection_pairs[[i]]
		
		#save first top 50 ontolgies sorted by adjusted pvalues
		write.table(res[1:50,],paste(prefix, "_Intersection_", pair[[1]], "_", pair[[2]], "_summary_topGO_analysis.csv", sep=""), quote=FALSE,row.names=FALSE)

		# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
		pdf(paste(prefix, "_Intersection_", pair[[1]], "_", pair[[2]], "_topGOPlot_fullnames.pdf", sep=""),
			height=12, width=12, paper='special', pointsize=18)
		showSigOfNodes(GOdata, score(fisher), useInfo = "all", sigForAll=FALSE, firstSigNodes=2, .NO.CHAR=50)
		dev.off()
	}
	
	for (i in 1:len)
	{
		fisher = union_fishers[[i]]
		res = union_enrichments[[i]]
		GOdata = union_GOdata[[i]]
		pair = union_pairs[[i]]
		
		#save first top 50 ontolgies sorted by adjusted pvalues
		write.table(res[1:50,],paste(prefix, "_Union_", pair[[1]], "_", pair[[2]], "_summary_topGO_analysis.csv", sep=""), quote=FALSE,row.names=FALSE)

		# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
		pdf(paste(prefix, "_Union_", pair[[1]], "_", pair[[2]], "_topGOPlot_fullnames.pdf", sep=""),
			height=12, width=12, paper='special', pointsize=18)
		showSigOfNodes(GOdata, score(fisher), useInfo = "all", sigForAll=FALSE, firstSigNodes=2, .NO.CHAR=50)
		dev.off()
	}
}


