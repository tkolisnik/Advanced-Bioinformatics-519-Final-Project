#pathwaysFromReads.R

#Created by Tyler Kolisnik
#With code taken from scripts produced by Tyler Kolisnik, Mark Bieda and Nathan Cormier.

#This is a program that takes in an annotated list file (i.e. from annotatePeaks.R) that has 
#raw sequence read counts in the 5th column (note: ensure not ENRICHMENT values)
#If using the Cormier et al. pipeline raw sequence read "tags" can be annotated to the 5th column by
#editing the peaks.xls output file to replace the enrichment column with the tags column, changing the format to .bed and running annotatePeaks.R

########################### Begin Parameter Block ###########################
inputFile <- ("/home/tylerk/MDSC519/mcf7/tagcount_mcf7_5000up_2000down_allData.txt") # Must be a text file with entrezIDs and log2 Fold Change values. 
outputDirectory <- ("/home/tylerk/MDSC519/mcf7/pathwayResults/readResultsApril25") # The absolute path files will be outputted to, with no trailing /
runName <- "RawReads" # The name that will be appended to all output files. 
pValueCutoff <- 0.2 # The pValue threshold for the Fold Change. 
KEGGspeciesCode <- "hsa" # Kegg species code, default = "hsa" (human), more available on bioconductor website. 

########################### End Parameter Block ###########################
library(gage)
library(pathview)
########################################## PathwayAnalysis Block ########################################## 

setwd(outputDirectory)

## gets the list of entrez genes and their associated peak scores from the input file
geneData <- unique(read.delim(inputFile)[,c("entrezID", "score")])
geneList <- unique(geneData$entrezID)
geneList <- geneList[!is.na(geneList)]
dataMatrix <- matrix(data=NA, nrow=length(geneList), ncol=1) 
rownames(dataMatrix) <- geneList
colnames(dataMatrix) <- "score"

## for each unique entrez gene, finds the raw sequence read number associated with that gene
for (i in 1:length(geneList)) {
  gene <- as.character(geneList[i])
  #print(gene)
  score <- max(geneData$score[which(geneData$entrezID == gene)])
  dataMatrix[gene,] <- score
}

data(kegg.gs) # load in the KEGG pathway data
gageOutput <- gage(log2(dataMatrix), gsets = kegg.gs, samp=NULL) # determines which pathways are enriched in the data set
pathwayData <- gageOutput$greater[,c("p.val", "q.val")]
pathCode <- substr(rownames(pathwayData), 1, 8)
pathName <- substr(rownames(pathwayData), 9, length(rownames(pathwayData)))
pathwayData <- cbind(pathCode, pathName ,pathwayData)
colnames(pathwayData) <- c("Pathway Code", "Pathway Name", "p Values", "q value")
#write.table(pathwayData, paste(runName, "full_pathway_list_new_Reads", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

write.table(gageOutput$greater,file="gageOutputGreater_Enrichment.xls", sep="\t")
## Use pathview to create output images

runPathview <- function(pid, upordown){
  # Output Kegg Native Graphs (.png files, works for all Pathways)
  pathview(gene.data=log2(dataMatrix), pathway.id=pid, kegg.native=T, species=KEGGspeciesCode, limit=(15), low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20),plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="mean", out.suffix=paste(runName, "KEGGnative_",upordown,sep=""))
  # Output Graphviz Graphs (.pdf files, doesn't work for all Pathways)
  pathview(gene.data=log2(dataMatrix), pathway.id=pid, kegg.native=F, species=KEGGspeciesCode, limit=(15), low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20), plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="mean", out.suffix=paste(runName, "graphviz_",upordown,sep=""))
}

# Select Up Regulated Pathways to be printed based on p value (p.val)
selUp <- gageOutput$greater[,"p.val"] < pValueCutoff & !is.na(gageOutput$greater[,"p.val"])
pathIdsUp <- substr(rownames(gageOutput$greater)[selUp], 1, 8)

# Print Selected Up Regulated Pathways 
# Try Catch allows for error handling if a Graphviz visualization does not exist for the pathway ID entered.
pathListUp <- sapply(pathIdsUp, function(pid) tryCatch(runPathview(pid, "upRegulated"), error=function(e) print("This pathway not found")))
########################################## End PathwayAnalysis Block ########################################## 
