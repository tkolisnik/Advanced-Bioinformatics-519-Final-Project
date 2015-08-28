#pathwaysFromBinaryData.R

#Created by Tyler Kolisnik
#With code taken from scripts produced by Tyler Kolisnik, Mark Bieda and Nathan Cormier.

#DESCRIPTION:
#This program takes in an annotated mapped ChIP-seq reads list text file, and assigns all of the genes with a peak a score of '1',
# it then obtains a complete gene list for the species, merges the two and assigns all the ones that 
# do not have a score of '1' a score of '0' and it then is inputted into gage for binary pathway analysis. 


########################### Begin Parameter Block ###########################
inputFile <- ("/home/tylerk/MDSC519/mcf7/mcf7_enrichment_up_2000down_allData.txt") # Must be a text file with entrezIDs and log2 Fold Change values. 
outputDirectory <- ("/home/tylerk/MDSC519/mcf7/pathwayResults") # The absolute path files will be outputted to, with no trailing /
runName <- "April25v2" # The name that will be appended to all output files. 
pValueCutoff <- 0.2 # The pValue threshold for the Fold Change. 
KEGGspeciesCode <- "hsa" # Kegg species code, default = "hsa" (human), more available on bioconductor website. 
pathwayDirectoryName <- "pathwaysApril25_binary"
########################### End Parameter Block ###########################

library(gage)
library(org.Hs.eg.db) #Must be changed to the db for the correct species if not using human
library(pathview)
########################################## gagePathwayAnalysis ########################################## 

setwd(outputDirectory)
pathwayDirectory <- paste(outputDirectory, pathwayDirectoryName, sep="/")
## gets the list of entrez genes and their associated peak scores from the input file
geneData <- unique(read.delim(inputFile)[,c("entrezID", "score")])
geneList <- unique(geneData$entrezID)
geneList <- geneList[!is.na(geneList)]
dataMatrix <- matrix(data=NA, nrow=length(geneList), ncol=1) 
rownames(dataMatrix) <- geneList
colnames(dataMatrix) <- "score"

## for each unique entrez gene with a peak, assign it a score of 1
for (i in 1:length(geneList)) {
  gene <- as.character(geneList[i])
  #print(gene)
  score <- 1
  dataMatrix[gene,] <- score
}
dataFrame<-as.data.frame(dataMatrix)

#data structure manipulation, messy but necessary, as merge is only supported in data frames, and 
#gage requires a specially formatted data matrix. 

allKeys<- org.Hs.egACCNUM #Must be changed if not using human
allKeys2<-mappedkeys(allKeys) #obtain all keys from the org db.
keyMatrix<-matrix(data=0,nrow=length(allKeys2),ncol=1)
rownames(keyMatrix)<-allKeys2
keyFrame <- as.data.frame(keyMatrix)
colnames(keyFrame) <- "score"
mergedDF<-merge(keyFrame,dataFrame,by="row.names",all.x="TRUE") #merges the data frames
mergedDF <- mergedDF[-2] #removes an unwanted byproduct column
mergedDF[is.na(mergedDF)]<-0 #assign all NA's a score of 0
rownames(mergedDF)<-mergedDF[,1] #change rownames
mergedDF <- mergedDF[-1] #remove unwanted byproduct column
mergedDM<- data.matrix(mergedDF) #convert it back to a data matrix (required by gage)

data(kegg.gs) # load in the KEGG pathway data
gageOutput <- gage(mergedDM, gsets = kegg.gs, samp=NULL) # determines which pathways are enriched in the data set
pathwayData <- gageOutput$greater[,c("p.val", "q.val")]
pathCode <- substr(rownames(pathwayData), 1, 8)
pathName <- substr(rownames(pathwayData), 9, length(rownames(pathwayData)))
pathwayData <- cbind(pathCode, pathName ,pathwayData)
colnames(pathwayData) <- c("Pathway Code", "Pathway Name", "p Values", "q value")
write.table(gageOutput$greater,file="gageOutputGreater_fromBinary.xls", sep="\t")

## Selects the pathways to be printed (the ones that are above the p value cut off) 
sel <- gageOutput$greater[,"p.val"] < pValueCutoff & !is.na(gageOutput$greater[,"p.val"])
pathIds <- substr(rownames(gageOutput$greater)[sel], 1, 8)

system(paste("mkdir ", pathwayDirectory, sep = ""))
setwd(pathwayDirectory)

runPathview <- function(pid, upordown){
  # Output Kegg Native Graphs (.png files, works for all Pathways)
  pathview(gene.data=mergedDM, pathway.id=pid, kegg.native=T, species=KEGGspeciesCode,  low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20),plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="mean", out.suffix=paste(runName, "KEGGnative_",upordown,sep=""))
  # Output Graphviz Graphs (.pdf files, doesn't work for all Pathways)
  pathview(gene.data=mergedDM, pathway.id=pid, kegg.native=F, species=KEGGspeciesCode, low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20), plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="mean", out.suffix=paste(runName, "graphviz_",upordown,sep=""))
}

# Print Selected Up Regulated Pathways 
# Try Catch allows for error handling if a Graphviz visualization does not exist for the pathway ID entered.
pathListUp <- sapply(pathIds, function(pid) tryCatch(runPathview(pid, "upRegulated"), error=function(e) print("This pathway not found")))


