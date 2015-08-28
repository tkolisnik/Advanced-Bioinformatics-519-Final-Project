## pathwayListFromGeneExpression.R

## With code taken from scripts produced by Tyler Kolisnik, Mark Bieda and Nathan Cormier.
## Created By: Tyler Kolisnik and Mark Bieda.

## Description: This is a program to analyze gene expression data from CEL files, outputted from a DNA or RNA microarray.
# It normalizes array data, generates differentially expressed genes, and produces a gage pathway list. 

## Run Time: Dependent on size of input files, but usually no more than 5 minutes. 

## Dependencies: 
library(simpleaffy) # Required for data normalization.
library(gage) # Required for Pathway Analysis.

##################################################################### Begin Parameter Block #####################################################################

inputDirectory <- "/home/tylerk/MDSC519/mcf7/GE" # Path of directory where input is located, must contain .CEL files and a covdesc file, no trailing /.
outputDirectory <- "/home/tylerk/MDSC519/mcf7/GE/output" # Path of directory where output will be located, no trailing /.
runName <- "April25GE" # A string that will be appended to all output files.
minFC <- 1.5 # This is the minimum log2 fold change for a gene to be differentially expressed.
ttestPVal <- 0.01 # This is the threshold p value for significance for differentially expressed genes.
hgCutoff <- 0.01 # This is the GOStats p value threshold.
pathwayPvalueCutoff <- 0.01 # Numeric PValue Threshold for Pathway analysis
KEGGspeciesCode <- "hsa" # A string of the KEGG code for the species you are using, default (human) is "hsa".
exonArrayPackageName <- "pd.huex.1.0.st.v2" # The package name for the exon array, only required if it is an exon array.
dbpackagename <- "huex10sttranscriptcluster.db" # The annotation database package, find the correct one for your array on the bioconductor website.
isItAnExonArray <- "yes" # String specifying if it is an exon (RNA) array, if it is not, set to "no" and it will assume it is a DNA array.
pathwayDirectoryName <- "pathway_imagesApril26" # String depicting the name of the folder within the output directory that will be created to contain the pathway images. 
orgDB <- "org.Hs.eg.db" # Organism Database ID (used by GOstats), must be changed if using an organism other than human (default).
##################################################################### End Parameter Block #######################################################################

##################################################################### Begin Coding Block #####################################################################
library(orgDB, character.only=TRUE)
library(dbpackagename, character.only=TRUE)
library(exonArrayPackageName, character.only=TRUE)
assign("dbpackagename", get(dbpackagename))
########################### Begin Normalization Block #########################

options(max.print=100)
setwd(inputDirectory)

if(isItAnExonArray == "yes"){
rawData <- read.celfiles(list.celfiles(),pkgname=exonArrayPackageName)
dataSet <- rma(rawData, target="core")
featureData(dataSet) <- getNetAffx(dataSet,"transcript")
column_adjust = read.table("covdesc", sep="\t", col.names="conditionName")
pData(dataSet) <- column_adjust

} else {
dataSet <- justRMA(phenoData="covdesc") # The input data is normalized using the justRMA function from the simpleaffy package here.
}

pdSet <- pData(dataSet)
pdSet[,"blankcol"] <- "c"
pData(dataSet) <- pdSet
exprset <- exprs(dataSet) 
conditions <- pData(dataSet)$conditionName
uniqueConditions <- unique(conditions)
mytestcond <- uniqueConditions[uniqueConditions !=baseline]
normalizedDataName <- paste(outputDirectory,"/", runName, "_ALLnormalized_", mytestcond[1], "_vs_",baseline, ".txt", sep="") # Sets file name.
write.table(exprset, file=normalizedDataName, row.names=FALSE, quote=FALSE,sep="\t")

########################### End Normalization Block ##########################

ids <- rownames(exprset)
linkedprobes <- select(dbpackagename, keys= ids, columns = "SYMBOL", keytype = "PROBEID") 

########################### Begin Gage Block #########################
pathwayDirectory <- paste(outputDirectory,"/", pathwayDirectoryName, sep="")
system(paste("mkdir ", pathwayDirectory, sep = ""))
setwd(pathwayDirectory)

ids <- rownames(exprset)

# SELECT BY ENTREZIDS
linkedprobes2 <- select(dbpackagename, keys= ids, columns = "ENTREZID", keytype = "PROBEID")
exprsetlinkedtogenes2 <- merge(exprset, linkedprobes2, by.x=0, by.y="PROBEID")
col_1st_exprset2 <- grep("ENTREZID", names(exprsetlinkedtogenes2))
exprsetlinkedtogenes2 <- exprsetlinkedtogenes2[,c(col_1st_exprset2, (1:ncol(exprsetlinkedtogenes2))[-col_1st_exprset2])]
exprsetlinkedtogenes2[2]<-NULL
expmatrix<-data.matrix(exprsetlinkedtogenes2)
rownames(expmatrix)<-expmatrix[,1]
dataMatrix<-expmatrix[,-1]
geneList <- rownames(dataMatrix)
data(kegg.gs)
#Calculate the differential regulation of pathways comparing the baseline against the comparator
gageOutput <- gage(dataMatrix, gsets = kegg.gs, ref = NULL, samp = NULL, compare="as.group")
write.table(gageOutput$less,file="gageOutputLess_GE.xls", sep="\t")
write.table(gageOutput$greater,file="gageOutputGreater_GE.xls", sep="\t")
########################### End Gage Block #########################