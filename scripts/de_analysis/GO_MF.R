# Rscript1.R
# eg. run: 
# Rscript Rscript1.R annotations.txt pvalues.txt output_file

args <- commandArgs(trailingOnly = TRUE)
universeFile = args[1]
interestingGenesFile = args[2]
output_file = args[3]
FDRinput = as.numeric(paste(args[4]))
	
# output_file = 'out_MF'
# FDRinput = 0.05

# set the output file

sink(output_file)

# load topGO
library("topGO")

# read in the 'gene universe' file - should only include all possible genes in the comparison of interest. It is a subset of all GO information available for all genes
geneID2GO <- readMappings(file = universeFile)
geneUniverse <- names(geneID2GO)
geneList <- read.table(interestingGenesFile, header=T)
geneUniverse_frame <- data.frame(geneUniverse)
merged1 = merge(geneList, geneUniverse_frame, by.x=colnames(geneList)[1], by.y="geneUniverse", all=F )
geneList_pvalues <- as.vector(merged1[,2])
names(geneList_pvalues) <- merged1[,1]
# names(geneList) <- geneUniverse

# function to identify interesting genes based on 'score' = pvalue (or could be logFC)
topDiffGenes <- function (score) {
    return(score < FDRinput)
}

# build the GOdata object in topGO

myGOdata <- new("topGOdata", description="My project", ontology="MF", nodeSize=5, allGenes=geneList_pvalues, geneSel=topDiffGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
# myGOdata

# run the Fisher's exact tests
resultFisherClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
# resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultFisherTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
# resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# run KS tests
resultKSClassic <- runTest(myGOdata, algorithm = "classic", statistic = "ks", scoreOrder = "increasing")
# resultKSElim <- runTest(myGOdata, algorithm = "elim", statistic = "ks", scoreOrder = "increasing")
resultKSTopgo <- runTest(myGOdata, algorithm = "weight01", statistic = "ks", scoreOrder = "increasing")


# see how many results we get where weight01 gives a P-value <= 0.01:
mysummary <- summary(attributes(resultKSTopgo)$score <= 0.01)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.01
mysummaryFisher <- summary(attributes(resultFisherTopgo)$score <= 0.01)
numsignifFisher <- as.integer(mysummaryFisher[[3]]) # how many terms is it true that P <= 0.01

# print out the top 'numsignif' results:
# allRes <- GenTable(myGOdata, classicFisher = resultFisherClassic, elimFisher = resultElim, topgoFisher = resultFisherTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
# allRes
# report classic Fisher, top go Fisher and topgo KS. Order by topgoFisher (best balance between simple and powerful approach) but report all KS significant p values, which tend to be more than Fisher.

if (numsignifFisher < 2) {numsignifFisher = 2}
if (numsignif < 2) {numsignif = 2}

allRes_Fisher <- GenTable(myGOdata, topgoFisher = resultFisherTopgo, topgoKS = resultKSTopgo, classicFisher = resultFisherClassic, orderBy = "topgoFisher", topNodes = numsignifFisher)
allRes_Fisher$cluster <- output_file
write.table(allRes_Fisher, file=paste("table_Fisher", output_file, sep=""), quote=F, row.names=F, sep="\t")

allRes_KS <- GenTable(myGOdata, topgoKS = resultKSTopgo, topgoFisher = resultFisherTopgo, classicFisher = resultFisherClassic, orderBy = "topgoKS", topNodes = numsignif)
allRes_KS$cluster <- output_file
write.table(allRes_KS, file=paste("table_KS", output_file, sep=""), quote=F, row.names=F, sep="\t")

 
# function to colour output
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
 
# pValue.KS.Elim <- score(resultKSElim)
pValue.KS.topgo <- score(resultKSTopgo)
pValue.Fisher.topGO <- score(resultFisherTopgo)[names(pValue.KS.topgo)]
gstat <- termStat(myGOdata, names(pValue.KS.topgo))
gSize <- gstat$Annotated/max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)

# comparison of KS with FISHER
# pdf(paste(output_file, 'MF_KS_vs_Fisher.pdf', sep=""))
# par(mar=c(5,5,4,3))
# plot(pValue.KS.topgo, pValue.Fisher.topGO, xlab = "Log p-value KS elim", ylab = "Log p-value Fisher topGO", pch = 19, cex = gSize, col = gCol, log='xy', main='Kolmogorov-Smirnov vs Fisher', cex.main=1.8, cex.lab=1.3)
# dev.off()

library("Rgraphviz")
# print a graph (to a pdf file) with the top 'numsignif' results:
pdf(paste(strsplit(output_file,'.txt')[[1]][1], '_KS.pdf', sep=""), width=20, height=20)
showSigOfNodes(myGOdata, score(resultKSTopgo), firstSigNodes = numsignif, useInfo = "all")
dev.off() 

pdf(paste(strsplit(output_file, '.txt')[[1]][1], '_Fisher.pdf', sep=""), width=20, height=20)
showSigOfNodes(myGOdata, score(resultFisherTopgo), firstSigNodes = numsignifFisher, useInfo = "all")
dev.off() 

# print out the genes that are annotated with the significantly enriched GO terms:
genesOfInterest <- as.character(geneList[,1][geneList[,2] < 0.05])

myterms <- allRes_Fisher$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
   myterm <- myterms[i]
   mygenesforterm <- mygenes[myterm][[1]] 
   myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
   mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
   mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
   print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
sink()


