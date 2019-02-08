# load libraries
library(edgeR)
library(gplots)
library(ggplot2)
library(dynamicTreeCut)

# define colours for plotting
cbred <- 2 # '#D55E00'
cbblue <- '#0072B2'
cbgreen <- '#009E73'
cblblue <- '#56B4E9'

FDR2use  <- 0.05

datapath <- "~/git/ma_genome/input/de_analysis"
outpath <- "~/git/ma_genome/output/de_analysis"
dir.create(file.path(outpath))
# annotation <- read.delim(file.path(datapath, "ma_annotation.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE) # using the annotation slows down all output, so it is best explored separately later.

count <- read.table(file.path(datapath, paste(sub_analyse,'_count.txt', sep="")), header=T, row.names=1) 
count <- round(count, digits=0)
design <- read.table(file.path(datapath, paste(sub_analyse,'_design.txt', sep="")), header=T)

model.formula <- as.formula("~0+group")
dmat <- model.matrix(model.formula, data=as.data.frame(design))
# dgl <- DGEList(counts=count, group=design$group, genes=annotation)
dgl <- DGEList(counts=count, group=design$group, genes=row.names(count))

# investigate the effect of methodology of filtering reads with low counts
filter_file <- file.path(paste(outpath, '/',sub_analyse, 'filtering_info.txt', sep=""))
ave0 <- dgl[aveLogCPM(dgl) >= 0,]
ave1 <- dgl[aveLogCPM(dgl) >= 1,]
sum2 <- dgl[rowSums(cpm(dgl)) >= 2,]
sum10 <- dgl[rowSums(cpm(dgl)) >= 10,]

write(paste("Comp\tRemaining\tMin\t1Q\tMedian\tMean\t3Q\tMax"), filter_file)
write(paste("All_"), filter_file, append=T)
write(paste(nrow(dgl), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(dgl)/ncol(dgl))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Ave0_"), filter_file, append=T)
write(paste(nrow(ave0), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(ave0)/ncol(ave0))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Ave1_"), filter_file, append=T)
write(paste(nrow(ave1), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(ave1)/ncol(ave1))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Sum2_"), filter_file, append=T)
write(paste(nrow(sum2), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(sum2)/ncol(sum2))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Sum10_"), filter_file, append=T)
write(paste(nrow(sum10), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(sum10)/ncol(sum10))), filter_file, append=T, sep='\t', ncol=6)

# compile nice table with perl -pe 's/_\n/\t/g' mafiltering_info.txt > filter.out
# chosen filtering is here
dgl <- dgl[aveLogCPM(dgl) >= 0,]

summary(aveLogCPM(dgl))
summary(rowSums(dgl$count))

y <- dgl
colnames(y) <- paste(colnames(y), design$group, sep="\n")

# MDS plot {
## colours
ma_red   <- rgb(206/255, 34/255, 43/255, 3/4)
ma_blue <- rgb(100/255, 179/255, 223/255, 3/4)
pdf(file.path(outpath,paste('MDS_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
cols = c(ma_red,ma_red,ma_red,ma_red,ma_red,ma_blue,ma_red,ma_red,ma_red,ma_red,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_red,ma_red,ma_red,ma_red,ma_red,ma_red,ma_red,ma_red,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_red,ma_red,ma_red,ma_red,ma_red,ma_red,ma_blue,ma_blue,ma_blue,ma_blue,ma_blue,ma_red,ma_red,ma_red,ma_red,ma_red,ma_red,ma_blue,ma_blue,ma_blue,ma_blue,ma_red,ma_blue,ma_blue,ma_blue)
pchs = c(design$group)
plotMDS(y, cex=0.5, col=cols, main=paste(sub_analyse,"MDS plot"), cex.main=1.8, cex.lab=1.3)
plotMDS(y, cex=2, pch=as.numeric(y$samples$group), col=cols, main=paste(sub_analyse,"MDS plot"), cex.main=1.8, cex.lab=1.3)
legend('bottomleft', inset=0.05, legend=levels(design$group), pch = pchs, col=c(ma_red, ma_blue) )
dev.off()
}

# estimate data normalisation factors and dispersion
xcpm <- mglmOneGroup(dgl$counts)     # computing a logCPM for making dispersion plot
dgl <- calcNormFactors(dgl)
dgl <- estimateGLMCommonDisp(dgl, dmat)
dgl <- estimateGLMTrendedDisp(dgl, dmat, min.n=1000)
dgl <- estimateGLMTagwiseDisp(dgl, dmat)

## dispersion plot
pdf(file.path(outpath, paste('Dispersion_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plot(xcpm, dgl$tagwise.dispersion, pch=16, cex=0.5, xlab="log2CPM", ylab="Dispersion", main=paste(sub_analyse," dispersion", sep=""))
if(!is.null(dgl$trended.dispersion)) points(xcpm,dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
abline(h=dgl$common.dispersion ,col=cblblue, lwd=2)
legend("topright", c("Common","Trended","Tagwise"), pch=16, col=c(cblblue, cbgreen, "black"), title="Dispersion")
dev.off()

##  fit the data model ------------------------------------------
fitres <- glmFit(dgl, dmat)
x <- read.delim(paste(datapath, sub_analyse, "_matrix.txt", sep=""), sep="\t", header=T)
sortedX <- data.frame(x[order(x$model_coefficients, decreasing=F),])
cmat <- as.matrix(sortedX[,-1])
colnames(cmat)[1] <- colnames(sortedX[2])
rownames(cmat) <- as.character(sortedX[,1])
# cmat # to check the contrast applied: male-biased genes have positive values

lrtres <- list()
for(k in 1:ncol(cmat)) lrtres[[k]] <- glmLRT(fitres, contrast=cmat[,k])	

logFC <- NULL
PV <- NULL
FDR <- NULL
for(k in 1:ncol(cmat)) {PV <- cbind(PV, lrtres[[k]]$table[,"PValue"])
	FDR= cbind(FDR, p.adjust(PV[,k],method="BH"))
	logFC= cbind(logFC, lrtres[[k]]$table[,"logFC"])
	}

xcpm <- lrtres[[1]]$table[,"logCPM"]
allzeros <- which(rowSums(cpm(dgl)) < 3,)
allused <- which(rowSums(cpm(dgl)) >= 3,)

cname <-colnames(cmat)
colnames(logFC) <- paste("logFC", cname,sep=".")
colnames(PV) <- paste("PV", cname, sep=".")
colnames(FDR) <- paste("FDR", cname, sep=".")
rownames(logFC) <- rownames(PV) <- rownames(FDR) <- rownames(fitres$coefficients)

## Make results table
idxzeros <- allzeros
# restab <- data.frame(rownames(dgl$counts),dgl$genes[,2], dgl$genes[,2], dgl$genes[,2], dgl$genes[,2],logCPM=xcpm,logFC,PV,FDR) # note start and end are the same in annotation file
restab <- data.frame(rownames(dgl$counts),logCPM=xcpm,logFC,PV,FDR) # note start and end are the same in annotation file

colnames(restab)[1] <- 'gid' # commented lines are used if more complicated annotation is loaded. In practice it slowed down the code to print a lot of annotation information, and the annotation is better explored separately.
# colnames(restab)[2] <- 'gname'
# colnames(restab)[3] <- 'chr'
# colnames(restab)[4] <- 'start'
# colnames(restab)[5] <- 'end'

# write.table(restab$gid, file=file.path(outpath, paste(sub_analyse, "_used_gene_names.txt", sep="")), quote=F, row.names=F, sep="\t")

## pvalue histogram
pdf(file.path(outpath, paste('p_hist_',FDR2use, '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
npanel <- ncol(logFC)
np <- ceiling(sqrt(npanel))
if(np*(np-1)>= npanel) mfcol <- c(np-1,np) else
{if((np-1)*(np+1)>= npanel) mfcol <- c(np+1,np-1) else mfcol <- c(np,np)}
par(mfcol=mfcol)
for(k in 1:npanel) {
hist(PV[,k],n=100,xlab="P-value",main=colnames(logFC)[k])
}
dev.off()

## whinin group pairwise scatter plot
wx <- dgl$counts
wg <- as.character(dgl$samples$group)
wug <- unique(wg)
wn <- length(wug)
for (k in 1:wn) {
	ix <- wg %in% wug[k]
	xmat <- log2(wx[,ix])
	pdf(file.path(outpath, paste('pairwise_raw_count_', wug[k], '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
if (sum(ix) > 1 ) {
	pairs(xmat,pch=16,cex=0.4,main=wug[k])
	}
#	dev.copy(pdf,file.path(outpath,'courtship_out',paste(sub_analyse,wug[k],'_pairwise_raw_count.pdf', sep="")), width=8, height=8)
	dev.off()
}	

# write results based on different logFC thresholds
restab_frame <- as.data.frame(restab)

for(logFC_use in c(3, 2, 1, 0) ) {
de.yes.no <- FDR < FDR2use & abs(logFC) > logFC_use
if (ncol(cmat) == 1) {
	de4 <- which((FDR[,1] < FDR2use) == 0 & abs(logFC[,1] > logFC_use) != 0) } else {
de4 <- which(rowSums(FDR[,c(1:ncol(FDR))] < FDR2use) == 0 & rowSums(abs(logFC[,c(1:ncol(FDR))]) > logFC_use) != 0 )
}

de.yes.no[de4,] <- FALSE
deidx <- ii <- rowSums(de.yes.no) > 0
delabel <- (sign(logFC)*de.yes.no)[ii,]
delabel[is.na(delabel)] <- 0
combinedp <- 1; for(k in 1:ncol(PV)) combinedp <- combinedp*PV[,k]
deidx <- ii

wtable_1 <- restab[deidx,] # this is the result table of DE genes.
demat <- as.matrix(logFC[deidx,])
write.table(wtable_1, file=file.path(outpath, paste('de_',FDR2use, '_',logFC_use, '_', sub_analyse,'.txt', sep="")), quote=F, row.names=F, sep='\t')

if (nrow(wtable_1) > 1) {
	if (ncol(cmat) == 1) {
		wtable_3 <- rbind(NUM_DE=sum(abs(delabel)), NUM_UP_DE =sum(delabel>0), NUM_DOWN_DE =sum(delabel<0)) } else {
		wtable_3 <- rbind(NUM_DE=colSums(abs(delabel)), NUM_UP_DE =colSums(delabel>0), NUM_DOWN_DE =colSums(delabel<0))
	}
	wtable_3 <- cbind(DE_numbers=rownames(wtable_3), wtable_3)
	write.table(wtable_3, file=file.path(outpath, paste('Number_de_',FDR2use, '_',logFC_use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
	}
}

fc <- logFC
fc[(fc)>10] <- 10
fc[ fc < -10] <- -10
par(mfcol=mfcol)

for(k in 1:ncol(cmat)) {
pdf(file.path(outpath, paste('FC-CPMplot_',FDR2use, '_', sub_analyse,'_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=8, height=8)
ylab <- colnames(logFC)[k]
deix <- which(de.yes.no[,k])
maPlot(x=NULL,y=NULL,logAbundance= xcpm, logFC = fc[,k], xlab = bquote(paste(log^2, CPM)), ylab = paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' - ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), de.tags= deix,pch = 19, cex = 0.3, smearWidth = 0.5, panel.first = grid(), smooth.scatter = FALSE, lowess = FALSE,  main = paste('LogFC plot ', strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""))
dev.off()
}

### Output subsets FDR05 for GO annotation

### 5% FDR 

list_de <- list()
list_nonde <- list()

for(k in 1:ncol(cmat)) {
FDR05 <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use)
FDR05_unbiased <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] > FDR2use)
UP <- subset(FDR05, FDR05[[paste('logFC.',colnames(cmat)[k], sep="")]] > 0)
DOWN <- subset(FDR05, FDR05[[paste('logFC.',colnames(cmat)[k], sep="")]] < 0)

# write.table(list_de[[k]]$gid, file=file.path(outpath, paste('gene_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

write.table(data.frame(gene=FDR05_unbiased$gid, logFC=FDR05_unbiased[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=FDR05_unbiased[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath, paste('unbiased_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')
write.table(data.frame(gene=UP$gid, logFC=UP[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=UP[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath, paste('UP_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')
write.table(data.frame(gene=DOWN$gid, logFC=DOWN[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=DOWN[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath, paste('DOWN_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')

gopath <- file.path(outpath, paste('GO_',FDR2use, '_', colnames(cmat)[k], sep=""))
dir.create(gopath)

write.table(data.frame(gene=FDR05$gid), file=file.path(gopath,'FDR_genes.txt'), quote=F, row.names=F, sep='\t')
write.table(data.frame(gene=UP$gid, logFC=UP[[paste('logFC.',colnames(cmat)[k], sep="")]]), file=file.path(gopath, paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], '_biased.txt', sep="")), quote=F, row.names=F, sep="\t")
write.table(data.frame(gene=DOWN$gid, logFC=DOWN[[paste('logFC.',colnames(cmat)[k], sep="")]]), file=file.path(gopath, paste(strsplit(colnames(cmat)[k], '\\.')[[1]][2], '_biased.txt', sep="")), quote=F, row.names=F, sep="\t")

GO_pvalues <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]], logFC=restab[[paste('logFC.',colnames(cmat)[k], sep="")]])
write.table(GO_pvalues, file=file.path(gopath, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

GO_pvalues_UP <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]])
GO_pvalues_UP[,2][restab[[paste('logFC.',colnames(cmat)[k], sep="")]] < 0 & restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use] <- 1
write.table(GO_pvalues_UP, file=file.path(gopath, "GO_pvalues_UP.txt"), quote=F, row.names=F, sep="\t")

GO_pvalues_DOWN <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]])
GO_pvalues_DOWN[,2][restab[[paste('logFC.',colnames(cmat)[k], sep="")]] > 0 & restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use] <- 1
write.table(GO_pvalues_DOWN, file=file.path(gopath, "GO_pvalues_DOWN.txt"), quote=F, row.names=F, sep="\t")
}



# Calculation of moderated log-counts-per-million for use in heatmaps
nc <- cpm(dgl, prior.count=2, log=T)
nc2 <- data.frame(row.names(nc), nc)

colnames(nc2) <- c('ID', as.character(design$sample))
restab_logCPM = merge(restab, nc2, by.x="gid", by.y="ID", all=F )
# restab_logCPM$bias <- 'empty'

write.table(restab_logCPM, file=file.path(outpath, paste('LogCPM_',FDR2use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

# Heatmaps 

colnames(restab_logCPM)
tail(restab_logCPM)
cmat <- data.frame(cmat)

for(k in 1:ncol(cmat)) {
	if (ncol(cmat) == 1) {
		cmat_subset <- cmat } else {
		cmat_subset <- subset(cmat, cmat[[colnames(cmat)[k]]]!=0) 
		}
design_subset <- design[design$group %in% row.names(cmat_subset),]

length(as.character(design_subset$sample))
rownames(restab_logCPM) <- restab_logCPM$gid
DE_counts <- subset(restab_logCPM, restab_logCPM[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use)
DE_counts_relevant <- subset(DE_counts, select=as.character(design_subset$sample))
# row.names(DE_counts_relevant) <- DE_counts$gname
colnames(DE_counts_relevant) <- paste(colnames(DE_counts_relevant), '\n', design_subset$group) 
if (nrow(DE_counts_relevant) > 2) {

pdf(file.path(outpath, paste('Heatmap_', FDR2use, '_DE_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)

cmethods <- c('average', 'ward.D') # options are complete ward.D average median ward.D2 mcquitty centroid single

for(m in 1:length(cmethods)) {
d <- as.matrix(DE_counts_relevant)
myheatcol <- colorpanel(100, cbred,'white', cbblue) # choose a color palette for the heat map
distmatrix <- as.dist(1-cor(t(d), method="pearson"))
hr <- hclust(distmatrix, method=cmethods[m])  # plot(hr)
mycl <- cutreeDynamic(hr, distM=as.matrix(distmatrix)) # dynamic tree cut
clusterCols <- rainbow(length(unique(mycl))) # get a color palette equal to the number of clusters
myClusterSideBar <- as.character(as.numeric(mycl))
heatmap.2(d, col=myheatcol, Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.5, cexRow=0.6, RowSideColors = myClusterSideBar, main=paste(sub_analyse, cmethods[m], strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use), srtCol=45, key.title=NA)
legend(x=0.15,y=1.12, legend = unique(mycl), col = unique(as.numeric(mycl)), lty= 1, lwd = 3, cex=.5, title="clusters", xpd=T)
DE_counts[[ cmethods[m] ]] = as.factor(mycl)
gopath <- file.path(outpath, paste('GO_',FDR2use, '_', colnames(cmat)[k], sep=""))
dir.create(gopath)
gopath_method <- file.path(outpath, paste('GO_',FDR2use, '_', colnames(cmat)[k],'/',cmethods[m], sep=""))
dir.create(gopath_method)

for(l in 1:length( levels (DE_counts[[cmethods[m]]] ) ) ) {
FDR05 <- data.frame(gene=DE_counts$gid, cluster=DE_counts[[cmethods[m]]], FDR=DE_counts[[paste('FDR.',colnames(cmat)[k], sep="")]] )
FDR05[,3][FDR05[,3] < FDR2use & FDR05[,2] != l] <- 1
write.table(data.frame(gene=FDR05[,1], FDR=FDR05[,3]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')
}
}
dev.off()

print(paste(colnames(cmat)[k], 'groups:  median', length(levels(DE_counts$median)), '  ward', length(levels(DE_counts$ward)), '  complete', length(levels(DE_counts$complete)), '  average', length(levels(DE_counts$average)), '  ward2', length(levels(DE_counts$ward.D2)), '  mcquitty', length(levels(DE_counts$mcquitty)), '  centroid', length(levels(DE_counts$centroid)), '  centroid', length(levels(DE_counts$centroid)),sep=" "))

write.table(DE_counts, file=file.path(outpath, paste('clusters',FDR2use, '_', sub_analyse, '_',colnames(cmat)[k], '.txt', sep="")), quote=F, row.names=F, sep='\t')
}
}

