DataSNP <- read.table("SNPDensity.R.Input.txt", header=TRUE)

## Calculating Heterozygosity and Sex-biased SNP Density ##
DataSNP$MaleTotal <- (DataSNP$Male1Total + DataSNP$Male2Total + DataSNP$Male3Total + DataSNP$Male4Total + DataSNP$Male5Total + DataSNP$Male6Total + DataSNP$Male7Total + DataSNP$Male8Total + DataSNP$Male9Total + DataSNP$Male10Total + DataSNP$Male11Total + DataSNP$Male12Total + DataSNP$Male13Total + DataSNP$Male14Total + DataSNP$Male15Total + DataSNP$Male16Total + DataSNP$Male17Total + DataSNP$Male18Total + DataSNP$Male19Total + DataSNP$Male20Total)
DataSNP$MaleHomo <- (DataSNP$Male1homo + DataSNP$Male2homo + DataSNP$Male3homo + DataSNP$Male4homo + DataSNP$Male5homo + DataSNP$Male6homo + DataSNP$Male7homo + DataSNP$Male8homo + DataSNP$Male9homo + DataSNP$Male10homo + DataSNP$Male11homo + DataSNP$Male12homo + DataSNP$Male13homo + DataSNP$Male14homo + DataSNP$Male15homo + DataSNP$Male16homo + DataSNP$Male17homo + DataSNP$Male18homo + DataSNP$Male19homo + DataSNP$Male20homo)
DataSNP$MaleHet <- (DataSNP$Male1het + DataSNP$Male2het + DataSNP$Male3het + DataSNP$Male4het + DataSNP$Male5het + DataSNP$Male6het + DataSNP$Male7het + DataSNP$Male8het + DataSNP$Male9het + DataSNP$Male10het + DataSNP$Male11het + DataSNP$Male12het + DataSNP$Male13het + DataSNP$Male14het + DataSNP$Male15het + DataSNP$Male16het + DataSNP$Male17het + DataSNP$Male18het + DataSNP$Male19het + DataSNP$Male20het)
DataSNP$MaleDens <- (DataSNP$MaleHet + 1)/(DataSNP$MaleTotal + 1)
DataSNP$MaleDensL10 <- log10(DataSNP$MaleDens)
DataSNP$MalObsHet <- (DataSNP$MaleHet)/(DataSNP$MaleTotal)

DataSNP$FemTotal <- (DataSNP$Fem1Total + DataSNP$Fem2Total + DataSNP$Fem3Total + DataSNP$Fem4Total + DataSNP$Fem5Total + DataSNP$Fem6Total + DataSNP$Fem7Total + DataSNP$Fem8Total + DataSNP$Fem9Total + DataSNP$Fem10Total + DataSNP$Fem11Total + DataSNP$Fem12Total + DataSNP$Fem13Total + DataSNP$Fem14Total + DataSNP$Fem15Total + DataSNP$Fem16Total + DataSNP$Fem17Total + DataSNP$Fem18Total + DataSNP$Fem19Total + DataSNP$Fem20Total)
DataSNP$FemHomo <- (DataSNP$Fem1homo + DataSNP$Fem2homo + DataSNP$Fem3homo + DataSNP$Fem4homo + DataSNP$Fem5homo + DataSNP$Fem6homo + DataSNP$Fem7homo + DataSNP$Fem8homo + DataSNP$Fem9homo + DataSNP$Fem10homo + DataSNP$Fem11homo + DataSNP$Fem12homo + DataSNP$Fem13homo + DataSNP$Fem14homo + DataSNP$Fem15homo + DataSNP$Fem16homo + DataSNP$Fem17homo + DataSNP$Fem18homo + DataSNP$Fem19homo + DataSNP$Fem20homo)
DataSNP$FemHet <- (DataSNP$Fem1het + DataSNP$Fem2het + DataSNP$Fem3het + DataSNP$Fem4het + DataSNP$Fem5het + DataSNP$Fem6het + DataSNP$Fem7het + DataSNP$Fem8het + DataSNP$Fem9het + DataSNP$Fem10het + DataSNP$Fem11het + DataSNP$Fem12het + DataSNP$Fem13het + DataSNP$Fem14het + DataSNP$Fem15het + DataSNP$Fem16het + DataSNP$Fem17het + DataSNP$Fem18het + DataSNP$Fem19het + DataSNP$Fem20het)
DataSNP$FemDens <- (DataSNP$FemHet + 1)/(DataSNP$FemTotal + 1)
DataSNP$FemDensL10 <- log10(DataSNP$FemDens)
DataSNP$FemObsHet <- (DataSNP$FemHet)/(DataSNP$FemTotal)

DataSNP$log10MFRatio <- (DataSNP$MaleDensL10 - DataSNP$FemDensL10)


myvars <- c("Transcript", "MaleTotal", "MaleHomo", "MaleHet", "MaleDens", "MaleDensL10", "FemTotal", "FemHomo", "FemHet", "FemDens", "FemDensL10", "log10MFRatio", "MalObsHet", "FemObsHet")
newSNP <- DataSNP[myvars]

## Filtering for scaffolds using FST information ##
DataFST <- read.table("FST.R.Input.txt", header=TRUE)
nrow(DataFST)
mydata <- na.omit(DataFST)
nrow(mydata) 
total_window <- nrow(mydata)
mydata <- subset(mydata, mydata$Length > 99)
nrow(mydata)
merge <- merge(newSNP, mydata, by=c("Transcript"))


##Assigning Transcripts to LGs###
Autosomes <- subset(merge, merge$LG !="1")
Chr1 <- subset(merge, merge$LG=="1")
Chr2 <- subset(merge, merge$LG=="2")
Chr3 <- subset(merge, merge$LG=="3")
Chr4 <- subset(merge, merge$LG=="4")
Chr5 <- subset(merge, merge$LG=="5")
Chr6 <- subset(merge, merge$LG=="6")
Chr7 <- subset(merge, merge$LG=="7")
Chr8 <- subset(merge, merge$LG=="8")

## Sort According to Position ####
Chr1_sort <- Chr1[order(Chr1$cf_female),] 
Chr2_sort <- Chr2[order(Chr2$cf_female),] 
Chr3_sort <- Chr3[order(Chr3$cf_female),] 
Chr4_sort <- Chr4[order(Chr4$cf_female),] 
Chr5_sort <- Chr5[order(Chr5$cf_female),] 
Chr6_sort <- Chr6[order(Chr6$cf_female),] 
Chr7_sort <- Chr7[order(Chr7$cf_female),] 
Chr8_sort <- Chr8[order(Chr8$cf_female),] 


### Sex-linked vs. Auto Statistical Tests ###

Sex <- subset(merge, merge$Sex.Linkage== "SL")
Auto <- subset(merge, merge$Sex.Linkage == "Au")
PAR <- subset(merge, merge$Sex.Linkage == "PAR")

wilcox.test(Sex$log10MFRatio,Auto$log10MFRatio)
wilcox.test(Sex$log10MFRatio,PAR$log10MFRatio)
wilcox.test(PAR$log10MFRatio,Auto$log10MFRatio)


### plots sex-linked vs. auto ##
library(ggplot2)
library(ggpubr)
theme_set(theme_bw(base_size=12)) ## No gray background

pdf("~/SNPDensity_Boxplot.pdf", width=7,height=5)
log10MFRatio_Boxplot <- ggplot(merge, aes(x = Sex.Linkage, y = log10MFRatio)) +
	geom_boxplot(outlier.shape=21, notch = TRUE) +
    scale_y_continuous(name = "log10MFRatio") + theme_bw() +
    theme(plot.title = element_blank(),  axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size=20), axis.ticks.x = element_blank()) + 
    theme(axis.text.x = element_text(color = "black", size=20)) + 
    theme(axis.text.y = element_text(color = "black", size=20)) 
log10MFRatio_Boxplot
dev.off()



###### Sliding Window Analysis ####
library(zoo)

Chr1RM<- rollmean(smooth(Chr1_sort$log10MFRatio),30)
Chr2RM<- rollmean(smooth(Chr2_sort$log10MFRatio),30)
Chr3RM<- rollmean(smooth(Chr3_sort$log10MFRatio),30)
Chr4RM<- rollmean(smooth(Chr4_sort$log10MFRatio),30)
Chr5RM<- rollmean(smooth(Chr5_sort$log10MFRatio),30)
Chr6RM<- rollmean(smooth(Chr6_sort$log10MFRatio),30)
Chr7RM<- rollmean(smooth(Chr7_sort$log10MFRatio),30)
Chr8RM<- rollmean(smooth(Chr8_sort$log10MFRatio),30)


testRM <- c(Chr2RM, Chr3RM, Chr4RM, Chr5RM, Chr6RM, Chr7RM, Chr8RM)


myfunction <- function(i){
Info <- sample(i,1,replace=FALSE)
return(Info)
}

my.perm <- c()
for(i in 1:10^3){ my.perm[i] <- myfunction(testRM) }
sorted.perm <- sort(my.perm)

## Construction Confidence Intervals ##
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]



RMpalette <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")

pdf("~/MA_LG1_log10MFRatio_RM.pdf", width=7,height=5)
Chr1RM <- rollmean(smooth(Chr1_sort$cf_female), 30)
MAlog10MFRatioChr1RM <- rollmean(smooth(Chr1_sort$log10MFRatio),30)
plot(Chr1_sort$cf_female, Chr1_sort$log10MFRatio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-1,1), xlab="Position(CM)", ylab="log10MFRatio",main="Chr1",cex.main=1.8,cex.lab=1.3)
lines(Chr1RM, MAlog10MFRatioChr1RM,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(v=46.44, lty=2)
abline(v=68.67, lty=2)
dev.off()


pdf("~/MA_LG1_ObsHet_RM.pdf", width=7,height=5)
Chr1RM <- rollmean(smooth(Chr1_sort$cf_female), 30)
plot(Chr1RM, MAFemObsHetChr1RM,col="red",type="l", lwd=5, ylim=c(0,0.025), xlab="Position(CM)", ylab="ObsHet",main="Chr1",cex.main=1.8,cex.lab=1.3)
lines(Chr1RM, MAMalObsHetChr1RM,type="l",lwd=5, col=RMpalette[5])
abline(v=46.44, lty=2)
abline(v=68.67, lty=2)
dev.off()
