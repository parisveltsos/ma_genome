data <- read.table("FST.R.Input.txt", header=TRUE)

nrow(data)

mydata <- na.omit(data)
nrow(mydata)

total_window <- nrow(mydata)

mydata <- subset(mydata, mydata$Length > 99)
total_window <- nrow(mydata)
nrow(mydata)

##Assigning Transcripts to Chromosomes###
MAChr1 <- subset(mydata, mydata$LG=="1")
MAChr2 <- subset(mydata, mydata$LG=="2")
MAChr3 <- subset(mydata, mydata$LG=="3")
MAChr4 <- subset(mydata, mydata$LG=="4")
MAChr5 <- subset(mydata, mydata$LG=="5")
MAChr6 <- subset(mydata, mydata$LG=="6")
MAChr7 <- subset(mydata, mydata$LG=="7")
MAChr8 <- subset(mydata, mydata$LG=="8")

MA_Autosomes <- subset(mydata, mydata$LG==2 | mydata$LG==3 | mydata$LG==4 | mydata$LG==5 | mydata$LG==6 | mydata$LG==7 | mydata$LG==8)

## Sort According to Position ####
MAChr1_sort <- MAChr1[order(MAChr1$cf_female),] 
MAChr2_sort <- MAChr2[order(MAChr2$cf_female),] 
MAChr3_sort <- MAChr3[order(MAChr3$cf_female),] 
MAChr4_sort <- MAChr4[order(MAChr4$cf_female),] 
MAChr5_sort <- MAChr5[order(MAChr5$cf_female),] 
MAChr6_sort <- MAChr6[order(MAChr6$cf_female),] 
MAChr7_sort <- MAChr7[order(MAChr7$cf_female),] 
MAChr8_sort <- MAChr8[order(MAChr8$cf_female),] 


####### Statistical Tests ########
wilcox.test(MAChr1$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr2$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr3$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr4$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr5$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr6$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr7$WeightFST,MA_Autosomes$WeightFST)
wilcox.test(MAChr8$WeightFST,MA_Autosomes$WeightFST)



### Sex-linked vs. Auto ###
Sex <- subset(mydata, mydata$Sex.Linkage== "SL")
Auto <- subset(mydata, mydata$Sex.Linkage == "Au")
PAR <- subset(mydata, mydata$Sex.Linkage == "PAR")

wilcox.test(Sex$WeightFST,Auto$WeightFST)
wilcox.test(Sex$WeightFST,PAR$WeightFST)
wilcox.test(PAR$WeightFST,Auto$WeightFST)


### plots sex-linked vs. auto ##
library(ggplot2)
library(ggpubr)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw(base_size=12)) ## No gray background

pdf("~/WeightFST_Boxplot.pdf", width=7,height=5)
WeightFST_Boxplot <- ggplot(mydata, aes(x = Sex.Linkage, y = WeightFST)) +
	geom_boxplot(outlier.shape=21, notch = TRUE) +
    scale_y_continuous(name = "FST") + theme_bw() + 
    theme(plot.title = element_blank(),  axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size=20), axis.ticks.x = element_blank()) + 
    theme(axis.text.x = element_text(color = "black", size=20)) + 
    theme(axis.text.y = element_text(color = "black", size=20)) +
    theme(text = element_text(size=20)) +
    theme(panel.grid = element_blank()) 
WeightFST_Boxplot
dev.off()



###### Sliding Window Analysis ####

library(zoo)

Chr1RM<- rollmean(smooth(MAChr1_sort$WeightFST),30)
Chr2RM<- rollmean(smooth(MAChr2_sort$WeightFST),30)
Chr3RM<- rollmean(smooth(MAChr3_sort$WeightFST),30)
Chr4RM<- rollmean(smooth(MAChr4_sort$WeightFST),30)
Chr5RM<- rollmean(smooth(MAChr5_sort$WeightFST),30)
Chr6RM<- rollmean(smooth(MAChr6_sort$WeightFST),30)
Chr7RM<- rollmean(smooth(MAChr7_sort$WeightFST),30)
Chr8RM<- rollmean(smooth(MAChr8_sort$WeightFST),30)


testRM <- c(Chr2RM, Chr3RM, Chr4RM, Chr5RM, Chr6RM, Chr7RM, Chr8RM)


myfunction <- function(i){
Info <- sample(i,1,replace=FALSE)
return(Info)
}

my.perm <- c()
for(i in 1:10^3){ my.perm[i] <- myfunction(testRM) }
sorted.perm <- sort(my.perm)
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]


RMpalette <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")



pdf("~/MA_LG1_WeightFST_RM.pdf", width=7,height=5)
MAChr1RM <- rollmean(smooth(MAChr1_sort$cf_female), 30)
MAWeightFSTChr1RM <- rollmean(smooth(MAChr1_sort$WeightFST),30)
plot(MAChr1_sort$cf_female, MAChr1_sort$WeightFST,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-0.05,0.3), xlab="Position(cM)", ylab="FST",main="Chr1",cex.main=1.8,cex.lab=1.3)
lines(MAChr1RM, MAWeightFSTChr1RM,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(v=46.44, lty=2)
abline(v=68.67, lty=2)
dev.off()