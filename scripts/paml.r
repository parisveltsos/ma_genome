datapath <- '~/git/ma_genome/input/paml'
outpath <- '~/git/ma_genome/output/paml'
dir.create(outpath)
m1data <- read.table(file.path(datapath, 'lnL.m1.test.txt'), header=T) # paml model 1 output
m2data <- read.table(file.path(datapath, 'lnL.m2.test.txt'), header=T) # paml model 2 output
str(m1data) 
str(m2data) 

seqkey <- read.table(file.path(datapath, 'names.txt'), header=T) # paml analysis renamed the ORF names, the paml and SEX-DETector names are in this file

merged1 <- merge(m1data, seqkey, by.x='dataset', by.y='DataSet', all=T)
nrow(merged1)
str(merged1)

merged2 <- merge(merged1, m2data, by.x='dataset', by.y='dataset', all=T)
str(merged2)

topo1_only <- subset(merged2, merged2$topo==1) # topology 1 is (((X,Y), M. huetii), R. communis)
# topology 2 is (((X, M. huetii), Y), R. communis)
# topology 3 is (((Y, M. huetii), X), R. communis)
# topology 4 is ((X,Y),(M. huetii, R. communis))

# generate a data.frame with model 1 and model 2 ds, dn, omega (w) for x, y, M. huetii (h) and R. communis (r) and log likelihoods of the models. Letters in parentheses refer to column names below.
pamldata <- data.frame(merged2$Sequence, merged2$ds7..1.x, merged2$ds7..2.x, merged2$ds6..3.x, merged2$ds5..4.x, merged2$dn7..1.x, merged2$dn7..2.x, merged2$dn6..3.x, merged2$dn5..4.x, merged2$w7..1.x, merged2$w7..2.x, merged2$w6..3.x, merged2$w5..4.x,  merged2$ds7..1.y, merged2$ds7..2.y, merged2$ds6..3.y, merged2$ds5..4.y, merged2$dn7..1.y, merged2$dn7..2.y, merged2$dn6..3.y, merged2$dn5..4.y, merged2$w7..1.y, merged2$w7..2.y, merged2$w6..3.y, merged2$w5..4.y, merged2$topo, merged2$lnL_m0.x, merged2$lnL_m1, merged2$lnL_m2)

colnames(pamldata) <- c('name', 'm1xds', 'm1yds', 'm1hds', 'm1rds', 'm1xdn', 'm1ydn', 'm1hdn', 'm1rdn','m1xw', 'm1yw', 'm1hw', 'm1rw', 'm2xds', 'm2yds', 'm2hds', 'm2rds', 'm2xdn', 'm2ydn', 'm2hdn', 'm2rdn', 'm2xw', 'm2yw', 'm2hw', 'm2rw', 'topo', 'lnL_m0', 'lnL_m1', 'lnL_m2')

# filter to reliable estimates
xy <- subset(pamldata, pamldata$m2yds < 2 & pamldata$m2xds < 2 & pamldata$m2yds > 0.001 & pamldata$m2xds > 0.001 & pamldata$topo==2)

par(mfrow=c(1,3)) 
par(mar=c(5,5,4,3))
boxplot(xy$m2xw, xy$m2yw, xy$m2hw, xy$m2rw, names=list("X", "Y", "M. huetii", "R. communis"), ylab="dN/dS", cex.lab=1.8, ylim=c(0,1.6)) # max(xy$m2yw) # 2.7035
boxplot(xy$m2xds, xy$m2yds, xy$m2hds, xy$m2rds, names=list("X", "Y", "M. huetii", "R. communis"), ylab="dS", cex.lab=1.8, ylim=c(0,0.11)) # max(xy$m2hds) # 0.2292
boxplot(xy$m2xdn, xy$m2ydn, xy$m2hdn, xy$m2rdn, names=list("X", "Y", "M. huetii", "R. communis"), ylab="dN", cex.lab=1.8)
dev.copy(pdf,file.path(outpath,'paml.pdf'), width=18, height=6)
dev.off()
# Fig 4

# write.table(pamldata, file=file.path(datapath, "pamldata.txt"), quote=F, row.names=F, sep="\t")

# stats on M2 - comparison of X and Y dN, dS, omega.
m2 <- subset(pamldata, pamldata$m2xds < 2 & pamldata$m2yds < 2 & pamldata$m2xds > 0.001 & pamldata$m2yds > 0.001 & pamldata$topo==1)
# m1 <- subset(pamldata, pamldata$m1xds < 2 & pamldata$m1yds < 2 & pamldata$m2xds > 0.001 & pamldata$m2yds > 0.001 & pamldata$topo==1) # using model 1 does not significantly alter the patterns
nrow(m2)
par(mfrow=c(1,3)) 
par(mar=c(5,5,4,3))
boxplot(m2$m2xw, m2$m2yw, outline=F, notch=T, ylab = "dN/dS", names=c("X","Y"))
boxplot(m2$m2xds, m2$m2yds, outline=F, notch=T, ylab = "dS", names=c("X","Y"))
boxplot(m2$m2xdn, m2$m2ydn, outline=F, notch=T, ylab = "dN", names=c("X","Y"))
wilcox.test(xy$m2xw, xy$m2yw)
wilcox.test(xy$m2xds, xy$m2yds)
wilcox.test(xy$m2xdn, xy$m2ydn)

# Likelihood Ratio Test for model 2 being better than null model
xy$likelihood <- 2*abs(xy$lnL_m2 - xy$lnL_m0)
xy$chisq <- pchisq(xy$likelihood, df=1, lower.tail=FALSE)

# Counting Significant genes 
M2signif <- subset(xy, xy$chisq < 0.05)
nrow(M2signif)
M2RP <- subset(M2signif, M2signif$m2yw > M2signif$m2xw)
nrow(M2RP)

# BH Correction
M2signif$BH <- p.adjust(M2signif$chisq, method = "hochberg", n =74) #use 74, i.e. the total number of genes before testing.
M2signif.BH <- subset(M2signif, M2signif$BH< 0.05)
 nrow(M2signif.BH)
M2.BH.RP <- subset(M2signif.BH, M2signif.BH$m2yw > M2signif.BH$m2xw)
nrow(M2.BH.RP)


# Plot data to see the effect of topology

m2 <- subset(merged2, merged2$ds7..2.x < 2 & merged2$ds7..1.x < 2 & merged2$ds7..2.x > 0.001 & merged2$ds7..1.x > 0.001 & merged2$topo==1)

par(mar=c(5,5,4,3))
plot(m2$ds7..2.x, m2$dn7..2.x, type="n", xlab="dS", ylab="dN", main="title", cex.main=1.8, cex.lab=1.3)
points(m2$ds7..1.x, m2$dn7..1.x, pch=1, col=4)
points(m2$ds7..2.x, m2$dn7..2.x, pch=1, col=1, cex=0.8)
points(m2$ds7..1.x, m2$dn7..1.x, pch=17, col=4)
points(m2$ds7..2.x, m2$dn7..2.x, pch=17, col=1, cex=0.8)
legend('topright', inset=0.05, legend=c('X', 'Y', 'topo != 1'), pch =c(1,1, 17), col=c(4,1,1) )
