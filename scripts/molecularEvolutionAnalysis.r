# Setup
library(zoo)
library(ggridges)
library(AID)
library(MASS)
library(car)
library(sjPlot)
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(grid) 
theme_set(theme_bw(base_size=20)) ## No gray background 

# Input data and merge in a single data.frame
datapath <- '~/git/ma_genome/input/dnds'
outpath <- '~/git/ma_genome/output/dnds'
dir.create(outpath)
xysl <- read.table(file.path(datapath, 'sex-linked_dnds.txt'), header=T) # pairwise (X,Y) dN, dS data
xysl$xy_dnds <- xysl$xy_dn/xysl$xy_ds #  calculation of dN/dS ratio
str(xysl) 
xy <- subset(xysl, xysl$xy_dnds!='-Inf') # keep non-infinite values
nrow(xy) 

an_ric <- read.table(file.path(datapath, 'an_ric.txt'), header=T) # pairwise comparison with R. communis
str(an_ric) 

an_hu <- read.table(file.path(datapath, 'an_hu.txt'), header=T) # pairwise comparison with M. huetii
str(an_hu)
an_ric_name <- read.table(file.path(datapath, 'an_ric_name.txt'), header=T) # name correspondence of paml and transcriptome homologue files
str(an_ric_name)
an_hu_name <- read.table(file.path(datapath, 'an_hu_name.txt'), header=T) # name correspondence of paml and transcriptome homologue files
str(an_hu_name)

an_ric_complete <- merge(an_ric, an_ric_name, by.x='an_ric', by.y='an_ric', all=F) # combine name and values
an_hu_complete <- merge(an_hu, an_hu_name, by.x='an_hu', by.y='an_hu', all=F)

an_hu_ric <- merge(an_hu_complete, an_ric_complete, by.x='gene', by.y='gene', all=T) # combine M. annua and M. huetii data

an_hu_ric_xy <- merge(an_hu_ric, xy, by.x='gene', by.y='gene', all=T) # combine with R. communis

anno <- read.table('~/git/ma_genome/input/ma_annotation.txt', header=T) # this is file S4 - Transcriptome annotation. It is used to identify candidate genes

allele_sp_data <- read.table(file.path(datapath, 'expression_stop_data.txt'), header=T) # allele-specific expression data from SEX-DETector

anno2 <- merge(anno, allele_sp_data, by.x="ID", by.y="id", all=T)
all <- merge(an_hu_ric_xy, anno2, by.x='gene', by.y='ID', all=T) # merge everything

write.table(all, file=file.path(datapath, 'all_out.txt'), quote=T, row.names=F, sep='\t')

malemap <- read.table(file.path(datapath, 'male_map_SD.txt'), header=T) # male map, from LepMap3
femalemap <- read.table(file.path(datapath, 'female_map_SD.txt'), header=T) # female map, from LepMap3
colnames(femalemap)[2] <- "cf" # rename the female recombination distance to 'cf'

map <- xy <- merge(malemap, femalemap, by.x="id", by.y="id", all=F)
colnames(map)[2] <- "lg"
colnames(map)[4] <- "sd"

map_data <- merge(map, all, by.x='id', by.y='gene', all=T)
expr_data <- read.table(file.path(datapath, 'mf_expr_data.txt'), header=T) # expression data from edgeR

plot_data <- merge(expr_data, map_data, by.x='gene', by.y='id', all=T)


# dN/dS vs expression - no obvious effect
dev.copy(pdf,file.path(outpath,'dnds_expression_w_outliers_male.pdf'), width=18, height=12)
par(mfrow=c(2,3))
par(mar=c(5,5,4,3))
plot(plot_data$logfc, plot_data$xy_dnds, xlab="M/F logFC", ylab="XY dN/dS", main="Sex bias vs sex chromosome omega", cex.main=1.6, cex.lab=1.3, pch=16)

plot(plot_data$logfc, plot_data$an_hu_dnds, xlab="M/F logFC", ylab="M. annua - M. huetii dN/dS", main="Sex bias vs omega (with M. huetii)", cex.main=1.6, cex.lab=1.3,  pch=16, col=1)
points(plot_data$logfc[plot_data$SEX.DETector=='Au'], plot_data$an_hu_dnds[plot_data$SEX.DETector=='Au'], col=4, pch=16)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5_f4m6'], plot_data$an_hu_dnds[plot_data$SEX.DETector=='f3m5_f4m6'], col=2, pch=16)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5'], plot_data$an_hu_dnds[plot_data$SEX.DETector=='f3m5'], col=2, pch=16)
points(plot_data$logfc[plot_data$SEX.DETector=='_f4m6'], plot_data$an_hu_dnds[plot_data$SEX.DETector=='_f4m6'], col=2, pch=16)
legend('topright', inset=0.05, legend=c('Expressed', 'Autosomal', 'Sex-linked'), pch =16, col=c(1,4,2) ) 

plot(plot_data$logfc, plot_data$an_ric_dnds, xlab="M/F logFC", ylab="M. annua - R. communis dN/dS", main="Sex bias vs omega (with R. communis)", cex.main=1.6, cex.lab=1.3, pch=16, col=1)
points(plot_data$logfc[plot_data$SEX.DETector=='Au'], plot_data$an_ric_dnds[plot_data$SEX.DETector=='Au'], pch=16, col=4)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5_f4m6'], plot_data$an_ric_dnds[plot_data$SEX.DETector=='f3m5_f4m6'], pch=16, col=16)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5'], plot_data$an_ric_dnds[plot_data$SEX.DETector=='f3m5'], pch=16, col=4)
points(plot_data$logfc[plot_data$SEX.DETector=='_f4m6'], plot_data$an_ric_dnds[plot_data$SEX.DETector=='_f4m6'], pch=16, col=4)
legend('topright', inset=0.05, legend=c('Expressed', 'Autosomal', 'Sex-linked'), pch =16, col=c(1,4,2) ) 

plot(plot_data$logfc, plot_data$xy_ds, xlab="M/F logFC", ylab="XY dS", main="Sex bias vs sex chromosome divergence (synonymous)", cex.main=1.6, cex.lab=1.3, pch=16)

plot(plot_data$logfc, plot_data$an_hu_ds, xlab="M/F logFC", ylab="M.annua - M.huetii dN", main="Sex bias vs divergence rate (synonymous)", cex.main=1.6, cex.lab=1.3, pch=16, col=1)
points(plot_data$logfc[plot_data$SEX.DETector=='Au'], plot_data$an_hu_ds[plot_data$SEX.DETector=='Au'], col=4, pch=16)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5_f4m6'], plot_data$an_hu_ds[plot_data$SEX.DETector=='f3m5_f4m6'], col=2, pch=16)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5'], plot_data$an_hu_ds[plot_data$SEX.DETector=='f3m5'], col=2, pch=16)
points(plot_data$logfc[plot_data$SEX.DETector=='_f4m6'], plot_data$an_hu_ds[plot_data$SEX.DETector=='_f4m6'], col=2, pch=16)
legend('topright', inset=0.05, legend=c('Expressed', 'Autosomal', 'Sex-linked'), pch =16, col=c(1,4,2) ) 

plot(plot_data$logfc, plot_data$an_ric_ds, xlab="M/F logFC", ylab="M.annua - R.communis dN", main="Sex bias vs divergence rate (synonymous)", cex.main=1.6, cex.lab=1.3, pch=16, col=1)
points(plot_data$logfc[plot_data$SEX.DETector=='Au'], plot_data$an_ric_ds[plot_data$SEX.DETector=='Au'], pch=16, col=4)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5_f4m6'], plot_data$an_ric_ds[plot_data$SEX.DETector=='f3m5_f4m6'], pch=16, col=16)
points(plot_data$logfc[plot_data$SEX.DETector=='f3m5'], plot_data$an_ric_ds[plot_data$SEX.DETector=='f3m5'], pch=16, col=4)
points(plot_data$logfc[plot_data$SEX.DETector=='_f4m6'], plot_data$an_ric_ds[plot_data$SEX.DETector=='_f4m6'], pch=16, col=4)
legend('topright', inset=0.05, legend=c('Expressed', 'Autosomal', 'Sex-linked'), pch =16, col=c(1,4,2) ) 
dev.off()

# LG1 plots
### merge pi data to the other data to plot
pidata <- read.table('~/git/ma_genome/input/dnds/pi_data.txt', header=T)
pidata2 <- data.frame(tapply(pidata$PI, pidata$CHROM, mean)) 
pidata2$id <- rownames(pidata2)
colnames(pidata2) <- c('pi', 'id')

plot_data2 <- merge(plot_data, pidata2, by.x='gene', by.y='id', all=T)
write.table(plot_data2, file=file.path(outpath, "plot_data2.txt"), quote=F, row.names=F, sep="\t")

## Map and dataset setup
### Identify positions of chromosome start on the male and female recombination map, when all LG lengths are added together in order of LG. This makes it easy to plot all LGs.
m_chr2_start <- max(plot_data2$cm[plot_data2$lg=="1"])
m_chr3_start <- max(plot_data2$cm[plot_data2$lg=="2"]) + m_chr2_start
m_chr4_start <- max(plot_data2$cm[plot_data2$lg=="3"]) + m_chr3_start
m_chr5_start <- max(plot_data2$cm[plot_data2$lg=="4"]) + m_chr4_start
m_chr6_start <- max(plot_data2$cm[plot_data2$lg=="5"]) + m_chr5_start
m_chr7_start <- max(plot_data2$cm[plot_data2$lg=="6"]) + m_chr6_start
m_chr8_start <- max(plot_data2$cm[plot_data2$lg=="7"]) + m_chr7_start

f_chr2_start <- max(plot_data2$cf[plot_data2$lg=="1"])
f_chr3_start <- max(plot_data2$cf[plot_data2$lg=="2"]) + f_chr2_start
f_chr4_start <- max(plot_data2$cf[plot_data2$lg=="3"]) + f_chr3_start
f_chr5_start <- max(plot_data2$cf[plot_data2$lg=="4"]) + f_chr4_start
f_chr6_start <- max(plot_data2$cf[plot_data2$lg=="5"]) + f_chr5_start
f_chr7_start <- max(plot_data2$cf[plot_data2$lg=="6"]) + f_chr6_start
f_chr8_start <- max(plot_data2$cf[plot_data2$lg=="7"]) + f_chr7_start

# subset to LG-specific data
chr1_data_pi  <- subset(plot_data2, plot_data2$lg=="1")
chr1_data_pi <- data.frame(chr1_data_pi[order(chr1_data_pi$cf, decreasing=T),])
chr2_data_pi  <- subset(plot_data2, plot_data2$lg=="2")
chr3_data_pi  <- subset(plot_data2, plot_data2$lg=="3")
chr4_data_pi  <- subset(plot_data2, plot_data2$lg=="4")
chr5_data_pi  <- subset(plot_data2, plot_data2$lg=="5")
chr6_data_pi  <- subset(plot_data2, plot_data2$lg=="6")
chr7_data_pi  <- subset(plot_data2, plot_data2$lg=="7")
chr8_data_pi  <- subset(plot_data2, plot_data2$lg=="8")

## Special gene setup
 # Attempt to identify candidate sex determining or sexually antagonistic genes based on their annotation. These are used for plotting and chisq tests of their enrichment amongst sex-biased genes or on the non-recombining region of the sex chromosome.
 # These gene names are identified by grep with the term of the subset name, e.g. for auxin
 # grep auxin plot_data2.txt | cut -f 5 | grep 1 | wc -l
 # they are summed up to identify the number of genes on the sex chromosome compared to the autosomes

chr1_data_pi_auxin_subset <- subset(chr1_data_pi, 
							  chr1_data_pi$gene=='comp16624_c0_seq1_m.9076' |
							  chr1_data_pi$gene=='comp18115_c0_seq1_m.9953' |
							  chr1_data_pi$gene=='comp6958_c1_seq1_m.28384' |
							  chr1_data_pi$gene=='comp12005_c0_seq1_m.29333' |
							  chr1_data_pi$gene=='comp12640_c0_seq1_m.6320' |
							  chr1_data_pi$gene=='comp13243_c0_seq1_m.27340' |
							  chr1_data_pi$gene=='comp13243_c0_seq2_m.27341' |
							  chr1_data_pi$gene=='comp17892_c0_seq3_m.35488' |
							  chr1_data_pi$gene=='comp18797_c0_seq1_m.29884' |
							  chr1_data_pi$gene=='comp18797_c0_seq3_m.29885' |
							  chr1_data_pi$gene=='comp18902_c0_seq32_m.21806' |
							  chr1_data_pi$gene=='comp6958_c1_seq1_m.28384' |
							  chr1_data_pi$gene=='comp9501_c1_seq1_m.36203')

chr1_data_pi_cytokinin_subset <- subset(chr1_data_pi,
										chr1_data_pi$gene=='comp11994_c0_seq1_m.5841' | 
										chr1_data_pi$gene=='comp11086_c0_seq1_m.15413')

chr1_data_pi_flower_subset <- subset(chr1_data_pi, 
							  chr1_data_pi$gene=='comp7940_c0_seq1_m.3233' |
							  chr1_data_pi$gene=='comp9002_c0_seq1_m.3853' |
							  chr1_data_pi$gene=='comp18438_c0_seq1_m.20340')

chr1_data_pi_pollen_subset <- subset(chr1_data_pi, chr1_data_pi$gene=='comp4768_c0_seq1_m.1436')

chr1_data_pi_stop_subset <- subset(chr1_data_pi, chr1_data_pi$gene=='comp16827_c0_seq1_m.24466')

chr1_data_pi_kelch_subset <- subset(chr1_data_pi, 
							  chr1_data_pi$gene=='comp4002_c0_seq1_m.877' |
							  chr1_data_pi$gene=='comp15179_c0_seq1_m.8114' |
							  chr1_data_pi$gene=='comp16827_c0_seq1_m.24466' | 
							  chr1_data_pi$gene=='comp14507_c0_seq1_m.7652' |
							  chr1_data_pi$gene=='comp53760_c0_seq1_m.11659')
		
chr1_data_pi_ethylene_subset <- subset(chr1_data_pi, 
							  chr1_data_pi$gene=='comp16624_c0_seq1_m.9076' |
							  chr1_data_pi$gene=='comp12481_c0_seq1_m.30717' |
							  chr1_data_pi$gene=='comp14102_c0_seq1_m.7345' |
							  chr1_data_pi$gene=='comp14102_c0_seq1_m.7346' |
							  chr1_data_pi$gene=='comp15893_c0_seq1_m.8587' |
							  chr1_data_pi$gene=='comp17370_c1_seq1_m.36993' |
							  chr1_data_pi$gene=='comp10518_c0_seq1_m.4836' |
							  chr1_data_pi$gene=='comp21517_c0_seq1_m.10731' |
							  chr1_data_pi$gene=='comp12777_c0_seq1_m.6425' |
							  chr1_data_pi$gene=='comp11613_c0_seq1_m.5564' |
							  chr1_data_pi$gene=='comp140695_c0_seq1_m.13813' |
							  chr1_data_pi$gene=='comp14976_c0_seq1_m.7969' |
							  chr1_data_pi$gene=='comp15893_c0_seq1_m.8589' |
							  chr1_data_pi$gene=='comp17334_c0_seq1_m.21089' |
							  chr1_data_pi$gene=='comp7519_c0_seq1_m.2984' |
							  chr1_data_pi$gene=='comp9185_c0_seq1_m.26430')

chr1_data_pi_swisnf_subset <- subset(chr1_data_pi, 
									 chr1_data_pi$gene=='comp13123_c0_seq1_m.666' | 
									 chr1_data_pi$gene=='comp14497_c0_seq1_m.7644')

chr1_data_pi_retinoldehydrogenase_subset <- subset(chr1_data_pi, chr1_data_pi$gene=='comp16799_c0_seq1_m.30095')

chr1_data_pi_calmodulin_subset <- subset(chr1_data_pi,
							  chr1_data_pi$gene=='comp16900_c0_seq1_m.9236' |
							  chr1_data_pi$gene=='comp7296_c0_seq1_m.2850' |
							  chr1_data_pi$gene=='comp9071_c0_seq1_m.3902' |
							  chr1_data_pi$gene=='comp10223_c0_seq2_m.29540' |
							  chr1_data_pi$gene=='comp11508_c0_seq2_m.36042' |
							  chr1_data_pi$gene=='comp14313_c1_seq1_m.34636' |
							  chr1_data_pi$gene=='comp14483_c0_seq1_m.7633' |
							  chr1_data_pi$gene=='comp17025_c0_seq1_m.9304' |
							  chr1_data_pi$gene=='comp9071_c0_seq1_m.3902' |
							  chr1_data_pi$gene=='comp9449_c0_seq1_m.4145')

chr1_data_pi_short.root_subset <- subset(chr1_data_pi, 
										chr1_data_pi$gene=='comp5548_c0_seq1_m.1902' |
										chr1_data_pi$gene=='comp14712_c0_seq1_m.38500' |
										chr1_data_pi$gene=='comp14712_c0_seq1_m.38500' )

chr1_data_pi_jasmonicacid_subset <- subset(chr1_data_pi, chr1_data_pi$gene=='comp6853_c0_seq1_m.20904') 

 # fdr < 0.05 indicating significant sex biased genes
chr1_data_pi_FDR <- subset(chr1_data_pi, chr1_data_pi$fdr < 0.05) 

# plotting
 ## these plots are first defined and then plotted as separate panels of a large figure.

## logFC plot
logfc_lg1_y_title <- expression(paste("Male/Female logFC"))

logfc_male_lg1_plot <- ggplot() +
  geom_point(data=chr1_data_pi, aes(x = cm, y = logfc), color = cbPalette[3], alpha = 0.3, size=10) +
  geom_point(data=chr1_data_pi_FDR, aes(x = cm, y = logfc), color = cbPalette[2], alpha = 0.3, size=10) + 
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cm, y = logfc), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cm, y = logfc), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cm, y = logfc), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cm, y = logfc), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cm, y = logfc), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cm, y = logfc), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cm, y = logfc), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cm, y = logfc), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cm, y = logfc), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cm, y = logfc), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cm, y = logfc), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cm, y = logfc), label="j", color = 1, alpha = 1, size=9) + 
  labs(x = "Male map (cM)", y = logfc_lg1_y_title) +
  geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) +
#  geom_smooth(data=chr1_data_pi, aes(x = cm, y = logfc), label="s", color = cbPalette[6] ) +
  geom_line(data=chr1_data_pi, aes(x=cm, y=rollmean(smooth(logfc), 50, na.pad=T)), color = cbPalette[6], size=2) 


### logFC plot all LGs
 # logfc_male_plot <- ggplot() + 
 # geom_point(data=chr1_data_pi, aes(x = cm, y = logfc), color = cbPalette[3], alpha = 0.3) + 
 # geom_smooth(data=chr1_data_pi, aes(x = cm, y = logfc), color = cbPalette[6] ) + 
 # labs(x = "Male map (cM)") +
 # geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) +
 # geom_point(data=chr2_data_pi, aes(x = cm + m_chr2_start, y = logfc), color = cbPalette[6], alpha = 0.3) + 
 # geom_smooth(data=chr2_data_pi, aes(x = cm + m_chr2_start, y = logfc), color = cbPalette[3]) +
 # geom_point(data=chr3_data_pi, aes(x = cm + m_chr3_start, y = logfc), color = cbPalette[3], alpha = 0.3) + 
 # geom_point(data=chr4_data_pi, aes(x = cm + m_chr4_start, y = logfc), color = cbPalette[6], alpha = 0.3) + 
 # geom_point(data=chr5_data_pi, aes(x = cm + m_chr5_start, y = logfc), color = cbPalette[3], alpha = 0.3) + 
 # geom_point(data=chr6_data_pi, aes(x = cm + m_chr6_start, y = logfc), color = cbPalette[6], alpha = 0.3) + 
 # geom_point(data=chr7_data_pi, aes(x = cm + m_chr7_start, y = logfc), color = cbPalette[3], alpha = 0.3) + 
 # geom_point(data=chr8_data_pi, aes(x = cm + m_chr8_start, y = logfc), color = cbPalette[6], alpha = 0.3) + 
 # geom_smooth(data=chr3_data_pi, aes(x = cm + m_chr3_start, y = logfc), color = cbPalette[6]) +
 # geom_smooth(data=chr4_data_pi, aes(x = cm + m_chr4_start, y = logfc), color = cbPalette[3]) +
 # geom_smooth(data=chr5_data_pi, aes(x = cm + m_chr5_start, y = logfc), color = cbPalette[6]) +
 # geom_smooth(data=chr6_data_pi, aes(x = cm + m_chr6_start, y = logfc), color = cbPalette[3]) +
 # geom_smooth(data=chr7_data_pi, aes(x = cm + m_chr7_start, y = logfc), color = cbPalette[6]) +
 # geom_smooth(data=chr8_data_pi, aes(x = cm + m_chr8_start, y = logfc), color = cbPalette[3]) 

logfc_female_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cf, y = logfc), color = cbPalette[8], alpha = 0.3, size=10) + 
  geom_point(data=chr1_data_pi_FDR, aes(x = cf, y = logfc), color = cbPalette[2], alpha = 0.3, size=10) + 
 # geom_point(data=chr1_data_pi_subset, aes(x = cf, y = logfc), color = 1, alpha = 1) + 
  geom_line(data=chr1_data_pi, aes(x=cf, y=rollmean(smooth(logfc), 50, na.pad=T)), color = cbPalette[3], size=2) +
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cf, y = logfc), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cf, y = logfc), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cf, y = logfc), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cf, y = logfc), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cf, y = logfc), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cf, y = logfc), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cf, y = logfc), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cf, y = logfc), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cf, y = logfc), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cf, y = logfc), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cf, y = logfc), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cf, y = logfc), label="j", color = 1, alpha = 1, size=9) + 
#  geom_smooth(data=chr1_data_pi, aes(x = cf, y = logfc), color = cbPalette[1] ) + 

  labs(x = "Female map (cM)", y = logfc_lg1_y_title) +
  geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
  geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) 

 #logfc_female_plot <- ggplot() + 
 # geom_point(data=chr1_data_pi, aes(x = cf, y = logfc), color = cbPalette[8], alpha = 0.3) + 
 # geom_smooth(data=chr1_data_pi, aes(x = cf, y = logfc), color = cbPalette[1] ) + 
 # labs(x = "Female map (cM)") +
 # geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
 # geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) +
 # geom_point(data=chr2_data_pi, aes(x = cf + f_chr2_start, y = logfc), color = cbPalette[1], alpha = 0.3) + 
 # geom_smooth(data=chr2_data_pi, aes(x = cf + f_chr2_start, y = logfc), color = cbPalette[8]) +
 # geom_point(data=chr3_data_pi, aes(x = cf + f_chr3_start, y = logfc), color = cbPalette[1], alpha = 0.3) + 
 # geom_point(data=chr4_data_pi, aes(x = cf + f_chr4_start, y = logfc), color = cbPalette[8], alpha = 0.3) + 
 # geom_point(data=chr5_data_pi, aes(x = cf + f_chr5_start, y = logfc), color = cbPalette[1], alpha = 0.3) + 
 # geom_point(data=chr6_data_pi, aes(x = cf + f_chr6_start, y = logfc), color = cbPalette[8], alpha = 0.3) + 
 # geom_point(data=chr7_data_pi, aes(x = cf + f_chr7_start, y = logfc), color = cbPalette[1], alpha = 0.3) + 
 # geom_point(data=chr8_data_pi, aes(x = cf + f_chr8_start, y = logfc), color = cbPalette[8], alpha = 0.3) + 
 # geom_smooth(data=chr3_data_pi, aes(x = cf + f_chr3_start, y = logfc), color = cbPalette[8]) +
 # geom_smooth(data=chr4_data_pi, aes(x = cf + f_chr4_start, y = logfc), color = cbPalette[1]) +
 # geom_smooth(data=chr5_data_pi, aes(x = cf + f_chr5_start, y = logfc), color = cbPalette[8]) +
 # geom_smooth(data=chr6_data_pi, aes(x = cf + f_chr6_start, y = logfc), color = cbPalette[1]) +
 # geom_smooth(data=chr7_data_pi, aes(x = cf + f_chr7_start, y = logfc), color = cbPalette[8]) +
 # geom_smooth(data=chr8_data_pi, aes(x = cf + f_chr8_start, y = logfc), color = cbPalette[1]) 

png(file.path(outpath, 'logFC_lg1.png'), width=900, height=550, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 6,6), "null"))))
grid.text("Sex bias LG1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(logfc_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(logfc_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
dev.off()

## Pi plot
pi_lg1_y_title <- expression(paste(italic("π")))
chr1_data_pi_sub <- subset(chr1_data_pi, chr1_data_pi$pi!='NA') 
pi_male_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cm, y = pi), color = cbPalette[3], alpha = 0.3, size=10) + 
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cm, y = pi), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cm, y = pi), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cm, y = pi), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cm, y = pi), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cm, y = pi), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cm, y = pi), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cm, y = pi), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cm, y = pi), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cm, y = pi), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cm, y = pi), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cm, y = pi), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cm, y = pi), label="j", color = 1, alpha = 1, size=9) + 
  geom_smooth(data=chr1_data_pi, aes(x = cm, y = pi), color = cbPalette[6] ) + 
  labs(x = "Male map (cM)", y = pi_lg1_y_title) +
  geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) 

pi_female_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi_sub, aes(x = cf, y = pi), color = cbPalette[8], alpha = 0.3, size=10) + 
  geom_line(data=chr1_data_pi_sub, aes(x=cf, y=rollmean(smooth(pi), 50, na.pad=T)), color = cbPalette[3], size=2) +
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cf, y = pi), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cf, y = pi), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cf, y = pi), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cf, y = pi), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cf, y = pi), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cf, y = pi), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cf, y = pi), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cf, y = pi), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cf, y = pi), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cf, y = pi), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cf, y = pi), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cf, y = pi), label="j", color = 1, alpha = 1, size=9) + 
  # geom_smooth(data=chr1_data_pi, aes(x = cf, y = pi), color = cbPalette[1] ) + 

  labs(x = "Female map (cM)", y = pi_lg1_y_title) +
  geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
  geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) 

png(file.path(outpath, 'pi_lg1.png'), width=900, height=650, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 6,6), "null"))))
grid.text("pi LG1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pi_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(pi_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
dev.off()


## log Y over X plot
YovX_lg1_y_title <- expression(paste("logFC (Y/X allele)"))
chr1_data_pi_sub2 <- subset(chr1_data_pi, chr1_data_pi$mean_m.Y_over_X!='NA') 

YovX_male_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cm, y = log(mean_m.Y_over_X)), color = cbPalette[3], alpha = 0.3, size=10) + 
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cm, y = log(mean_m.Y_over_X)), label="j", color = 1, alpha = 1, size=9) + 
  geom_smooth(data=chr1_data_pi, aes(x = cm, y = log(mean_m.Y_over_X)), color = cbPalette[6] ) + 
  labs(x = "Male map (cM)", y = YovX_lg1_y_title) +
  geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) 

YovX_female_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cf, y = log(mean_m.Y_over_X)), color = cbPalette[8], alpha = 0.3, size=10) + 
  geom_line(data=chr1_data_pi_sub2, aes(x=cf, y=rollmean(smooth(log(mean_m.Y_over_X)), 50, na.pad=T)), color = cbPalette[3], size=2) +
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cf, y = log(mean_m.Y_over_X)), label="j", color = 1, alpha = 1, size=9) + 
#  geom_smooth(data=chr1_data_pi, aes(x = cf, y = log(mean_m.Y_over_X)), color = cbPalette[1] ) + 
  labs(x = "Female map (cM)", y = YovX_lg1_y_title) +
  geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
  geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) 

png(file.path(outpath, 'YovX_lg1.png'), width=900, height=650, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 6,6), "null"))))
grid.text("YovX LG1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(YovX_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(YovX_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
dev.off()

## xy dnds plot
xy_dnds_male_lg1_y_title <- expression(paste("dN/dS (X vs Y)"))
chr1_data_pi_sub3 <- subset(chr1_data_pi, chr1_data_pi$xy_dnds!='NA') 

xy_dnds_male_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cm, y = xy_dnds), color = cbPalette[3], alpha = 0.3, size=10) + 
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cm, y = xy_dnds), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cm, y = xy_dnds), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cm, y = xy_dnds), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cm, y = xy_dnds), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cm, y = xy_dnds), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cm, y = xy_dnds), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cm, y = xy_dnds), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cm, y = xy_dnds), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cm, y = xy_dnds), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cm, y = xy_dnds), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cm, y = xy_dnds), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cm, y = xy_dnds), label="j", color = 1, alpha = 1, size=9) + 
  geom_smooth(data=chr1_data_pi, aes(x = cm, y = xy_dnds), color = cbPalette[6] ) + 
  labs(x = "Male map (cM)", y = xy_dnds_male_lg1_y_title) +
  geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) 

xy_dnds_female_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cf, y = xy_dnds), color = cbPalette[8], alpha = 0.3, size=10) + 
  geom_line(data=chr1_data_pi_sub3, aes(x=cf, y=rollmean(smooth(xy_dnds), 50, na.pad=T)), color = cbPalette[3], size=2) +
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cf, y = xy_dnds), label="a", color = 1, alpha = 1, size=9) +
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cf, y = xy_dnds), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cf, y = xy_dnds), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cf, y = xy_dnds), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cf, y = xy_dnds), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cf, y = xy_dnds), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cf, y = xy_dnds), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cf, y = xy_dnds), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cf, y = xy_dnds), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cf, y = xy_dnds), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cf, y = xy_dnds), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cf, y = xy_dnds), label="j", color = 1, alpha = 1, size=9) + 
#  geom_smooth(data=chr1_data_pi, aes(x = cf, y = xy_dnds), color = cbPalette[1] ) + 
  labs(x = "Female map (cM)", y = xy_dnds_male_lg1_y_title) +
  geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
  geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) 

png(file.path(outpath, 'xy_dnds.png'), width=900, height=650, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 6,6), "null"))))
grid.text("xy_dnds LG1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(xy_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(xy_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
dev.off()

## M annua Ricinus dnds plot
an_ric_dnds_male_lg1_y_title <- expression(paste("dN/dS", italic("M. annua"), " vs ", italic("R. communis")))
chr1_data_pi_sub4 <- subset(chr1_data_pi, chr1_data_pi$an_ric_dnds!='NA') 

an_ric_dnds_male_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cm, y = an_ric_dnds), color = cbPalette[3], alpha = 0.3, size=10) + 
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cm, y = an_ric_dnds), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cm, y = an_ric_dnds), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cm, y = an_ric_dnds), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cm, y = an_ric_dnds), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cm, y = an_ric_dnds), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cm, y = an_ric_dnds), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cm, y = an_ric_dnds), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cm, y = an_ric_dnds), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cm, y = an_ric_dnds), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cm, y = an_ric_dnds), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cm, y = an_ric_dnds), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cm, y = an_ric_dnds), label="j", color = 1, alpha = 1, size=9) + 
  geom_smooth(data=chr1_data_pi, aes(x = cm, y = an_ric_dnds), color = cbPalette[6] ) + 
  labs(x = "Male map (cM)", y = an_ric_dnds_male_lg1_y_title) +
  geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) 

an_ric_dnds_female_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cf, y = an_ric_dnds), color = cbPalette[8], alpha = 0.3, size=10) + 
  geom_line(data=chr1_data_pi_sub4, aes(x=cf, y=rollmean(smooth(an_ric_dnds), 50, na.pad=T)), color = cbPalette[3], size=2) +
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cf, y = an_ric_dnds), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cf, y = an_ric_dnds), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cf, y = an_ric_dnds), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cf, y = an_ric_dnds), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cf, y = an_ric_dnds), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cf, y = an_ric_dnds), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cf, y = an_ric_dnds), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cf, y = an_ric_dnds), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cf, y = an_ric_dnds), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cf, y = an_ric_dnds), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cf, y = an_ric_dnds), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cf, y = an_ric_dnds), label="j", color = 1, alpha = 1, size=9) + 
#  geom_smooth(data=chr1_data_pi, aes(x = cf, y = an_ric_dnds), color = cbPalette[1] ) + 
  labs(x = "Female map (cM)", y = an_ric_dnds_male_lg1_y_title) +
  geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
  geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) 

png(file.path(outpath, 'an_ric_dnds.png'), width=900, height=650, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 6,6), "null"))))
grid.text("an_ric_dnds LG1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(an_ric_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(an_ric_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
dev.off()

## M annua M huetii dnds plot
an_hu_dnds_male_lg1_y_title <- expression(paste("dN/dS", italic("M. annua"), " - ", italic("M. huetii")))
chr1_data_pi_sub5 <- subset(chr1_data_pi, chr1_data_pi$an_hu_dnds!='NA') 

an_hu_dnds_male_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cm, y = an_hu_dnds), color = cbPalette[3], alpha = 0.3, size=10) + 
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cm, y = an_hu_dnds), label="a", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cm, y = an_hu_dnds), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cm, y = an_hu_dnds), label="f",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cm, y = an_hu_dnds), label="p", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cm, y = an_hu_dnds), label="s", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cm, y = an_hu_dnds), label="k",color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cm, y = an_hu_dnds), label="e", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cm, y = an_hu_dnds), label="w", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cm, y = an_hu_dnds), label="r", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cm, y = an_hu_dnds), label="d", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cm, y = an_hu_dnds), label="h", color = 1, alpha = 1, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cm, y = an_hu_dnds), label="j", color = 1, alpha = 1, size=9) + 
  geom_smooth(data=chr1_data_pi, aes(x = cm, y = an_hu_dnds), color = cbPalette[6] ) + 
  labs(x = "Male map (cM)", y = an_hu_dnds_male_lg1_y_title) +
  geom_vline(xintercept=53.85, linetype=2, colour=cbPalette[7]) 

an_hu_dnds_female_lg1_plot <- ggplot() + 
  geom_point(data=chr1_data_pi, aes(x = cf, y = an_hu_dnds), color = cbPalette[8], alpha = 0.3, size=10) + 
  geom_line(data=chr1_data_pi_sub5, aes(x=cf, y=rollmean(smooth(an_hu_dnds), 50, na.pad=T)), color = cbPalette[3], size=2) +
  geom_text(data=chr1_data_pi_auxin_subset, aes(x = cf, y = an_hu_dnds), label="a", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_cytokinin_subset, aes(x = cf, y = an_hu_dnds), label="c", color = 1, alpha = 1, size=9) +  
  geom_text(data=chr1_data_pi_flower_subset, aes(x = cf, y = an_hu_dnds), label="f",color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_pollen_subset, aes(x = cf, y = an_hu_dnds), label="p", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_stop_subset, aes(x = cf, y = an_hu_dnds), label="s", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_kelch_subset, aes(x = cf, y = an_hu_dnds), label="k",color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_ethylene_subset, aes(x = cf, y = an_hu_dnds), label="e", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_swisnf_subset, aes(x = cf, y = an_hu_dnds), label="w", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_retinoldehydrogenase_subset, aes(x = cf, y = an_hu_dnds), label="r", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_calmodulin_subset, aes(x = cf, y = an_hu_dnds), label="d", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_short.root_subset, aes(x = cf, y = an_hu_dnds), label="h", color = 1, alpha = 1, size=9, size=9) + 
  geom_text(data=chr1_data_pi_jasmonicacid_subset, aes(x = cf, y = an_hu_dnds), label="j", color = 1, alpha = 1, size=9, size=9) + 
#  geom_smooth(data=chr1_data_pi, aes(x = cf, y = an_hu_dnds), color = cbPalette[1] ) + 
  labs(x = "Female map (cM)", y = an_hu_dnds_male_lg1_y_title) +
  geom_vline(xintercept=46.4, linetype=2, colour=cbPalette[4]) +
  geom_vline(xintercept=68.67, linetype=2, colour=cbPalette[4]) 

png(file.path(outpath, 'an_hu_dnds.png'), width=900, height=650, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(1, 6,6), "null"))))
grid.text("an_hu_dnds LG1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(an_hu_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(an_hu_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
dev.off()

## Plot all
png(file.path(outpath, 'all.png'), width=2100, height=2970, units="px", pointsize=12, bg="white")
pushViewport(viewport(layout = grid.layout(6, 2, heights = unit(c(1,1,1,1,1,1), "null"))))
print(logfc_male_lg1_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(logfc_female_lg1_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(YovX_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(YovX_female_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(xy_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(xy_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(pi_male_lg1_plot, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
print(pi_female_lg1_plot, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
print(an_hu_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 5, layout.pos.col = 1))
print(an_hu_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 5, layout.pos.col = 2))
print(an_ric_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
print(an_ric_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
dev.off()

dev.copy(pdf,file.path(outpath,'females_all.pdf'), width=21, height=29.7)
pushViewport(viewport(layout = grid.layout(6, 2, heights = unit(c(1,1,1,1,1,1), "null"))))
print(logfc_female_lg1_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pi_female_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(an_hu_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
print(an_ric_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))
print(YovX_female_lg1_plot, vp = viewport(layout.pos.row = 5, layout.pos.col = 1:2))
print(xy_dnds_female_lg1_plot, vp = viewport(layout.pos.row = 6, layout.pos.col = 1:2))
dev.off()

dev.copy(pdf,file.path(outpath,'males_all.pdf'), width=21, height=29.7)
pushViewport(viewport(layout = grid.layout(6, 2, heights = unit(c(1,1,1,1,1,1), "null"))))
print(logfc_male_lg1_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pi_male_lg1_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
print(an_hu_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
print(an_ric_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))
print(YovX_male_lg1_plot, vp = viewport(layout.pos.row = 5, layout.pos.col = 1:2))
print(xy_dnds_male_lg1_plot, vp = viewport(layout.pos.row = 6, layout.pos.col = 1:2))
dev.off()


# Chisq SL, PAR, LG1 regions for candidate genes
 # no enrichment found

sd_lm <- read.table(file.path(datapath, 'sd_lm_data.txt'), header=T) # This is information on sex-linkage inference from SEX-DETector or LepMap3. A gene might be identified as sex-linked from neither, one or both techniques. The final analysis used sex-linkage as estimated from LepMap3. 
str(sd_lm)
boxdata <- merge(sd_lm, plot_data2, by.x='id', by.y='gene', all=F) # Add other information to sex detector genes only 
boxdata2 <- merge(sd_lm, plot_data2, by.x='id', by.y='gene', all=T) # Add sex detector genes to all other information
summary(boxdata$sd_lm)

boxdata2_lg1 <- subset(boxdata2, boxdata2$lg=='1') 

candidate.genes <- summary(grepl(("flower|auxin|cytokinin|pollen|kelch|ethylene|swi|retinol|calmodulin|short-root|jasmonic"), boxdata2$gene.y))
candidate.genes.lg1 <- summary(grepl(("flower|auxin|cytokinin|pollen|kelch|ethylene|swi|retinol|calmodulin|short-root|jasmonic"), boxdata2_lg1$gene.y))
all.genes <- nrow(boxdata2)
all.genes.lg1 <- nrow(boxdata2_lg1)

summary(boxdata2_lg1$cm=="53.85") # PAR vs SL ORFs

boxdata2$islg1 <- boxdata2$lg=='1'
boxdata2$iscg <- grepl(("flower|auxin|cytokinin|pollen|kelch|ethylene|swi|retinol|calmodulin|short-root|jasmonic"), boxdata2$gene.y)
boxdata2$islg1_issl <- boxdata2$lg=='1' & boxdata2$cm=='53.85'

boxdata3 <- subset(boxdata2, boxdata2$lg!='NA') # Mapped data only by LepMap3 (SDonly genes missing) 
boxdata3$sd_lm <- as.character(boxdata3$sd_lm)
boxdata3$sd_lm[boxdata3$sd_lm!="lmonly" & boxdata3$sd=="Au"] <- "Au"
boxdata3$sd_lm[is.na(boxdata3$sd_lm)] <- "Au"
boxdata3$sd_lm <- as.factor(boxdata3$sd_lm)

sdonly <- subset(boxdata2, boxdata2$sd_lm=="sdonly")
sdonly2 <- subset(sdonly, boxdata$sd!="NA") 

boxdata4 <- as.data.frame(rbind(boxdata3, sdonly2)) # Mapped data including SD only genes
boxdata4$isAu <- boxdata4$sd_lm=="Au"
boxdata4$isFDR <- boxdata4$fdr < 0.05
boxdata4$isFDRlogfc <- boxdata4$fdr < 0.05 & abs(boxdata4$logfc > 1)
boxdata4$isMB <- boxdata4$logfc > 0 & boxdata4$fdr < 0.05 # needs to have unbiased genes removed from dataset
boxdata4$isMBlogfc <- boxdata4$logfc > 1 & boxdata4$fdr < 0.05

boxdata4$Au <- "Au"
boxdata4$Au[boxdata4$lg==1] <- "PAR"
 # boxdata4$Au[boxdata4$sd_lm=="both" | boxdata4$sd_lm=="both_somerec" | boxdata4$sd_lm=="lmonly" & boxdata4$lg==1] <- "PAR"
boxdata4$Au[boxdata4$cm=="53.85"] <- "SL"
boxdata4$Au2 <- boxdata4$Au
 # boxdata4$Au2[boxdata4$sd_lm=="sdonly"] <- "PAR" # or SL? # adding these genes does not qualitatively change the result
 # boxdata4$islg1[boxdata4$sd_lm=="sdonly"] <- "TRUE" # same
summary(factor(subset(boxdata4, boxdata4$sd_lm=="sdonly")$logfc))
boxdata4$SB <- "anbiased"
boxdata4$SB[boxdata4$fdr < 0.05 & boxdata4$logfc > 0 ] <- "male"
boxdata4$SB[boxdata4$fdr < 0.05 & boxdata4$logfc < 0 ] <- "female"

write.table(boxdata4, file=file.path(outpath, 'boxdata4_out.txt'), quote=F, row.names=F, sep='\t') # All merged information for other use 

box1vsAu <- subset(boxdata4, boxdata4$lg!='NA')
box1vsAu_FDR <- subset(box1vsAu, box1vsAu$fdr < 0.05)
box1vsAu_FDR_logfc <- subset(box1vsAu, box1vsAu$fdr < 0.05 & abs(box1vsAu$logfc) > 1)

boxSLvsPAR <- subset(box1vsAu, box1vsAu$Au!="Au")
boxSLvsPAR_FDR <- subset(boxSLvsPAR, boxSLvsPAR$fdr < 0.05)
boxSLvsPAR_FDR_logfc <- subset(boxSLvsPAR, boxSLvsPAR$fdr < 0.05 & abs(boxSLvsPAR$logfc) > 1)

boxSLvsAu <- subset(box1vsAu, box1vsAu$Au!="PAR")
boxSLvsAu_FDR <- subset(boxSLvsAu, boxSLvsAu$fdr < 0.05)
boxSLvsAu_FDR_logfc <- subset(boxSLvsAu, boxSLvsAu$fdr < 0.05 & abs(boxSLvsAu$logfc) > 1)

sjt.xtab(box1vsAu$islg1, box1vsAu$iscg, title='Chisq CG LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_cg_LG1_vs_Au.html'))

sjt.xtab(box1vsAu$islg1, box1vsAu$isFDR, title='Chisq sex bias (FDR) LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_SB_FDR_lg1vsau.html'))

sjt.xtab(box1vsAu$islg1, box1vsAu$isFDRlogfc, title='Chisq sex bias (FDR + logfc) LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_SB_FDRlogfc_lg1vsau.html'))

sjt.xtab(box1vsAu_FDR$islg1, box1vsAu_FDR$isMB, title='Chisq male vs female bias (FDR) LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_MFB_FDR_lg1vsau.html'))

sjt.xtab(box1vsAu_FDR_logfc$islg1, box1vsAu_FDR_logfc$isMB, title='Chisq male vs female bias (FDR + logfc) LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_MFB_FDRlogfc_lg1vsau.html'))



sjt.xtab(boxSLvsAu$Au, boxSLvsAu$iscg, title='Chisq CG LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_cg_SL_vs_Au.html'))

sjt.xtab(boxSLvsAu$Au, boxSLvsAu$isFDR, title='Chisq sex bias (FDR) SL vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_SB_FDR_SLvsAu.html'))

sjt.xtab(boxSLvsAu$Au, boxSLvsAu$isFDRlogfc, title='Chisq sex bias (FDR + logfc) SL vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_SB_FDRlogfc_SLvsAu.html'))

sjt.xtab(boxSLvsAu_FDR$Au, boxSLvsAu_FDR$isMB, title='Chisq male vs female bias (FDR) SL vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_MFB_FDR_SLvsAu.html'))

sjt.xtab(boxSLvsAu_FDR_logfc$Au, boxSLvsAu_FDR_logfc$isMB, title='Chisq male vs female bias (FDR + logfc) SL vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_MFB_FDRlogfc_SLvsAu.html'))


sjt.xtab(boxSLvsPAR$Au, boxSLvsPAR$iscg, title='Chisq CG LG1 vs Autosomes', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_cg_SL_vs_PAR.html'))

sjt.xtab(boxSLvsPAR$Au, boxSLvsPAR$isFDR, title='Chisq sex bias (FDR) SL vs PAR', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_SB_FDR_SLvsPAR.html'))

sjt.xtab(boxSLvsPAR$Au, boxSLvsPAR$isFDRlogfc, title='Chisq sex bias (FDR + logfc) SL vs PAR', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_SB_FDRlogfc_SLvsPAR.html'))

sjt.xtab(boxSLvsPAR_FDR$Au, boxSLvsPAR_FDR$isMB, title='Chisq male vs female bias (FDR) SL vs PAR', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_MFB_FDR_SLvsPAR.html'))

 # too few to work
 sjt.xtab(boxSLvsPAR_FDR_logfc$Au, boxSLvsPAR_FDR_logfc$isMB, title='Chisq male vs female bias (FDR + logfc) SL vs PAR', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'chisq_MFB_FDRlogfc_SLvsPAR.html'))


# box cox transformations
 # https://stackoverflow.com/questions/33999512/how-to-use-the-box-cox-power-transformation-in-r

### Without SD only
par(mfrow=c(2,3)) 
par(mar=c(5,5,4,3))

box1vsAu$SB <- factor(box1vsAu$SB)
box1vsAu$Au <- factor(box1vsAu$Au)

# ridge plots

ridge_logfc <- ggplot(subset(box1vsAu, !is.na(logfc)), aes(x=abs(logfc), y=Au:SB, fill = Au:SB)) + geom_density_ridges(aes(point_fill = SB), alpha = .6) + xlab("male/female |logFC|") + scale_fill_manual(values=c("gray", "red", "#0072B2", "gray", "red", "#0072B2", "gray", "red", "#0072B2")) + theme(legend.position="none") + coord_cartesian(ylim = c(1, 10)) +
theme(axis.text.x = element_text(size=26), text = element_text(size=34)) +
annotate(geom="text", x=4.6, y=1.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="anbiased")), size=10, hjust=0) + 
annotate(geom="text", x=4.6, y=2.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=3.3, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=4.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=5.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=6.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=7.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=8.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=4.6, y=9.25, label=nrow(subset(box1vsAu, box1vsAu$logfc!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="male")), size=10, hjust=0)

dev.copy(pdf, file.path(outpath,'ridge_logfc.pdf'), width=12, height=10)
ridge_logfc
dev.off()

ridge_pi <- ggplot(subset(box1vsAu, !is.na(pi)), aes(x=pi, y=Au:SB, fill = Au:SB)) + geom_density_ridges(aes(point_fill = SB), alpha = .6) + xlab("Pi") + scale_fill_manual(values=c("gray", "red", "#0072B2", "gray", "red", "#0072B2", "gray", "red", "#0072B2")) + theme(legend.position="none") + coord_cartesian(ylim = c(1, 10), xlim = c(-0.005, 0.06)) +
theme(axis.text.x = element_text(size=26), text = element_text(size=34)) +
annotate(geom="text", x=.054, y=1.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="anbiased")), size=10, hjust=0) + 
annotate(geom="text", x=.054, y=2.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=3.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=4.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=5.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=6.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=7.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=8.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.054, y=9.25, label=nrow(subset(box1vsAu, box1vsAu$pi!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="male")), size=10, hjust=0)
dev.copy(pdf, file.path(outpath,'ridge_pi.pdf'), width=12, height=10)
ridge_pi
dev.off()


ridge_an_hu_dnds <- ggplot(subset(box1vsAu, !is.na(an_hu_dnds)), aes(x=an_hu_dnds, y=Au:SB, fill = Au:SB)) + geom_density_ridges(aes(point_fill = SB), alpha = .6) + xlab("dN/dS (M. huetii)") + scale_fill_manual(values=c("gray", "red", "#0072B2", "gray", "red", "#0072B2", "gray", "red", "#0072B2")) + theme(legend.position="none") + coord_cartesian(ylim = c(1, 11),  xlim = c(-0.1, .8)) +
theme(axis.text.x = element_text(size=26), text = element_text(size=34)) +
annotate(geom="text", x=.74, y=1.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="anbiased")), size=10, hjust=0) + 
annotate(geom="text", x=.74, y=2.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=3.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=4.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=5.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=6.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=7.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=8.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.74, y=9.25, label=nrow(subset(box1vsAu, box1vsAu$an_hu_dnds!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="male")), size=10, hjust=0)
dev.copy(pdf, file.path(outpath,'ridge_an_hu_dnds.pdf'), width=12, height=10)
ridge_an_hu_dnds
dev.off()

ridge_an_ric_dnds <- ggplot(subset(box1vsAu, !is.na(an_ric_dnds)), aes(x=an_ric_dnds, y=Au:SB, fill = Au:SB)) + geom_density_ridges(aes(point_fill = SB), alpha = .6) + xlab("dN/dS (R. communis)") + scale_fill_manual(values=c("gray", "red", "#0072B2", "gray", "red", "#0072B2", "gray", "red", "#0072B2")) + theme(legend.position="none") + coord_cartesian(ylim = c(1, 10.5),  xlim = c(-0.1, .8)) + 
theme(axis.text.x = element_text(size=26), text = element_text(size=34)) +
annotate(geom="text", x=.7, y=1.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="anbiased")), size=10, hjust=0) + 
annotate(geom="text", x=.7, y=2.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=3.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="Au" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=4.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=5.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=6.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="PAR" & box1vsAu$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=7.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=8.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=.7, y=9.25, label=nrow(subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' & box1vsAu$Au=="SL" & box1vsAu$SB=="male")), size=10, hjust=0)
dev.copy(pdf, file.path(outpath,'ridge_an_ric_dnds.pdf'), width=12, height=10)
ridge_an_ric_dnds
dev.off()


box1vsAu_lg1 <- subset(box1vsAu, box1vsAu$lg==1)

ridge_YoverX <- ggplot(subset(box1vsAu_lg1, !is.na(mean_m.Y_over_X)), aes(x=log(mean_m.Y_over_X), y=Au:SB, fill = Au:SB)) + geom_density_ridges(aes(point_fill = SB), alpha = .6) + xlab("log Y/X allele logFC") + scale_fill_manual(values=c("gray",  "gray", "red", "#0072B2")) + theme(legend.position="none") + coord_cartesian(ylim = c(1, 8.6)) +
theme(axis.text.x = element_text(size=26), text = element_text(size=34)) +
annotate(geom="text", x=3.9, y=1.3, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA' & box1vsAu_lg1$Au=="PAR" & box1vsAu_lg1$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=3.9, y=2.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA' & box1vsAu_lg1$Au=="PAR" & box1vsAu_lg1$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=3.9, y=3.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA' & box1vsAu_lg1$Au=="PAR" & box1vsAu_lg1$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=3.9, y=4.3, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA' & box1vsAu_lg1$Au=="SL" & box1vsAu_lg1$SB=="anbiased")), size=10, hjust=0) + 
annotate(geom="text", x=3.9, y=5.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA' & box1vsAu_lg1$Au=="SL" & box1vsAu_lg1$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=3.9, y=6.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA' & box1vsAu_lg1$Au=="SL" & box1vsAu_lg1$SB=="male")), size=10, hjust=0)
dev.copy(pdf, file.path(outpath,'ridge_YoverX.pdf'), width=12, height=8)
ridge_YoverX
dev.off()


ridge_xy_dnds <- ggplot(subset(box1vsAu_lg1, !is.na(xy_dnds)), aes(x=xy_dnds, y=Au:SB, fill = Au:SB)) + geom_density_ridges(aes(point_fill = SB), alpha = .6) + xlab("X - Y dN/dS") + scale_fill_manual(values=c("gray",  "gray", "red", "#0072B2")) + theme(legend.position="none") + coord_cartesian(ylim = c(1, 8.6), xlim=c(-0.14,1.6)) +
theme(axis.text.x = element_text(size=26), text = element_text(size=34)) +
annotate(geom="text", x=1.5, y=1.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA' & box1vsAu_lg1$Au=="PAR" & box1vsAu_lg1$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=1.5, y=2.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA' & box1vsAu_lg1$Au=="PAR" & box1vsAu_lg1$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=1.5, y=3.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA' & box1vsAu_lg1$Au=="PAR" & box1vsAu_lg1$SB=="male")), size=10, hjust=0) +
annotate(geom="text", x=1.5, y=4.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA' & box1vsAu_lg1$Au=="SL" & box1vsAu_lg1$SB=="anbiased")), size=10, hjust=0) +
annotate(geom="text", x=1.5, y=5.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA' & box1vsAu_lg1$Au=="SL" & box1vsAu_lg1$SB=="female")), size=10, hjust=0) +
annotate(geom="text", x=1.5, y=6.25, label=nrow(subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA' & box1vsAu_lg1$Au=="SL" & box1vsAu_lg1$SB=="male")), size=10, hjust=0)

dev.copy(pdf, file.path(outpath,'ridge_xy_dnds.pdf'), width=12, height=8)
ridge_xy_dnds
dev.off()

box1vsAu2 <- subset(box1vsAu, box1vsAu$logfc!='NA') # remove NA for boxcox test to work
box1vsAu_sbonly <- subset(box1vsAu, box1vsAu$SB!='anbiased') 

mod_logfc_Au <- lm(box1vsAu2$logfc ~ box1vsAu2$Au * box1vsAu2$SB)
anova(mod_logfc_Au)
mod_logfc_Au_sb <- lm(box1vsAu_sbonly$logfc ~ box1vsAu_sbonly$Au * box1vsAu_sbonly$SB )
anova(mod_logfc_Au_sb)
summary(mod_logfc_Au_sb)

mod_abslogfc_Au <- lm(abs(box1vsAu2$logfc) ~ box1vsAu2$Au * box1vsAu2$SB)
summary(mod_abslogfc_Au)
anova(mod_abslogfc_Au)
box1vsAu2$bclfc <- boxcoxnc(abs(box1vsAu2$logfc), method="ad")$tf.data # NOT normal
mod_abslogfc_Au_bc <- lm(box1vsAu2$bclfc ~ box1vsAu2$Au * box1vsAu2$SB)
summary(mod_abslogfc_Au_bc)
anova(mod_abslogfc_Au_bc)

mod_abslogfc_Au_sb <- lm(abs(box1vsAu_sbonly$logfc) ~ box1vsAu_sbonly$SB * box1vsAu_sbonly$Au)
anova(mod_abslogfc_Au_sb)
summary(mod_abslogfc_Au_sb)
box1vsAu_sbonly$bclfc <- boxcoxnc(abs(box1vsAu_sbonly$logfc), method="ad")$tf.data #  normal
mod_abslogfc_Au_sb_bc <- lm(box1vsAu_sbonly$bclfc ~ box1vsAu_sbonly$SB * box1vsAu_sbonly$Au)
anova(mod_abslogfc_Au_sb_bc)
summary(mod_abslogfc_Au_sb_bc)

box1vsAu_pi <- subset(box1vsAu, box1vsAu$pi!='NA') # if remove pi =0, does not affect outcome 
mod_pi_Au <- lm(box1vsAu_pi$pi ~ box1vsAu_pi$Au * box1vsAu_pi$SB)
anova(mod_pi_Au, test='F')
box1vsAu_pi$bcpi <- boxcoxnc(unname(box1vsAu_pi$pi) + 0.0000000001, method="ad")$tf.data # NOT normal
mod_pi_Au_bc <- lm(box1vsAu_pi$bcpi ~ box1vsAu_pi$Au * box1vsAu_pi$SB)
anova(mod_pi_Au_bc, test='F')
summary(mod_pi_Au_bc)
box1vsAu_sbonly_pi <- subset(box1vsAu_sbonly, box1vsAu_sbonly$pi!='NA' & box1vsAu_sbonly$pi!=0)
mod_pi_Au_sb <- lm(box1vsAu_sbonly_pi$pi ~ box1vsAu_sbonly_pi$SB * box1vsAu_sbonly_pi$Au)
anova(mod_pi_Au_sb, test='F')
summary(mod_pi_Au_sb)
box1vsAu_sbonly_pi$bcpi <- boxcoxnc(unname(box1vsAu_sbonly_pi$pi), method="ad")$tf.data #  normal
mod_pi_Au_sb_pi_bc <- lm(box1vsAu_sbonly_pi$bcpi ~ box1vsAu_sbonly_pi$SB * box1vsAu_sbonly_pi$Au)
anova(mod_pi_Au_sb_pi_bc)
summary(mod_pi_Au_sb_pi_bc)

box1vsAu_an_hu_dnds <- subset(box1vsAu, box1vsAu$an_hu_dnds!='NA') # If remove an_hu_dnds =0, does not affect outcome 
mod_an_hu_dnds_Au <- lm(box1vsAu_an_hu_dnds$an_hu_dnds ~ box1vsAu_an_hu_dnds$Au * box1vsAu_an_hu_dnds$SB)
anova(mod_an_hu_dnds_Au, test='F')
box1vsAu_an_hu_dnds$an_hu_dnds_bc <- boxcoxnc(unname(box1vsAu_an_hu_dnds$an_hu_dnds + 0.0000000001), method="ad")$tf.data #  normal
mod_an_hu_dnds_Au_bc <- lm(box1vsAu_an_hu_dnds$an_hu_dnds_bc ~ box1vsAu_an_hu_dnds$Au * box1vsAu_an_hu_dnds$SB)
anova(mod_an_hu_dnds_Au_bc)
summary(mod_an_hu_dnds_Au_bc)

box1vsAu_sbonly_an_hu_dnds <- subset(box1vsAu_sbonly, box1vsAu_sbonly$an_hu_dnds!='NA')
mod_an_hu_dnds_Au_sb <- lm(box1vsAu_sbonly_an_hu_dnds$an_hu_dnds ~ box1vsAu_sbonly_an_hu_dnds$SB * box1vsAu_sbonly_an_hu_dnds$Au)
anova(mod_an_hu_dnds_Au_sb, test='F')
summary(mod_an_hu_dnds_Au_sb)
box1vsAu_sbonly_an_hu_dnds$an_hu_dnds_bc <- boxcoxnc(unname(box1vsAu_sbonly_an_hu_dnds$an_hu_dnds + 0.0000000001), method="ad")$tf.data #  normal
mod_an_hu_dnds_Au_sb_bc <- lm(box1vsAu_sbonly_an_hu_dnds$an_hu_dnds_bc ~ box1vsAu_sbonly_an_hu_dnds$SB * box1vsAu_sbonly_an_hu_dnds$Au)
anova(mod_an_hu_dnds_Au_sb_bc, test='F')
summary(mod_an_hu_dnds_Au_sb_bc)

box1vsAu_an_ric_dnds <- subset(box1vsAu, box1vsAu$an_ric_dnds!='NA' )
mod_an_ric_dnds_Au <- lm(box1vsAu_an_ric_dnds$an_ric_dnds ~ box1vsAu_an_ric_dnds$Au * box1vsAu_an_ric_dnds$SB)
anova(mod_an_ric_dnds_Au, test='F')
box1vsAu_an_ric_dnds$an_ric_dnds_bc <- boxcoxnc(unname(box1vsAu_an_ric_dnds$an_ric_dnds  + 0.000000000001), method="ad")$tf.data #  normal
mod_an_ric_dnds_Au_bc <- lm(box1vsAu_an_ric_dnds$an_ric_dnds_bc ~ box1vsAu_an_ric_dnds$Au * box1vsAu_an_ric_dnds$SB)
anova(mod_an_ric_dnds_Au_bc, test='F')
summary(mod_an_ric_dnds_Au_bc, test='F')

box1vsAu_sbonly_an_ric_dnds <- subset(box1vsAu_sbonly, box1vsAu_sbonly$an_ric_dnds!='NA' )
mod_an_ric_dnds_Au_sb <- lm(box1vsAu_sbonly_an_ric_dnds$an_ric_dnds ~ box1vsAu_sbonly_an_ric_dnds$SB * box1vsAu_sbonly_an_ric_dnds$Au)
anova(mod_an_ric_dnds_Au_sb, test='F')
summary(mod_an_ric_dnds_Au_sb)
box1vsAu_sbonly_an_ric_dnds$an_ric_dnds_bc <- boxcoxnc(unname(box1vsAu_sbonly_an_ric_dnds$an_ric_dnds), method="ad")$tf.data #  normal
mod_an_ric_dnds_Au_sb_bc <- lm(box1vsAu_sbonly_an_ric_dnds$an_ric_dnds_bc ~ box1vsAu_sbonly_an_ric_dnds$SB * box1vsAu_sbonly_an_ric_dnds$Au)
anova(mod_an_ric_dnds_Au_sb_bc, test='F')
summary(mod_an_ric_dnds_Au_sb_bc)


box1vsAu_lg1 <- subset(box1vsAu, box1vsAu$lg==1)
box1vsAu_lg1_YovX <- subset(box1vsAu_lg1, box1vsAu_lg1$mean_m.Y_over_X!='NA')

box1vsAu_lg1_sbonly <- subset(box1vsAu_lg1, box1vsAu_lg1$SB!='anbiased') 
box1vsAu_lg1_sbonly_YovX <- subset(box1vsAu_lg1_sbonly, box1vsAu_lg1_sbonly$mean_m.Y_over_X!='NA')

mod_xylogfc <- lm(box1vsAu_lg1_YovX$mean_m.Y_over_X ~ box1vsAu_lg1_YovX$SB * box1vsAu_lg1_YovX$Au)
anova(mod_xylogfc, test='F')
summary(mod_xylogfc)
box1vsAu_lg1_YovX$mean_m.Y_over_X_bc <- boxcoxnc(unname(box1vsAu_lg1_YovX$mean_m.Y_over_X), method="ad")$tf.data # NOT normal
mod_xylogfc_bc <- lm(box1vsAu_lg1_YovX$mean_m.Y_over_X_bc ~ box1vsAu_lg1_YovX$SB * box1vsAu_lg1_YovX$Au)
anova(mod_xylogfc_bc, test='F')
summary(mod_xylogfc_bc)

mod_xylogfc_sb <- lm(box1vsAu_lg1_sbonly_YovX$mean_m.Y_over_X ~ box1vsAu_lg1_sbonly_YovX$SB * box1vsAu_lg1_sbonly_YovX$Au)
anova(mod_xylogfc_sb, test='F')
summary(mod_xylogfc_sb)
box1vsAu_lg1_sbonly_YovX$mean_m.Y_over_X_bc <- boxcoxnc(unname(box1vsAu_lg1_sbonly_YovX$mean_m.Y_over_X), method="ad")$tf.data #  normal
mod_xylogfc_sb_bc <- lm(box1vsAu_lg1_sbonly_YovX$mean_m.Y_over_X_bc ~ box1vsAu_lg1_sbonly_YovX$SB * box1vsAu_lg1_sbonly_YovX$Au)
anova(mod_xylogfc_sb_bc, test='F')
summary(mod_xylogfc_sb_bc)

box1vsAu_lg1_xydnds <- subset(box1vsAu_lg1, box1vsAu_lg1$xy_dnds!='NA')
box1vsAu_lg1_sbonly_xydnds <- subset(box1vsAu_lg1_sbonly, box1vsAu_lg1_sbonly$xy_dnds!='NA')
mod_xydnds <- lm(box1vsAu_lg1_xydnds$xy_dnds ~ box1vsAu_lg1_xydnds$SB * box1vsAu_lg1_xydnds$Au)
anova(mod_xydnds, test='F')
summary(mod_xydnds)
box1vsAu_lg1_xydnds$xy_dnds_bc <- boxcoxnc(unname(box1vsAu_lg1_xydnds$xy_dnds + 0.0000000000001), method="ad")$tf.data # NOT normal
mod_xydnds_bc <- lm(box1vsAu_lg1_xydnds$xy_dnds_bc ~ box1vsAu_lg1_xydnds$SB * box1vsAu_lg1_xydnds$Au)
anova(mod_xydnds_bc, test='F')
summary(mod_xydnds_bc)

mod_xydnds_sb <- lm(box1vsAu_lg1_sbonly_xydnds$xy_dnds ~ box1vsAu_lg1_sbonly_xydnds$SB * box1vsAu_lg1_sbonly_xydnds$Au)
anova(mod_xydnds_sb, test='F')
summary(mod_xydnds_sb)
box1vsAu_lg1_sbonly_xydnds$xy_dnds_bc <- boxcoxnc(unname(box1vsAu_lg1_sbonly_xydnds$xy_dnds + 0.0000000000001), method="ad")$tf.data #  normal
mod_xydnds_sb_bc <- lm(box1vsAu_lg1_sbonly_xydnds$xy_dnds_bc ~ box1vsAu_lg1_sbonly_xydnds$SB * box1vsAu_lg1_sbonly_xydnds$Au)
anova(mod_xydnds_sb_bc, test='F')
summary(mod_xydnds_sb_bc)

tab_model(mod_abslogfc_Au_bc, mod_pi_Au_bc, mod_an_hu_dnds_Au_bc, mod_an_ric_dnds_Au_bc, file=file.path(outpath, 'bclm_inclAu.html'), show.aic=T)

tab_model(mod_abslogfc_Au_sb_bc, mod_pi_Au_sb_pi_bc, mod_an_hu_dnds_Au_sb_bc, mod_an_ric_dnds_Au_sb_bc, show.aic=T, file=file.path(outpath, 'bclm_inclAu_sbonly.html'))

tab_model(mod_xylogfc_bc, mod_xydnds_bc, show.aic=T, file=file.path(outpath, 'bclm_exclAu.html'))

tab_model(mod_xylogfc_sb_bc, mod_xydnds_sb_bc, show.aic=T, file=file.path(outpath, 'bclm_exclAu_sbonly.html'))

# Plot sex bias vs allele specific expression
 # one ethylene is slightly Y biased and male biased, but it is mapped to the wrong LG (2).
par(mar=c(5,5,4,3))
plot(c(-5,5), c(-1,2.6), type="n", xlab="logFC (Y/X)", ylab="logFC (male/female)", main="Sex bias vs allele specific expression", cex.main=1.8, cex.lab=1.3)
points(log(boxdata4$mean_m.Y_over_X), boxdata4$logfc, pch=1,  col=1)
dev.copy(pdf, file.path(outpath,'logfc_vs_xylogfc.pdf'), width=10, height=10)
dev.off()


