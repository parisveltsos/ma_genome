anno <- read.table('~/git/ma_genome/input/de_analysis/ma_annotation.txt', header=T)
str(anno)
de <- read.table('~/git/ma_genome/output/de_analysis/de_0.05_0_ma.txt', header=T)
str(de)
merged1 <- merge(anno, de, by.x="ID", by.y="gid", all=F )
write.table(merged1, file='~/git/ma_genome/output/de_analysis/annotated_de.txt', quote=F, row.names=F, sep="\t")
