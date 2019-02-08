#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o merc2-output.txr
#BSUB -e merc2-error.txt
#BSUB -J merc2
#BSUB -u parisveltsos@gmail.com
#BSUB -N
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -R "rusage[mem=20000]"
#BSUB -M 20194304

module add UHTS/Analysis/vcftools/0.1.14;
module add UHTS/Analysis/GenomeAnalysisTK/3.7;
module add Development/java/1.8.0_121;
module add UHTS/Analysis/picard-tools/2.2.1;
module add UHTS/Aligner/bwa/0.7.13;
module add UHTS/Analysis/samtools/1.3;
module add R/3.3.2;


REFERENCE=/scratch/beegfs/monthly/pveltsos/Mercurialis/reference/SDtranscriptome.fasta
DICTIONARY=/scratch/beegfs/monthly/pveltsos/Mercurialis/reference/SDtranscriptome.fasta.dict
READS1=/scratch/beegfs/monthly/pveltsos/Mercurialis/Parents/500G1_WTCHG_25562_04_1.fq.gz
READS2=/scratch/beegfs/monthly/pveltsos/Mercurialis/Parents/500G1_WTCHG_25562_04_2.fq.gz
NAME=G1

cd /scratch/cluster/monthly/pveltsos/Mercurialis/

GenomeAnalysisTK -T GenotypeGVCFs -R $REFERENCE --variant G1/G1.realigned.snps.indels.g.vcf --variant G2/G2.realigned.snps.indels.g.vcf --variant G4/G4.realigned.snps.indels.g.vcf --variant M1/M1.realigned.snps.indels.g.vcf --variant M2/M2.realigned.snps.indels.g.vcf --variant M3/M3.realigned.snps.indels.g.vcf -o allParents.raw.snps.indels.vcf

echo "vcftools to split to indels and snps, to use for better filtering by GATK. Note that the output filenames are longer than defined in the command"
vcftools --vcf allParents.raw.snps.indels.vcf --keep-only-indels --recode --recode-INFO-all --out allParents.raw.indels
vcftools --vcf allParents.raw.snps.indels.vcf --remove-indels --recode --recode-INFO-all --out allParents.raw.snps

# GenomeAnalysisTK -T SelectVariants -R $REFERENCE -V allParents.raw.snps.indels.vcf -selectType SNP -o raw_snps.vcf

echo "filter out indels and overlapping SNPs"
GenomeAnalysisTK -T VariantFiltration -R $REFERENCE -o allParents.indelFiltered.snp.vcf --mask allParents.raw.indels.recode.vcf --maskExtension 3 --maskName InDel --variant allParents.raw.snps.recode.vcf 

echo "snp quality filtering"

java -jar /software/UHTS/Analysis/GenomeAnalysisTK/3.7/bin/GenomeAnalysisTK.jar -T VariantFiltration -R $REFERENCE -o allParents.fullyFiltered.snp.vcf --variant allParents.indelFiltered.snp.vcf --clusterWindowSize 10 --clusterSize 3 --filterExpression "QUAL <  30" --filterName "QualFilter" --filterExpression "DP <30" --filterName "DepthFilter" --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) >0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "QD < 5.0" --filterName "QualByDepth"

echo "indel quality filtering"
java -jar /software/UHTS/Analysis/GenomeAnalysisTK/3.7/bin/GenomeAnalysisTK.jar -T VariantFiltration -R $REFERENCE -o allParents.fullyFiltered.indels.vcf --variant allParents.raw.indels.recode.vcf --clusterWindowSize 10 --clusterSize 3 --filterExpression "QUAL <  30" --filterName "QualFilter" --filterExpression "DP < 30" --filterName "DepthFilter" --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "QD < 5.0" --filterName "QualByDepth"

echo "extract high quality snps and indels"
cat allParents.fullyFiltered.snps.vcf | grep 'PASS\|^#' > allParents.highQual.snps.vcf

cat allParents.fullyFiltered.indels.vcf | grep 'PASS\|^#' > allParents.highQual.indels.vcf


echo "snp quality filtering - modern"

java -jar /software/UHTS/Analysis/GenomeAnalysisTK/3.7/bin/GenomeAnalysisTK.jar -T VariantFiltration -R $REFERENCE -o allParents.fullyFiltered.modern.snp.vcf --variant allParents.indelFiltered.snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "modern"


allParents.fullyFiltered.modern.snp.vcf | grep 'PASS\|^#' > allParents.highQual.modern.snps.vcf


