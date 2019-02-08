#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o G2xM1_f1-output.txt
#BSUB -e G2xM1_f1-error.txt
#BSUB -J G2xM1_f1
#BSUB -u parisveltsos@gmail.com
#BSUB -N
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10194304

module add UHTS/Aligner/bowtie2/2.3.0;
module add UHTS/Analysis/samtools/1.3;

# READS=G2xM1_f1_WTCHG_25561_07
READS=${1}

cd /scratch/beegfs/weekly/pveltsos/LMapping

# bowtie2-build -f SDtranscriptome.fasta SDtranscriptome.db

bowtie2 -x SDtranscriptome.db -1 <(gzip -dc ${READS}_1_sequence.txt.gz) -2 <(gzip -dc ${READS}_2_sequence.txt.gz) -S ${READS}.sam -p 8 --no-unal

rm ${READS}_1_sequence.txt.gz
rm ${READS}_2_sequence.txt.gz

samtools view -Sb ${READS}.sam > ${READS}.bam

rm ${READS}.sam

samtools sort ${READS}.bam >  ${READS}_sorted.bam

rm ${READS}.bam

samtools index ${READS}_sorted.bam



