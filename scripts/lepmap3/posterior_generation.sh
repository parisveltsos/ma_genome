	module add UHTS/Aligner/bowtie2/2.3.0;
	module add UHTS/Analysis/samtools/1.3;

	bowtie2-build -f SDtranscriptome.fasta SDtranscriptome.db

	bowtie2 -x SDtranscriptome.db -1 <(gzip -dc 500M1_WTCHG_25562_11_1_sequence.txt.gz) -2 <(gzip -dc 500M1_WTCHG_25562_11_2_sequence.txt.gz) -S 500M1_WTCHG_25562_11.sam -p 4 --no-unal
 
	samtools view -Sb 500M1_WTCHG_25562_11.sam > 500M1_WTCHG_25562_11.bam

	samtools sort 500M1_WTCHG_25562_11.bam 500M1_WTCHG_25562_11_sorted.bam

	samtools index 500M1_WTCHG_25562_11_sorted.bam

	rm 500M1_WTCHG_25562_11.sam

	rm 500M1_WTCHG_25562_11.bam

# Put F2 reads in same folder and tun the bowtie.sh script in parallel for each individual

	for i in $(ls | grep txt.gz | grep _1_ | awk -F '_1_' '{print $1}'); do bsub <<< "/scratch/beegfs/weekly/pveltsos/LMapping/F2/bowtie.sh $i"; done

# make file Lep-Map3 needs

# the awk files are in LepMap3/scripts [LM3 website](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/)

	bsub -q dee-hugemem <<< "module add UHTS/Analysis/samtools/1.3; samtools mpileup -q 10 -Q 10 -s $(cat sorted_bams) | awk -f pileupParser2.awk | awk -f pileup2posterior.awk | gzip > posterior.gz"

	
