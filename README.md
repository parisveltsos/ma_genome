# Scripts used for *Mercurialis* genome analysis

## Preliminaries

The scripts should work if you put the analysis files in `~/git/ma_genome/`. If not you will need to change the working directory in the scripts.

The output folder is not synced and you will need to create it

	mkdir `~/git/ma_genome/output/'

Prerequisite packages (they appear in the beginning of R scripts) can be installed with

	install.packages("package_name", dependencies=TRUE)

or for edgeR and topGO

	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install("package_name", version = "3.8")

## Differential expression analysis

Run the script `scripts/de_analysis/edgeR_main.r`. It generates 152 Mb output of the differential expression analysis in `~/git/ma_genome/output/de_analysis` once GO analysis is also run.

### Annotate DE genes

Run script `~/git/ma_genome/scripts/de_analysis/annotation_export.r`.

The output file is `~/git/ma_genome/output/de_analysis/annotated_de.txt`

### GO analysis

The GO analysis is run separately after the main script finishes. Output is produced for all DE genes (`GO_05` folder), upregulated genes (`GO_UP`) and downregulated genes (folder `GO_DOWN`).

	cd ~/git/ma_genome/output/de_analysis

	for i in $(ls | grep GO_0.05); do cd $i; source ~/git/ma_genome/scripts/de_analysis/GO_analysis_05.sh; cd ..; done

This takes a while to run, because there is large numbers of DE genes.

## Molecular evolution analysis

The script `molecularEvolutionAnalysis.r` is documented. Briefly:

lines 1-60: Merge dN,dS data from pairwise comparisons (X-Y, M. annua-M. huetii, M. annua-R. communis, genetic map location in males and females, transcript annotation with GO terms and gene names).

lines 63-98: Investigation of dN/dS correlation with gene expression. No effect was found and the result was not reported.

lines 100-109: Merging of pi data with all other data.

lines 110-137: Subsetting of data to specific linkage groups, and addition of the genetic map distances of all LGs for easier plotting.

lines 139-221: Definition of candidate genes based on grep-ing of their names in the transcriptome annotation.

lines 229-606: Plotting of metrics associated with the transcripts across LG1 for the male and female recombination map. Figure 5 is produced from this.

lines 609-616: Addition of sex-linkage inference from SEX-DETector into all data. This information is not used.

lines 620-705 :Chisq tests for enrichment in sex-biased, male- vs female- biased and candidate genes for pairs of Au, SL, PAR regions. Table 3 is produced from this.

lines 708-933: Box-Cox transformations and glms for the effects of genomic position and sex-bias for transcriptome metrics. Ridge plots (Figure 6) starting at line 712 on non-transformed data to show their distribution. glms are on box-cox transformed data.

lines 935-941: Plot of sex-bias vs allele-specific expression. 

## Stop codon detection

### Install bioperl

[Instructions](https://bioperl.org/INSTALL.html)

on a Mac

	brew install cpanm

	cpanm Bio::Perl

	cpanm Bio::SeqIO

I needed to run the perl Build.PL part of the instructions. Bioperl is needed by the `translate.pl` script, running it on a cluster where BioPerl is already installed might be easier than installing locally.

### Calculate protein length of X and Y sequences from SEX-DETector

	cd ~/git/ma_genome/input/stop_codons
	mkdir temp

Traslate the sequences to protein, and keep the sequence up to the first stop codon

	perl ~/git/ma_genome/scripts/translate.pl f4m6_sex-linked_sequences.fasta | sed 's/\*.*//' > temp/f4m6stop.fa

	perl ~/git/ma_genome/scripts/translate.pl f3m5_sex-linked_sequences.fasta | sed 's/\*.*//' > temp/f3m5stop.fa

Calculate the length of each translated sequence. This is run for the 2 families used by SEX-DETector separately

f4m6:

	perl ~/git/ma_genome/scripts/length_fasta.pl temp/f4m6stop.fa > temp/f4m6_lengths.txt

	cat <(echo -e 'id\tX1length') <(grep X1 temp/f4m6_lengths.txt | perl -pe 's/\>// ; s/\|// ; s/_X1//') > temp/f4m6X1.txt

	cat <(echo -e 'id\tYlength') <(grep Y temp/f4m6_lengths.txt | perl -pe 's/\>// ; s/\|// ; s/_Y//') > temp/f4m6Y.txt

	paste temp/f4m6X1.txt temp/f4m6Y.txt > temp/f4m6_stop.txt

f3m5:

	perl ~/git/ma_genome/scripts/length_fasta.pl temp/f3m5stop.fa > temp/f3m5_lengths.txt

	cat <(echo -e 'id\tX1length') <(grep X1 temp/f3m5_lengths.txt | grep -v 12275 | perl -pe 's/\>// ; s/\|// ; s/_X1//') > temp/f3m5X1.txt

	cat <(echo -e 'id\tYlength') <(grep Y temp/f3m5_lengths.txt | perl -pe 's/\>// ; s/\|// ; s/_Y//') > temp/f3m5Y.txt

	paste temp/f3m5X1.txt temp/f3m5Y.txt > temp/f3m5_stop.txt

### Identify shorter Y protein in R

	f3m5 <- read.table('~/git/ma_genome/input/stop_codons/temp/f3m5_stop.txt', header=T)
	str(f3m5)

	subf3m5 <- subset(f3m5, f3m5$X1length!=f3m5$Ylength) 
	subf3m5

 # comp14079_c0_seq1m.7322
 # comp18849_c0_seq1m.25736

f4m6 <- read.table('~/git/ma_genome/input/stop_codons/temp/f4m6_stop.txt', header=T)
str(f4m6)

subf4m6 <- subset(f4m6, f4m6$X1length!=f4m6$Ylength) 
subf4m6
 
 # comp16827_c0_seq1m.24466 
 # comp17541_c0_seq10m.26313

## Stop codon tree

Script `stop_codon_tree.sh` aligns and constructs a tree based on the sanger sequence of recovered clones (File S2).

## PAML analysis

The script `~/git/ma_genome/scripts/paml.r` is documented and compares dN, dS, dN/dS between X and Y, and tests how many Y linked sequeces are supported to be evolving faster.

## Genome capture pipeline

Scripts for running on a cluster are provided in `~/git/ma_genome/scripts/Genome capture pipeline`

## Lep-Map3

Download Lep-Map3 from its [website](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/). 

### Generate posterior input file

The pipeline at `~/git/ma_genome/scripts/lepmap3/posterior_generation.sh` describes how to make the posterior.gz file that LepMap3 uses as input. The pipeline requires the F2 reads, male genomic reads, and transcriptome.

The posterior.gz used is provided in the `large_intermediate` folder.

### Run LM3

The pipeline at `~/git/ma_genome/scripts/lepmap3/LM3.sh` generates the linkage maps with LM3.

## SNP calling from 6 sequenced individuals

The scripts used are at `~/git/ma_genome/scripts/SNPcalling`. First run the `*bjob*` files and then the `part2.sh` script. It outputs high confidence SNPs and high confidence indels in vcf format.

These vcf files are provided in the `large_intermediate` folder.

## kmer analysis

Generate a `.histo` file on a computer cluster

	module add UHTS/Analysis/jellyfish/2.2.6; gzip -dc M1_1.fa.gz M1_2.fa.gz | jellyfish count -m 31 -o fastqM1_31.counts -C -s 10000000000 -U 500 -t 30 /dev/fd/0 

	module add UHTS/Analysis/jellyfish/2.2.6; jellyfish histo fastqM1_31.counts > fastqM1_31_counts.histo

Upload the resulting file in [Genomescope](http://qb.cshl.edu/genomescope/).

## gff files

`merc_blat.txt` generates a gff on the genome based on the transcriptome.

`augustus mercurialis.sh` generates a gff by training augustus gene prediction based on the *M. annua* transcriptome, after masking the genome using the Mercurialis repeat library.


