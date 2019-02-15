##Move to correct WD
cd /path/to/files

module load repeatmasker/4.0.7

module add SequenceAnalysis/GenePrediction/augustus/3.2.3;

module load blat

# mask genome with Mercurialis repeats

RepeatMasker -lib MercurialislRepeats.fa -gff Mannua_M1_v1.4.fa -dir out

# blat transcriptome to genome

blat -minIdentity=92 Mannua_M1_v1.4.fa.masked SexDETector_orfs.fa -out=psl ORFs.psl

# sort file otherwise blat2hints.pl complains
# plat2hints.pl from https://github.com/nextgenusfs/augustus/blob/master/scripts/blat2hints.pl

cat ORFs.psl | sort -n -k 16,16 | sort -s -k 14,14 > ORFs.sorted.psl

# make hints file

./blat2hints.pl --in=ORFs.sorted.psl --out=hints.E.gff

# run augustus

augustus --species=arabidopsis --hintsfile=hints.E.gff --extrinsicCfgFile=extrinsic.ME.cfg Mannua_M1_v1.4.fa.masked > augustus_predictions.gff


