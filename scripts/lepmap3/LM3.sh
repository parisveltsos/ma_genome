# replace individual names...

	(cat header.txt;zcat posterior.gz | awk '(NR>1)') | gzip > posterior2.gz

# calculate ibd
	
	zcat posterior2.gz | awk '(NR==1){print;print;print;print;print}(NR==1||rand()<0.02)' | java -cp Lep-MAP3/bin IBD data=- >ibd.txt

# Find suspicious individuals
	grep -w m6 ibd.txt |sort -n -k 3,3|grep -w -f <(awk '($3=="m6" && $4=="f4")' pedigree.txt|cut -f 2)
	#fB_2	m6	0.09	0.819	0.181	0
	#fL_2	m6	0.147	0.707	0.293	0
	#fJ_2	m6	0.182	0.636	0.364	0
	#...

# comment individuals fB_2,fL_2,fJ_2

# create the final pedigree
	
	grep -v "#" pedigree.txt |./transpose_tab|awk '{print "CHR\tPOS\t" $0}' >ped_t.txt

# ParentCall2

	zcat posterior2.gz | java -cp Lep-MAP3/bin ParentCall2 data=ped_t.txt posteriorFile=- halfSibs=1 removeNonInformative=1 | gzip >data.call.gz

# Filtering2

	zcat data.call.gz | java -cp Lep-MAP3/bin Filtering2 dataTolerance=0.001 data=- | gzip >data.filt.gz 

# split data for each transcript to folder tmp

	zcat data.filt.gz|awk '{fn="tmp/" $1; print $0 >fn}'

# Run OrderMarkers2 on each transcript

	for i in $(cat tr.txt)
	do
		awk '{print 1}' tmp/$i >map_tmp.txt
		cat tmp/CHR tmp/$i|java -cp Lep-MAP3/bin OrderMarkers2 data=- map=map_tmp.txt scale=1 1 outputPhasedData=4 recombination1=0 recombination2=0 improveOrder=0 hyperPhaser=1 >tmp/order$i.txt 2>>err
	done

# SeparateChromosomes2

	awk -vpedigree=1 -f order2data.awk tmp/ordercomp* | java -cp Lep-MAP3/bin SeparateChromosomes2 data=- numThreads=8 lodLimit=10 >map10.txt

# JoinSingles2All

	awk -vpedigree=1 -f order2data.awk tmp/ordercomp* | java -cp Lep-MAP3/bin JoinSingles2All data=- numThreads=8 lodLimit=8 map=map10.txt lodDifference=2 iterate=1 >map10_js.tx

# OrderMarkers2

	for i in {1..8}; do awk -vpedigree=1 -f order2data.awk tmp/ordercomp*|java -cp Lep-MAP3/bin OrderMarkers2 data=- map=map10_js.txt chromosome=$i outputPhasedData=1 >order$i.txt 2>order$i.err; done

# create file of SNP names

	for i in tmp/ordercomp*; do awk -vn=$i 'END{if (NR>1) print substr(substr(n,1,length(n)-4), 10)}' $i; done|awk 'BEGIN{print "#"}1' >snps_tr.txt

# create final map
	for i in {1..8}
	do
		awk -vn=$i '{print $1 "\t" n "\t "$2 "\t" $3 "\t" $7 "\t" $10 "\t" $13 "\t" $16 "\t" $19}' order$i.txt
	done|awk -vOFS="\t" '(NR==FNR){m[NR-1]=$1}(NR!=FNR){if ($1 in m) {$1=m[$1];print}}' snps_tr.txt - >final_map.txt

# find the sex region

	grep -P "48.30\t53.85" final_map.txt
