nproc=7
genome="~/bin/genomes/CoPRO/CoPRO_hg19_dm6"

for len in 19 20 21 22 23 24
do
	echo $len;
	# add -f to specify FASTA format
	bowtie2 -p $nproc -x $genome -f --very-sensitive \
	-U "../data/${len}mers.fa.gz" |
	# convert to BAM
	samtools view -bS -@ 2 '-' > ../data/CoPRO_hg19_dm6/${len}mers.bam;
done

