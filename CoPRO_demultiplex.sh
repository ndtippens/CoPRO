#!/bin/bash
# USAGE:
# Provide a list of fastq or fastq.gz files to be aligned. Only specify one of each read pair; the other filename is assumed (R2 if R1, or R1 if R2.)
# Example:
# ./CoPRO_demux.sh ./data/K562_CoPRO*R1.fastq.gz

# number of processors
nproc=8
# path to bowtie2 genome index (don't include .1.bt2)
genome="~/bin/genomes/CoPRO/CoPRO_hg19_dm6"


RNA5p="CTGTCTCTTATACACATCTCCGAGCCCACGAGACAT"
RNA3p="GATCGTCGGACTGTAGAACTCTGAACGTGTAG"

# sample barcodes are in 5' of Read1
BC_ADAPTORS="\
	-g BC1=CGTGATC \
	-g BC2=ACATCGC \
	-g BC3=GCCTAAC \
	-g BC4=TGGTCAC \
	-g BC5=CACTGTC \
	-g BC6=ATTGGCC \
	-g BC8=TCAAGTC \
	-g BC9=CTGATCC \
	-g BC10=AAGCTAC \
	-g BC17=CTCTACC \
	-g BC19=ACTAGC"

for fastq in "$@"
do
	if [ -z `basename $fastq | grep -i .fastq` ]
	then
		echo `basename $fastq` "does not have .fastq suffix - aborting";
		exit 1;
	fi
done


for fastq in "$@"
do
	dname=`dirname $fastq`;
	fname=`basename $fastq`;
	fpath=$dname/${fname%_R[12].fastq*};
	ext=${fname#*_R[12].};

	echo "**************   BEGIN $fname   **************"

	if [ ! -e ${fpath}_BC*_R1.fastq.gz ]
	then
		# -O 5 = require 5 matching bases to demultiplex
		cutadapt --no-trim $BC_ADAPTORS -O 5 \
		--untrimmed-o ${fpath}_noBC_R1.fastq.gz --untrimmed-p ${fpath}_noBC_R2.fastq.gz \
		-o ${fpath}_{name}_R1.fastq.gz -p ${fpath}_{name}_R2.fastq.gz \
		${fpath}_R1.$ext ${fpath}_R2.$ext;
	fi  

	echo "************** FINISHED $fname **************"
done

