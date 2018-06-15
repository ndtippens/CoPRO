#!/bin/bash
# USAGE:
# Provide a list of fastq or fastq.gz files to be aligned. Only specify one of each read pair; the other filename is assumed (R2 if R1, or R1 if R2.)
# Example:
# ./CoPRO_align.sh ./data/K562*R1.fastq.gz

# number of processors
nproc=8
# path to bowtie2 genome index (don't include .*.bt2 extension)
genome="~/bin/genomes/CoPRO/CoPRO_hg19_dm6"


# lowercase a = trim from 3' of R1
# uppercase A = trim from 3' of R2
RNA5p="CTGTCTCTTATACACATCTCCGAGCCCACGAGACAT"
RNA3p="GATCGTCGGACTGTAGAACTCTGAACGTGTAG"
BCs[1]="-a CGTGATC...$RNA5p -A GATCACG$RNA3p"
BCs[2]="-a ACATCGC...$RNA5p -A GCGATGT$RNA3p" 
BCs[3]="-a GCCTAAC...$RNA5p -A GTTAGGC$RNA3p"
BCs[4]="-a TGGTCAC...$RNA5p -A GTGACCA$RNA3p"
BCs[5]="-a CACTGTC...$RNA5p -A GACAGTG$RNA3p"
BCs[6]="-a ATTGGCC...$RNA5p -A GGCCAAT$RNA3p"
BCs[8]="-a TCAAGTC...$RNA5p -A GACTTGA$RNA3p"
BCs[9]="-a CTGATCC...$RNA5p -A GGATCAG$RNA3p"
BCs[10]="-a AAGCTAC...$RNA5p -A GTAGCTT$RNA3p"
BCs[17]="-a CTCTACC...$RNA5p -A GGTAGAG$RNA3p"
BCs[19]="-a ACTAGC...$RNA5p -A GCTAGT$RNA3p"


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

    dname=`dirname $fastq`; # output directory
    fname=`basename $fastq`; # output file
    fpath=$dname/${fname%_R[12].fastq*}; # file path
    gname=`basename $genome`; # genome name
    ext=${fname#*_R[12].}; # extension of input files
    BCnum=${fname#*BC}; # Barcode number
    BCnum=${BCnum%_*};
    # base filename used for output
    outbase=$dname/$gname/${fname%_R[12].fastq*};

    # don't process unbarcoded reads
    [ $BCnum == "" ] && continue;
    # get appropriate adaptor seqs for this barcode
    ADSEQ=${BCs[$BCnum]};

    # make a directory for the specified genome
    [ -d $dname/$gname ] || mkdir $dname/$gname;


    echo "**************   BEGIN $fname   **************"

    if [ ! -e ${fpath}_R1.trim.fastq ] && [ ! -e ${outbase}.bam ]
    then
        # ADAPTOR TRIMMING
        # -m 18 = require at least 18 bp for alignment
        # -O 2 = trim if at least 2 bases match adaptor sequence
        cutadapt $ADSEQ -m 18 -O 2 --cores=$nproc \
        -o ${fpath}_R1.trim.fastq -p ${fpath}_R2.trim.fastq ${fpath}_R1.$ext ${fpath}_R2.$ext;
        
        # exit with error if file wasn't made
        [ -e ${fpath}_R1.trim.fastq ] || exit 1;

        # bowtie2 alignment options:
        # --very-sensitive alignment parameters
        # -X 1000 = maximum insert size
        # --no-mixed = no unpaired alignments
        # --no-discordant = no alignments > 1000 bp away from eachother
        # --no-unal = no unaligned reads
        bowtie2 -p $nproc -x $genome \
        --very-sensitive -X 1000 \
        --no-mixed --no-discordant --no-unal \
        -1 ${fpath}_R1.trim.fastq -2 ${fpath}_R2.trim.fastq \
        > ${outbase}.sam;
        
        # exit with error if file wasn't made
        [ -e ${outbase}.sam ] || exit 1;
        
        # optional: remove trimmed fastq
        rm ${fpath}_R[12].trim.fastq;
        
        # sort & convert to BAM
        samtools sort -@ $nproc -n -o ${outbase}.bam ${outbase}.sam;
        # exit with error if file wasn't made
        [ -e ${outbase}.bam ] || exit 1;
        
        # recommended: remove SAM
        rm ${outbase}.sam;
    fi
    
    # Spike-in input (BC19) used a single-end sequencing adaptor
    #echo "Aligning BC19 using R1 only..."
    #(bowtie2 -p $nproc -x $genome --no-unal -U ${fpath}_BC19_R1.fastq.gz |
    #awk '/^[[:space:]]*@/{ print; next }{ print | "sort -k3,3 -k4,4n" } END{print}') |
    #samtools view -bS -@ 8 '-' > ${outbase}_BC19.bam;
    #./CountSpikeIns.r -i ${outbase}_BC19.bam
    
    echo "************** FINISHED $fname **************";
done

