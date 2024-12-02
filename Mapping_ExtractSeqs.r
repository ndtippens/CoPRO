#!/usr/bin/Rscript

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFiles))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
registerDoParallel(cores=7);
source("./CoPRO_Functions.r");


# Normally, people rely on mapping quality scores (MapQ)
# to decide whether to keep a read. However, our data
# captures nascent elongating RNAs. Therefore many
# unmappable reads can be saved by looking for longer
# RNAs sharing a 5' end.

load('../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Pooled.Rdata');

# Get 5' sites from mapped reads
good = as.data.table( mcols(Pooled[Pooled$C+Pooled$R>0]) );
good = good[, .(totC = sum(C), totR = sum(R)), by=ID5 ];
good = good[totC >= 5,];

CheckMap = Pooled[ Pooled$ID5 %in% good$ID5 ];
CheckMap = CheckMap[ !duplicated(CheckMap$ID5) ];
CheckMap$ID3 = NULL;
CheckMap$C = NULL;
CheckMap$U = NULL;
CheckMap$R = NULL;
CheckMap$E = NULL;
rm(Pooled);
invisible(gc());


lens = 19:24;
genomef = "~/bin/genomes/hg19/hg19.2bit";
CheckMap = resize(CheckMap, width=max(lens));
seqs = import.2bit(genomef, which=CheckMap);
mstr = strand(CheckMap) == '-';
seqs[mstr] = reverseComplement(seqs[mstr]);
seqs = unique(seqs);
for( len in lens ) {
	writeXStringSet(subseq(seqs, start=1, width=len), paste0("../data/", len, "mers.fa.gz"), compress=T );
}

