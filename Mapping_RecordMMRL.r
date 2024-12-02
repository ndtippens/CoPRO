#!/usr/bin/Rscript

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFiles))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyr))
register(MulticoreParam(7));
source("./CoPRO_Functions.r");

for( len in 19:24 ) {
    fpath = paste0("../data/CoPRO_hg19_dm6/", len, "mers.bam");
    mappable = GRanges( readGAlignments( fpath, use.names=F, param=ScanBamParam(mapqFilter=30)) );
    mappable = compute_IDs(mappable);
    
    # each mappable N-mer identifies 2 RNAs (one on each strand)
    # ID5 = position of 5'RNA on + strand
    # ID3 = position of 5'RNA on - strand
    MapID5 = unique(c( mappable$ID5, -mappable$ID3 ));
    MapID5 = MapID5[order(MapID5)];
    save(MapID5, file=paste0("../data/Mappable", len, "mers.Rdata"));
}

