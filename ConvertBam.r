#!/usr/bin/Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(rtracklayer)
    library(GenomicRanges)
});
register(MulticoreParam(6));
source("./CoPRO_Functions.r");

froot  = '../data/CoPRO_hg19_dm6/CoPRO_AllMerge_';
Barcodes = list(
    BC2="Uncapped_Rep2",
    BC3="Capped_Rep2",
    BC4="RppH_Rep2",
    BC6="Uncapped_Rep1",
    BC8="Capped_Rep1",
    BC9="RppH_Rep1",
    BC19="SpikeIns"
);

collapse_reads = function(reads) {
    uniq = unique(reads);
    # sum over duplicates to get a count for each unique 5'/3' end
    uniq$count = countOverlaps( uniq, reads, type="equal" );
    return( uniq );
}

yield.bam = function(X) {
    y = GRanges(readGAlignmentPairs(X, use.names=T, param=ScanBamParam(mapqFilter=30, flag=scanBamFlag(isSecondaryAlignment = F))));
    names(y) = NULL;
    return(y);
}

map.bam = function(al.gr) {
    isHs = seqnames(al.gr) %in% hg19chr;
    isDm = substr(seqnames(al.gr), 0, 3) == "dm6";
    isChr = substr(seqnames(al.gr), 0, 3) == "chr";
    human = al.gr[ isHs ];

    SI = al.gr[ !isDm & !isChr ];
    SIlen = count(width(SI));
    hit = paste(seqnames(SI), end(SI));
    correct = SI[(hit %in% sizes) & start(SI) < 3];

    spikeins = c(
        dm6 = sum(isDm),
        ppp1 = sum(seqnames(correct) == "ppp1"),
        ppp2 = sum(seqnames(correct) == "ppp2"),
        ppp3 = sum(seqnames(correct) == "ppp3"),
        ppp4 = sum(seqnames(correct) == "ppp4"),
        ppp5 = sum(seqnames(correct) == "ppp5"),
        Cap1 = sum(seqnames(correct) == "Cap1"),
        Cap2 = sum(seqnames(correct) == "Cap2"),
        Cap3 = sum(seqnames(correct) == "Cap3"),
        Cap4 = sum(seqnames(correct) == "Cap4"),
        Cap5 = sum(seqnames(correct) == "Cap5")
    );

    return(list(human, spikeins));
}

reduce.bam = function(x,y) {
    # merge GRanges
    x[[1]] = append(x[[1]], y[[1]]);
    x[[2]] = x[[2]] + y[[2]];
    # clean up
    rm(list=c('y'));
    gc();

    # print the number of readpairs processed
    #msgout(pn(length(x[[1]])), 'mapped human reads');
    return(x);
}



for( bc in names(Barcodes) ) {
    ctFile = paste0(froot, bc, '.bam');

    # load sizes of spike-ins (only count reads mapping near first & last base)
    header = scanBamHeader(ctFile, what='targets')[[1]]$targets;
    cnames= names( header );
    csizes= as.numeric(header);
    sizes = c( paste(cnames, csizes), paste(cnames, csizes-1), paste(cnames, csizes-2) );

    msgout( "Processing ", ctFile );
    infile = BamFile(ctFile, yieldSize=0.3 * 10^6);
    aligned = reduceByYield( infile, yield.bam, map.bam, reduce.bam, parallel=T );
    
    # remove extra seqinfo from human reads & apply standard chr ordering (1-22, X, Y, M)
    sinfo = seqinfo( aligned[[1]] );
    newchr = match( hg19chr, names(sinfo) );
    seqinfo(aligned[[1]], new2old=newchr ) = sinfo[names(sinfo)[newchr]];
    msgout(pn(length(aligned[[1]])), 'mapped human reads');
    print(aligned[[2]]);

    # compute coverage from identical reads => 'count' column
    seqlib = collapse_reads(aligned[[1]]);

    # Swap strand to match synthesized RNA
    strand(seqlib) = ifelse( strand(seqlib) == "+", '-', '+' );
    RNApl = strand(seqlib) == "+";
    RNA5 = ifelse( RNApl, start(seqlib), end(seqlib) );
    # remove biotin run-on base
    RNA3 = ifelse( RNApl, end(seqlib)-1, start(seqlib)+1 );
     start(seqlib) = ifelse( RNApl, RNA5, RNA3 );
       end(seqlib) = ifelse( RNApl, RNA3, RNA5 );
    rm(RNA5, RNA3);
    invisible(gc());

    # compute universal 5p and 3p IDs for each alignment
    seqlib = compute_IDs(seqlib);
    # sort by 5p and 3p coordinates
    seqlib = seqlib[ order(seqlib$ID5), ];

    aligned[[1]] = seqlib;
    # save in Rdata format
    save(aligned, file=paste0('../data/CoPRO_hg19_dm6/CoPRO_AllMerge_', Barcodes[[bc]], '.Rdata'));
    rm(aligned, seqlib);
    invisible(gc());
}

