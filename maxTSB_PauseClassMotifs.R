#!/usr/bin/Rscript
suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(GenomicRanges)
    library(lattice)
    library(foreach)
    library(doParallel)
    library(ggplot2)
    library(scales)
    library(rtfbsdb)
});

setwd("~/Rmain/CoPRO/CoPRO/")

hsdb = CisBP.extdata("Homo_sapiens")
genome = "~/bin/genomes/hg19/hg19.2bit"
hs.gencode = "~/bin/genomes/hg19/gencode.v19.transcripts.gtf"

plogP = function(x) {
    return(round(-log10(x), 2));
}

drawHeatmap = function(motifs, TSBs, tfname) {
    mWidth = motifs$chromEnd[1] - motifs$chromStart[1];
    tWidth = motifs$peakEnd[1] - motifs$peakStart[1];
    tIDs = paste0(motifs$chrom, ':', motifs$peakStart);
    tNames = paste0(TSBs[,1], ':', TSBs[,2]-14);

    isMinus = (motifs$strand == '-');
    uTSBs = unique(tIDs);
    nTSBs = nrow(TSBs);

    # require at least 100 elements with a motif for analysis
    if(length(uTSBs) < 100 | sum(isMinus) == 0 | sum(!isMinus) == 0)
        return();
    tMap = match(tIDs,tNames);
    motifs$Peak=TSBs[tMap, 'Peak'];
    motifs$Class=TSBs[tMap, 'Class'];
    motifs$Tstrand=TSBs$strand[tMap];
    motifs$mStart = motifs$chromStart - motifs$peakStart + 1;


    pmat = matrix(0, nrow=length(uTSBs), ncol=tWidth);
    uIdx = match(uTSBs, tIDs);
    rownames(pmat) = uTSBs;

    # label pause base
    idx = cbind( 1:nrow(pmat), 114+motifs[uIdx,'Peak'] );
    pmat[idx] = 2;

    UpBases = 100*nTSBs;
    PsBases = sum(TSBs[,'Peak']);
    DownBases = tWidth*nTSBs - UpBases - PsBases;

    mpos = ifelse( motifs$Tstrand == '+', motifs$mStart, tWidth-motifs$mStart-mWidth );
    isMinus = motifs$strand == as.character(motifs$Tstrand);
    idx = cbind( match(tIDs, uTSBs), mpos );
    posMtfs = idx[!isMinus,2];
    negMtfs = idx[isMinus,2];
    PosEles = length(unique(tMap[!isMinus]));
    NegEles = length(unique(tMap[isMinus]));

#    fwdTSSi = motifs[!isMinus,'fwdTSS'] - motifs[!isMinus,'peakStart'];
#    revTSSi = motifs[!isMinus,'revTSS'] - motifs[!isMinus,'peakStart'];
#    PosOutMtf = sum( posMtfs < (revTSSi-60) | posMtfs > (fwdTSSi+60) );
#    PosPerMtf = sum( posMtfs >= (revTSSi-60) & posMtfs < revTSSi );
#    PosPerMtf = sum( posMtfs <= (fwdTSSi+60) & posMtfs > fwdTSSi ) + PosPerMtf;
#    PosCtrMtf = sum( posMtfs <= fwdTSSi & posMtfs >= revTSSi );
#    ctable = matrix( c(PosPerMtf, PerBases-PosPerMtf, PosOutMtf, OutBases-PosOutMtf), nrow=2, ncol=2 );
#    SigPosPeriph = fisher.test( ctable, alternative='g' )$p.value;
#    ctable = matrix( c(PosCtrMtf, CtrBases-PosCtrMtf, PosOutMtf, OutBases-PosOutMtf), nrow=2, ncol=2 );
#    SigPosCtr = fisher.test( ctable, alternative='g' )$p.value;
#
#    fwdTSSi = motifs[isMinus,'fwdTSS'] - motifs[isMinus,'peakStart'];
#    revTSSi = motifs[isMinus,'revTSS'] - motifs[isMinus,'peakStart'];
#    NegOutMtf = sum( negMtfs < (revTSSi-60) | negMtfs > (fwdTSSi+60) );
#    NegPerMtf = sum( negMtfs >= (revTSSi-60) & negMtfs < revTSSi );
#    NegPerMtf = sum( negMtfs <= (fwdTSSi+60) & negMtfs > fwdTSSi ) + NegPerMtf;
#    NegCtrMtf = sum( negMtfs <= fwdTSSi & negMtfs >= revTSSi );
#    ctable = matrix( c(NegPerMtf, PerBases-NegPerMtf, NegOutMtf, OutBases-NegOutMtf), nrow=2, ncol=2 );
#    SigNegPeriph = fisher.test( ctable, alternative='g' )$p.value;
#    ctable = matrix( c(NegCtrMtf, CtrBases-NegCtrMtf, NegOutMtf, OutBases-NegOutMtf), nrow=2, ncol=2 );
#    SigNegCtr = fisher.test( ctable, alternative='g' )$p.value;
#    write( paste( tfname, plogP(SigPosPeriph), plogP(SigPosCtr), plogP(SigNegPeriph), plogP(SigNegCtr) ), file='' );

    for( i in 2:8 ) {
        idx[,2] = idx[,2] + 1;
        pmat[idx[!isMinus,]] = 3;
        pmat[idx[isMinus,]] = 4;
    }

    xats  = c(1, 61, 121);
    xlabs = c(-60, 'maxTSB', 60);
    pdf(file=paste0("../plots/motifs/", tfname, ".pdf"), width=4, height=6);
    image(1:121, 1:nrow(pmat), t(pmat[,(21:141)+20+14]), ylab=paste(nrow(pmat), 'max TSBs'),
        xlab='Position (bp)', yaxt='n', xaxt='n', main=tfname, useRaster=T,
        col=c('white', 'black', '#C93A27', '#2D5DC9'),
    );
    axis(1, at=xats, labels=xlabs);
    dev.off();
    
    print(summary(factor(motifs$Class)));
    return();
}


load("../data/newK562_CoPRO_maxTSB.Rdata");
load("../data/Mappable19mers.Rdata");
#pTSB = import.bed("../out/newK562_CoPRO_PauseClasses.bed");
pTSB = mTSB[ mTSB$ID5 %in% MapID5 & mTSB$C >= 2 ];
#mcols(pTSB) = mcols(pTSB)[,2]
#colnames(mcols(pTSB)) = c('PsClass');
# plot all motifs relative to the initiation base (and sorted by distance to the max pause)
pTSB = pTSB[ order(pTSB$PsClass, width(pTSB), decreasing=T) ];
TSBprom = promoters( pTSB, upstream=100, downstream=120 );

hsmotifs = tfbs.createFromCisBP( hsdb );
bw.plus  = "./data/CoPRO_hg19_dm6/CoPRO_AllMerge_Uncapped_3p_pl.bw";
bw.minus = "./data/CoPRO_hg19_dm6/CoPRO_AllMerge_Uncapped_3p_mn.bw";

xmtfs = tfbs.selectExpressedMotifs( hsmotifs, genome, hs.gencode, bw.plus, bw.minus, seq.datatype="GRO-seq", ncores=7 );
tfIDs = apply(xmtfs@tf_info[,c('TF_Name', 'Motif_ID')], 1, paste, collapse=' ')
tfIDs = substr(tfIDs, 0, sapply(tfIDs, nchar, USE.NAMES=F)-5)

TSBdf = as.data.frame(TSBprom)[,-4];
colnames(TSBdf)[1] = 'chr';
TSBdf$Peak = width(pTSB);
TSBdf$Class = pTSB$PsClass;
TSBdf = TSBdf[,c(1:3, 5:6, 4)];
TSBmot = tfbs.scanTFsite( xmtfs, genome, TSBdf, ncores=7 );
for( i in 1:length(tfIDs) ) {
    drawHeatmap( TSBmot$result[[i]], TSBdf, tfIDs[i] );
}
