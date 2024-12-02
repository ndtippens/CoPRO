#!/usr/bin/Rscript

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFiles))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
registerDoParallel(cores=7);
source("./CoPRO_Functions.r");

load("../data/newK562_CoPROcalls.Rdata");
load("../data/newK562_CoPRO_maxTSB.Rdata");

eTRE = TREs[ TREs$nTSS & TREs$pTSS ];
unst = eTRE$nS == 0;
rm(TREs);
invisible(gc());

TSS = TSS[ order(-TSS$Cap) ];
hits = findOverlaps( TSS, eTRE );
# find first & second max TSS
maxt = hits@from[!duplicated(hits@to)];
TSS1 = TSS[  maxt ];
TSS  = TSS[ -maxt ];
hits = findOverlaps( TSS, eTRE );
maxt = hits@from[!duplicated(hits@to)];
TSS2 = TSS[  maxt ];

hits = findOverlaps( mTSB, TSS1 );
TSB1 = mTSB[ hits@from ];
hits = findOverlaps( mTSB, TSS2 );
TSB2 = mTSB[ hits@from ];

hits = findOverlaps( TSS1, eTRE );
TSS1 = TSS1[ order(hits@to) ];
hits = findOverlaps( TSS2, eTRE );
TSS2 = TSS2[ order(hits@to) ]; # now 1:1 correspondence with eTRE ordering
hits = findOverlaps( TSB1, eTRE );
TSB1 = TSB1[ order(hits@to) ];
hits = findOverlaps( TSB2, eTRE );
TSB2 = TSB2[ order(hits@to) ];

pclass = mTSB$PsClass;
summary(pclass);
PsVstab = chisq.test(pclass, mTSB$Stability);
print(PsVstab);
print(PsVstab$observed / PsVstab$expected);
PsVloc  = chisq.test(pclass, mTSB$ID5 %in% c(TSB1$ID5, TSB2$ID5));
print(PsVloc);
print(PsVloc$observed / PsVloc$expected);
PsVloc  = chisq.test(pclass, mTSB$ID5 %in% TX$maxTSS);
print(PsVloc);
print(PsVloc$observed / PsVloc$expected);


pTSB = c(TSB1, TSB2);
summary(pTSB$PsClass);
temp = pTSB[pTSB$PsClass != 'Both' & pTSB$Stability != 'None'];
PsVstab = chisq.test(temp$PsClass, temp$Stability);
print(PsVstab);
print(PsVstab$observed / PsVstab$expected);
summary(pTSB$Location);
temp = pTSB[pTSB$PsClass != 'Both' & pTSB$Location == 'Promoter'];
PsVloc  = chisq.test(temp$PsClass, temp$Location);
print(PsVloc);
print(PsVloc$observed / PsVloc$expected);

pTSB = TSB1[ TSB1$Location == "Promoter" ];

binSize=1;
nBins=201;
xat = 101 + c(-100, 0, 100);
xlabels = c(-100, 0, 100);
pcols = colorRampPalette(c("white", "firebrick"))(501);
mkcols = c('gray', 'dodgerblue', 'forestgreen');

negstr = strand(pTSB) == "-";
# center on pause base
pauseWindow = GenomicRanges::shift( pTSB, ifelse(negstr, -pTSB$Peak, pTSB$Peak ));
pauseWindow = promoters( pauseWindow, upstream=floor(binSize*nBins/2), downstream=ceiling(binSize*nBins/2) );
negstr = which(negstr);
earlyc = which(pTSB$Peak <= 32);
latec  = which(pTSB$Peak >  32);

Factors = c('unphospho_rep1', 'unphospho_rep2', 'Ser2_rep1', 'Ser2_rep2', 'Ser5_rep1', 'Ser5_rep2', 'total_rep1');
rawm = matrix(0, nrow=length(pauseWindow), ncol=nBins);
EarlyM = matrix(0, nrow=nBins, ncol=length(Factors));
colnames(EarlyM) = Factors;
LateM  = EarlyM;

pdf("../plots/maxTSB_mNETseq.pdf", width=5, height=5);
    for( fn in Factors ) {
        msgout(fn);
        bwf = BigWigFile( paste0('/media/nate/Diatom/Sequencing/ChIP/HeLa-mNETseq/HeLa_ANET_', fn, '_F.bw') );
        rawm = summary( bwf, pauseWindow, size=nBins, type='mean', defaultValue=0, as='matrix' );
        bwf = BigWigFile( paste0('/media/nate/Diatom/Sequencing/ChIP/HeLa-mNETseq/HeLa_ANET_', fn, '_R.bw') );
        rawm[ negstr,] = -summary( bwf, pauseWindow, size=nBins, type='mean', defaultValue=0, as='matrix' )[negstr,];
        rawm[ negstr, ] = rawm[ negstr, nBins:1 ];

        pAll = log2(1+rawm[order(-pTSB$Peak),]);
        EarlyM[,fn] = colSums(1+rawm[earlyc,]);
        LateM[,fn] = colSums(1+rawm[latec,]);

        print( levelplot( x=t(pAll),
            main=fn, xlab="Distance from maxPause (nt)", ylab="",
            col.regions=pcols,
            useRaster=T, aspect="fill", at=seq(0, max(pAll), length.out=500),
            scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
        ) );
    }

    unprat = (EarlyM[,'unphospho_rep1']/EarlyM[,'total_rep1'])/(LateM[,'unphospho_rep1']/LateM[,'total_rep1']);
    Ser2rat = (EarlyM[,'Ser2_rep2']/EarlyM[,'total_rep1'])/(LateM[,'Ser2_rep2']/LateM[,'total_rep1']);
    Ser5rat = (EarlyM[,'Ser5_rep2']/EarlyM[,'total_rep1'])/(LateM[,'Ser5_rep2']/LateM[,'total_rep1']);
    Signal = cbind(unprat, Ser2rat, Ser5rat);
    matplot( x=seq(-100, 100, binSize), y=Signal, type='l', lty=1, lwd=2, ylim=c(0.25, 4),
        log='y', xlab='Distance from maxPause', ylab='Early / Late', col=mkcols, main="HeLa mNET-seq");
    legend( 'topleft', legend=c('unphospho', 'Ser2', 'Ser5'), lty=1, col=mkcols, lwd=2, bty='n');
dev.off();
