#!/usr/bin/Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicFiles)
    library(GenomicRanges)
    library(data.table)
    library(rtracklayer)
    library(RColorBrewer)
    library(foreach)
    library(doParallel)
    library(ggplot2)
});
registerDoParallel(cores=7);
source("./CoPRO_Functions.r");
genomef = "~/bin/genomes/hg19/hg19.2bit";



all_pause = function( h ) {
    h = unname(h);
    # order positions from max height -> min
    p = order(h, decreasing = T);
    # final position required to capture 50% of signal
    endp = which(cumsum(h[p]) > 0.5)[1];
    p = p[ 1 : endp ];

    # cluster positions within 5 bases of eachother

    return( p );
}

pause_stats = function( x ) {
       x = unname(x);
    cumx = cumsum(x);
    for( w in 1:15 ) {
        pkSum = diff( cumx[1:61], lag=w );
        if( max(pkSum) >= 0.5 )
            return( c(1+which.max(pkSum), w) );
    }
    return( c(NA,NA) );
}


load("../data/newK562_CoPROcalls.Rdata");



length(unique(sTSB$TSS));
maxTSB = as.data.table( mcols(sTSB) );
# determine which site has max capped reads at each TSS
maxTSB[, isMax := (which.max(C) == seq_len(.N)), by=TSS ];
mTSB = sTSB[ maxTSB[,isMax] ];
mTSB = mTSB[ order(mTSB$ID5) ];
mTSB$Stability = TSS$Stability[ match(mTSB$ID5, TSS$ID5) ];
length(mTSB);

TSSid = mTSB$ID5;
nTSS = length(mTSB);

# gather reads for plotting
windowSize = 400;
pReads = TSB[ width(TSB) > 18 & width(TSB) <= windowSize & TSB$ID5 %in% TSSid ];
ID5long = unique(TSB$ID5[width(TSB)>=100]);
rm(TSB);
invisible(gc());

cmat = matrix(data=0, nrow=nTSS, ncol=windowSize);
rownames(cmat) = TSSid;
colnames(cmat) = as.character(1:windowSize);
umat=cmat;
rmat=cmat;
tmat=cmat;
cmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$C;
umat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$U;
rmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$R;
tmat = cmat + umat + rmat;

csum = rowSums(cmat);
usum = rowSums(umat);
rsum = rowSums(rmat);
tsum = rowSums(tmat);
maxC = apply(cmat, 1, max);
maxU = apply(umat, 1, max);
maxR = apply(rmat, 1, max);
maxT = apply(tmat, 1, max);
posT = apply(tmat, 1, which.max);
count_TSS = function(x) { return( sum(x>0) ); };
capN = apply(cmat, 1, count_TSS);
uncN = apply(umat, 1, count_TSS);
rppN = apply(rmat, 1, count_TSS);
totN = apply(tmat, 1, count_TSS);

mTSB$PsOcc  = rowSums(cmat[,20:55]);
mTSB$EscOcc = rowSums(cmat[,100:400]);
mTSB$LateOcc = rowSums(cmat[,33:55]);

cmat = cmat/csum;
rmat = rmat/rsum;
umat = umat/usum;
tmat = tmat/tsum;


SCS = function(x) { sum(cumsum(x)) };
classifyPause = function( rdmat ) {
    totPs = rowSums(rdmat[,20:60]);
    totPs[!totPs] = 1;
    pscore = rowSums(rdmat[,20:32])/totPs;
    # flate = rowSums(rdmat[,36:60])/totPs;
    cutoffs = quantile(pscore);
    cls = ifelse( pscore <= cutoffs[2], 2, 3 );
    cls[ pscore >= cutoffs[4] ] = 1;
    return(cls);
}


PsClasses=c('Early', 'Late', 'Both');
mTSB$PsClass = factor(PsClasses[classifyPause(cmat)]);
mTSB$Peak = 0;
mTSB[mTSB$PsClass=='Early']$Peak = 19+apply(cmat[as.character(mTSB[mTSB$PsClass=='Early']$ID5),][,20:32], 1, which.max);
mTSB[mTSB$PsClass=='Late']$Peak = 31+apply(cmat[as.character(mTSB[mTSB$PsClass=='Late']$ID5),][,32:55], 1, which.max);
mTSB[mTSB$PsClass=='Both']$Peak = 19+apply(cmat[as.character(mTSB[mTSB$PsClass=='Both']$ID5),][,20:55], 1, which.max);

save(mTSB, file="../data/newK562_CoPRO_maxTSB.Rdata");

eTRE = TREs[ TREs$nTSS & TREs$pTSS & TREs$Cap >= 3 ];
TSS = TSS[ order(-TSS$Cap) ];
hits = findOverlaps( TSS, eTRE );
# find first & second max TSS
maxt = hits@from[!duplicated(hits@to)];
TSS1 = TSS[  maxt ];
TSS  = TSS[ -maxt ];
hits = findOverlaps( TSS, eTRE );
maxt = hits@from[!duplicated(hits@to)];
TSS2 = TSS[  maxt ];

mTSB = mTSB[ mTSB$ID5 %in% c(TSS1$mTSB, TSS2$mTSB) ];
mTSB = mTSB[ seqnames(mTSB) != 'chrX' ];

# require 21mer mappability (20 nt + 1 run-on), + Capped & RppH expression
load("../data/Mappable21mers.Rdata");
mTSB = mTSB[ mTSB$ID5 %in% MapID5 ];

pdf("../plots/maxTSB_PauseClasses.pdf", width=3, height=4);

with( mTSB[mTSB$PsOcc>1],
plot( y=10*PsOcc, x=(EscOcc+.1)/PsOcc, ylab="Pause Reads",
    xlab="Escape Index", pch=19, cex=0.25, col=tblack,
    ylim=c(3,1000), xlim=c(0.1,4) )
);

boxplot(log2(1+mTSB$C) ~ mTSB$PsClass, xlab="Pause Class", ylab="Total Reads (log2)",
    outline=F, notch=T
);
t.test(log2(1+mTSB$C) ~ mTSB$PsClass, subset=mTSB$PsClass != 'Both');

boxplot(log2(1+mTSB$EscOcc) ~ mTSB$PsClass, xlab="Pause Class", ylab="Reads >60 nt (log2)",
    outline=F, notch=T
);
t.test(log2(1+mTSB$EscOcc) ~ mTSB$PsClass, subset=mTSB$PsClass != 'Both');

dev.off();


PsBs = GenomicRanges::shift(mTSB, shift=ifelse(strand(mTSB) == '+', mTSB$Peak, -mTSB$Peak));
#PsBs = GenomicRanges::resize(TSB[width(TSB)<=55], width=1, fix="end");
PsBed = mTSB;
start(PsBed) = ifelse(strand(PsBed) == '+', start(mTSB),   end(PsBs));
  end(PsBed) = ifelse(strand(PsBed) == '+',   end(PsBs), start(mTSB));
mcols(PsBed) = as.integer(mcols(PsBed)[,c('PsClass')]);
colnames(mcols(PsBed))='score';
export(PsBed, "../out/K562_CoPRO_PauseClasses.bed");


rdwin = width(TSB) >= 20 & width(TSB) <= 55;
hiTSB = sTSB$ID5[ sTSB$C >= 3 & sTSB$ID5 %in% TSB$ID5[rdwin] ];
hiTSB = TSB[ TSB$ID5 %in% hiTSB & rdwin ];
hiTSB = hiTSB[ order(hiTSB$ID5, hiTSB$C, decreasing=T) ];
PsBs = GRanges();
for(i in 1:3) {
    PsBs = append(PsBs, hiTSB[ !duplicated(hiTSB$ID5) ]);
    hiTSB = hiTSB[ duplicated(hiTSB$ID5) ];
}

PsBs = GenomicRanges::resize(PsBs, width=1, fix="end");

# get sequence
flank3 = GenomicRanges::promoters(PsBs, upstream=5, downstream=6);
flank3 = flank3[ sample(1:length(flank3), 10^5, prob=flank3$C/sum(flank3$C), replace=T) ];
## NOTE: import.2bit ignores strand, so we must manually complement the sequence.

TSS3 = TSS[TSS$ID5 %in% flank3$TSS];
flank3 = flank3[flank3$TSS %in% TSS3$ID5];
TSS3 = TSS3[match(flank3$TSS, TSS3$ID5)];
out = cbind( Location=as.character(TSS3$Location), RNAclass=as.character(TSS3$Stability), Name=TSS3$TUname, Occupancy=flank3$C, PsClass=as.character(flank3$PsClass), CapR=TSS3$CapR );
write.csv(out, "../out/AllPauseSeq_info.csv");

seq3 = import.2bit(genomef, which=flank3);
seq3[ strand(flank3) == '-' ] = reverseComplement(seq3[ strand(flank3) == '-' ]);
export(seq3, "../out/K562_AllPause_seqs.fasta.gz")
export(seq3[flank3$Peak >= 20 & flank3$Peak < 36], "../out/K562_EarlyPause_seqs.fasta.gz");
export(seq3[flank3$Peak > 36], "../out/K562_LatePause_seqs.fasta.gz");
export(seq3[flank3$PsClass == 'Both'], "../out/K562_BothPause_seqs.fasta.gz");



rdwin = width(TSB) >= 80 & width(TSB) <= 115;
hiTSB = sTSB$ID5[ sTSB$C >= 3 & sTSB$ID5 %in% TSB$ID5[rdwin] ];
hiTSB = TSB[ TSB$ID5 %in% hiTSB & rdwin ];
hiTSB = hiTSB[ order(hiTSB$ID5, hiTSB$C, decreasing=T) ];
PsBs = GRanges();
for(i in 1:3) {
    PsBs = append(PsBs, hiTSB[ !duplicated(hiTSB$ID5) ]);
    hiTSB = hiTSB[ duplicated(hiTSB$ID5) ];
}
PsBs = GenomicRanges::resize(PsBs, width=1, fix="end");
flank3 = GenomicRanges::promoters(PsBs, upstream=5, downstream=6);
flank3 = flank3[ sample(1:length(flank3), 10^5, prob=flank3$C/sum(flank3$C), replace=T) ];
seq3 = import.2bit(genomef, which=flank3);
seq3[ strand(flank3) == '-' ] = reverseComplement(seq3[ strand(flank3) == '-' ]);
export(seq3, "../out/K562_All3pSeqs.fasta.gz");

TSSid = mTSB$ID5;
nTSS = length(mTSB);
print(nTSS);


# gather reads for plotting
windowSize = 400;
pReads = pReads[ pReads$ID5 %in% TSSid ];

cmat = matrix(data=0, nrow=nTSS, ncol=windowSize);
rownames(cmat) = TSSid;
colnames(cmat) = as.character(1:windowSize);
umat=cmat;
rmat=cmat;
tmat=cmat;
cmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$C;
umat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$U;
rmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$R;
tmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$C + pReads$U + pReads$R + pReads$E;

csum = rowSums(cmat);
usum = rowSums(umat);
rsum = rowSums(rmat);
tsum = rowSums(tmat);
maxC = apply(cmat, 1, max);
maxU = apply(umat, 1, max);
maxR = apply(rmat, 1, max);
maxT = apply(tmat, 1, max);
posT = apply(tmat, 1, which.max);
count_TSS = function(x) { return( sum(x>0) ); };
capN = apply(cmat, 1, count_TSS);
uncN = apply(umat, 1, count_TSS);
rppN = apply(rmat, 1, count_TSS);
totN = apply(tmat, 1, count_TSS);
cmat = cmat/csum;
rmat = rmat/rsum;
umat = umat/usum;
tmat = tmat/tsum;

pcc = mTSB$PsClass;
pscols = c('dodgerblue3', 'forestgreen', 'firebrick', 'darkorchid4', 'navajowhite4', 'gray');
levels(pcc) = adjustcolor( pscols, alpha.f = 0.2 );

#plot(pfrac[,1], pfrac[,2], col=as.character(pcc), xlab='Early Pause Fraction', ylab="Late Pause Fraction", pch=19, cex=0.5);
#legend("topright", legend=PsClasses, title='Pause Classes', pch=19, col=levels(pcc))
#scatterplot3d(x=pfrac, color=as.character(pcc), xlab='Early Pause Fraction', ylab='Late Pause Fraction', zlab="Late Pause Fraction", pch=19, cex.symbols=0.5, angle=180);


# check for ChIP factor binding around each TSB
TSBprom = GenomicRanges::promoters(mTSB, upstream=120, downstream=120);
fnames = list.files('~/Rmain/TREs/sources/k562/ENCODE/pk/');
ChIPmat = foreach( fname=fnames, .combine=`cbind` ) %dopar% {
    pkfile = paste0('~/Rmain/TREs/sources/k562/ENCODE/pk/', fname);
    bedf = as.data.frame( read.delim( gzfile(pkfile), header=F ) )[,1:3];
    colnames(bedf) = c('seqname', 'start', 'end');
    fname=sub('.bed.gz', '', fname);
    bedg = as( bedf, 'GRanges' );
    hits = countOverlaps( TSBprom, bedg, type='any', ignore.strand=T );
    hits=as.data.frame(hits);
    colnames(hits) = fname;
    return(hits);
}

fit = pamk(t(ChIPmat), krange=2, usepam=F, critout=F);
clust = fit$pamobject$clustering;
clord = order(unname(clust));
tord  = order(mTSB$PsClass);

acols = colorRampPalette(c("steelblue", "white", "firebrick"))(1000);
bcols = colorRampPalette(c("white", "#333333"))(1000);

plotA = rbind(  Cls = as.numeric(mTSB$PsClass)[tord] / 6,
                Loc = as.numeric(mTSB$Location)[tord] / 3
        );
amap = levelplot( x=plotA,
    main="max TSBs", ylab="Site", xlab="",
    cuts=1000, colorkey=F, col.regions=acols,
    aspect="fill", at=seq(0, 1, length.out=1000), useRaster=T,
    scales=list(
        x=list(draw=F),
        y=list(draw=F)
    )
);
bmap = levelplot( x=t(ChIPmat[tord,clord]),
    main="Transcription Factors", xlab="", ylab="",
    cuts=1000, colorkey=list(col=bcols), col.regions=bcols,
    useRaster=T, aspect="fill", at=seq(0, 1, length.out=1000),
    scales=list(x=list(draw=F), y=list(draw=F))
);
grid.arrange(amap, bmap, widths=c(0.2, 0.8), padding=unit(0.2, 'line'));




ftest = foreach( bound = iter(ChIPmat, by='col'), .combine=`rbind` ) %dopar% {
    # skip empty ChIP profiles
    if(!any(bound))
        return(c(1,1));

    pvals = rep(0, 3);
    for( c in PsClasses ) {
        if( sum(bound & mTSB$PsClass == c) > 3 ) {
            #msgout(sum(bound & myTRE$class == c), sum(!bound & myTRE$class == c), sum(bound & myTRE$class != c), sum(!bound & myTRE$class != c) )
            pvals[which(PsClasses==c)] = fisher.test( x=as.factor(mTSB$PsClass == c), y=as.factor(bound>0), alt='greater' )$p.value;
        } else
            pvals[which(PsClasses==c)] = 1;
    }
    return(pvals);
}
rownames(ftest) = colnames(ChIPmat);
colnames(ftest) = PsClasses;

write.csv(ftest, file="../out/maxTSB_PauseClass_FactorEnrich.csv");

sigf = rowMins(ftest) <= 10^-12;
fmat = ftest[sigf,];
fmat[fmat>.01] = 1;
ford = order(-fmat[,2], -fmat[,1]);

fmat = -log10(fmat[ford,]);
fmat[fmat>40] = 40;

ccols = colorRampPalette(c("white", "gold", "firebrick"))(500);
levelplot( x=t(fmat),
    main="Factor Enrichment", xlab="", ylab="",
    cuts=500, colorkey=list(col=ccols), col.regions=ccols,
    useRaster=T, aspect="fill", at=seq(5, 40, length.out=500),
    scales=list(x=list(tck=0), y=list(tck=0, cex=0.5))
);
levelplot( x=t(fmat[,1:2]),
    main="Factor Enrichment", xlab="", ylab="",
    cuts=500, colorkey=list(col=ccols), col.regions=ccols,
    useRaster=T, aspect="fill", at=seq(5, 40, length.out=500),
    scales=list(x=list(tck=0), y=list(tck=0))
);
dev.off();



PsList  = apply(tmat, 1, all_pause);
PsLen = sapply(PsList, length);
bimod = PsList[PsLen==2];
bidiff = sapply(bimod, diff);
#hist(abs(bidiff), breaks=max(bidiff), xlim=c(0,20), xlab="Distance between Pause Peaks");



PsStats = t(apply(tmat, 1, pause_stats));
PsNum = apply(rmat[,17:71], 1, SCS);
PsNum = (PsNum - apply(cmat[,17:71], 1, SCS)) / PsNum;
mTSB$Pdist  = PsStats[as.character(mTSB$ID5),1];
mTSB$Pwidth = PsStats[as.character(mTSB$ID5),2];
pcidx = cbind( as.character(mTSB$ID5), mTSB$Peak );
mTSB$Pcapr  = cmat[ pcidx ] / (cmat+umat)[pcidx];

pTSB = as.data.table( mcols(pReads) );
pTSB[, Pdist  := PsStats[as.character(ID5),1]];
pTSB[, Pwidth := PsStats[as.character(ID5),2]];
pTSB[, w := abs(ID5-ID3)+1 ];
RCapR = 100*round(mTSB$CapR,1);
mTSB$TermR = TSS$TermR[match(mTSB$ID5, TSS$ID5)];


pReads$PPeak = PsStats[as.character(pReads$ID5),1];
pReads$PWid  = PsStats[as.character(pReads$ID5),2];
#AllDist = na.omit(as.vector(PsStats));
AllDist = unlist(PsList);
AllDiff = abs(unlist( sapply(PsList[PsLen>1], diff) ));

pdf("../plots/Pause_Stats.pdf", width=4, height=4);
densityplot(mTSB$Pdist, plot.points=F, col="steelblue3", bw=1, xlab="Peak Distance (bases)", ylab="Fraction of maxTSBs", xlim=c(15,60), lwd=2);
densityplot(sapply(PsList, `[`, 1), plot.points=F, col="steelblue3", bw=1, xlab="Strongest Pause Base", ylab="Fraction of maxTSBs", xlim=c(15,60), lwd=2);
densityplot(sapply(PsList, `[`, 2), plot.points=F, col="steelblue3", bw=1, xlab="2nd Strongest Pause Base", ylab="Fraction of maxTSBs", xlim=c(15,60), lwd=2);
densityplot(sapply(PsList, `[`, 3), plot.points=F, col="steelblue3", bw=1, xlab="3rd Strongest Pause Base", ylab="Fraction of maxTSBs", xlim=c(15,60), lwd=2);
hist(AllDiff, breaks=32, col="steelblue3", xlab="Distance Between Pause Peaks", xlim=c(0,20), main="");

Pcenter = mTSB$Pdist + floor(mTSB$Pwidth/2);
boxplot(100*ElongR ~ Pcenter, data=mTSB, xlab="Pause Position (bases)", ylab="Released Polymerase (%)",
        ylim=c(0,30), col="steelblue3", subset=!is.na(Pwidth) & !is.na(ElongR), outline=F);

boxplot(100*ElongR ~ RCapR, data=mTSB, xlab="Percent Capped at Pause", ylab="Released Polymerase (%)",
        ylim=c(0,50), col="steelblue3", subset=!is.na(CapR) & !is.na(ElongR), outline=F);

boxplot(100*TermR ~ Peak, data=mTSB, xlab="Pause Position (bases)", ylab="Terminating Polymerase (%)",
        subset=!is.na(TermR), notch=T, col="steelblue3", outline=F);

boxplot(100*TermR ~ RCapR, data=mTSB, xlab="Percent Capped at Pause", ylab="Terminating Polymerase (%)",
        subset=!is.na(TermR), notch=T, col="steelblue3", xlim=c(0.8,11.2), ylim=c(0,100), outline=F);

boxplot(100*Pcapr ~ Peak, data=mTSB, xlab="Pause Position", ylab="Percent Capped Reads",
        subset=C>5 & U>=2, notch=T, col="steelblue3", ylim=c(0,100), outline=F);


toPlot = mTSB$U >= 3;
ggplot( data=data.frame(x=Pcenter, y=100*mTSB$CapR)[toPlot,], aes(x=x, y=y)) + stat_bin2d(binwidth=c(1,5)) +
    scale_x_continuous(name = "Pause Position (bases)", limits=c(15,70)) +
     scale_y_continuous(name = "Capping Ratio (%)", limits=c(0,105)) +
    scale_fill_gradientn(name = "Number TSBs", colors=heatcols ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position=c(0.8, 0.25),
        axis.line = element_line(colour = "black")
);

Pcenter = pTSB$Pdist + floor(pTSB$Pwidth/2);
ggplot( data=data.frame(x=Pcenter, y=pTSB$Pwidth), aes(x=x, y=y)) + stat_bin2d(binwidth=1) +
    scale_x_continuous(name = "Pause Position (bases)", limits=c(15,75)) +
     scale_y_continuous(name = "Pause Width", limits=c(-1,17)) +
    scale_fill_gradientn(name = "Number TSBs", colors=heatcols ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position=c(0.8, 0.8),
        axis.line = element_line(colour = "black")
);

dev.off();


# read dREG HD enh/pro calls
#enhpro = import.bed("../../K562_Heatshock/data/K562_FC_NHS_BRs_dREG_HD.bed");
#hits = findOverlaps(pReads, enhpro, ignore.strand=T);
#pReads$dREG = factor('None', levels=c('None', 'Enh', 'Pro'));
#pReads$dREG[hits@from] = ifelse(as.integer(enhpro$name[hits@to]), 'Pro', 'Enh');
#pClass = ifelse(is.na(mTSB$Location), 0, mTSB$Location == 'Promoter');
pClass = mTSB$PsClass;





#rows=round((1:500)*(nrow(cmat)/500));
acols = colorRampPalette(c("white", "firebrick"))(1000);
bcols = colorRampPalette(c("white", "deepskyblue4"))(1000);
ccols = colorRampPalette(c('deepskyblue4', 'white', 'firebrick'))(500);

pdf(file="../plots/maxTSB_Pausing.pdf", width=6, height=6);

plotByClass = function(signal, class) {
    a = colMeans(signal[class==PsClasses[1],], na.rm=T);
    b = colMeans(signal[class==PsClasses[2],], na.rm=T);
    return(t(rbind(a,b)));
}
matplot( x=1:100, y=colMeans(tmat[,1:100]),
    xlim=c(20,65), type='l', col='gray', lty=1, lwd=1,
    main='Total (C+U+R+E)', xlab='RNA Length', ylab="Mean Signal"
);
matplot( x=1:100, y=plotByClass(tmat[,1:100], pClass),
    xlim=c(20,65), type='l', col=pscols, lty=1, lwd=1,
    main='Total (C+U+R+E)', xlab='RNA Length', ylab="Mean Signal"
);
matplot( x=1:100, y=plotByClass(rmat[,1:100], pClass),
    xlim=c(20,65), type='l', col=pscols, lty=1, lwd=1,
    main='RppH', xlab='RNA Length', ylab="Mean Signal"
);
matplot( x=1:100, y=plotByClass(cmat[,1:100], pClass),
    xlim=c(20,65), type='l', col=pscols, lty=1, lwd=1,
    main='Capped', xlab='RNA Length', ylab="Mean Signal"
);
matplot( x=1:100, y=plotByClass(umat[,1:100], pClass),
    xlim=c(20,65), type='l', col=pscols, lty=1, lwd=1,
    main='Uncapped', xlab='RNA Length', ylab="Mean Signal"
);

rows = (1:1000)*floor(nTSS/1000);
pOrder = order(pClass, mTSB$Peak, decreasing=T);
#pOrder = order( colSums( apply(tmat, 1, cumsum) ) );
plotR = t(rmat[,1:100][pOrder,][rows,]);
plotC = t(cmat[,1:100][pOrder,][rows,]);
plotU = t(umat[,1:100][pOrder,][rows,]);
plotT = t(tmat[,1:100][pOrder,][rows,]);
plotR[plotR>0.5] = 0.5;
plotC[plotC>0.5] = 0.5;
plotU[plotU>0.5] = 0.5;
plotT[plotT>0.5] = 0.5;

xat = c(20, 35, 60, 100);
xlabels = xat;



edata = cbind(  Tot = log(tsum[pOrder[rows]]),
                Esc = mTSB$ElongR[pOrder[rows]],
                Loc = mTSB$Location[pOrder[rows]],
                Cls =   pClass[pOrder[rows]]
        );
plotA = t(edata) / apply(edata, 2, max);
amap = levelplot( x=plotA,
    main="Reads", ylab="TSS", xlab="",
    cuts=1000, colorkey=F, col.regions=acols,
    aspect="fill", at=seq(0, 1, length.out=500), useRaster=T,
    scales=list(
        x=list(labels=colnames(edata), rot=90, tck=c(1,0), cex=0.7),
        y=list(draw=F)
    )
);
cmap = levelplot( x=plotT,
    main="Total RNA", xlab="Distance (nt)", ylab="",
    cuts=1000, colorkey=list(col=bcols), col.regions=bcols,
    useRaster=T, aspect="fill", at=seq(0, 0.5, length.out=500),
    scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
);
grid.arrange(amap, cmap, widths=c(0.2, 0.8), padding=unit(0.2, 'line'));


edata = cbind(   Tot = log(rsum[pOrder[rows]]),
                 Esc = mTSB$ElongR[pOrder[rows]],
                 Loc = mTSB$Location[pOrder[rows]],
                 Cls =   pClass[pOrder[rows]]
        );
plotA = t(edata) / apply(edata, 2, max);
amap = levelplot( x=plotA,
    main="Reads", ylab="TSS", xlab="",
    cuts=1000, colorkey=F, col.regions=acols,
    aspect="fill", at=seq(0, 1, length.out=1000), useRaster=T,
    scales=list(
        x=list(labels=colnames(edata), rot=90, tck=c(1,0), cex=0.7),
        y=list(draw=F)
    )
);
cmap = levelplot( x=plotR,
    main="RppH", xlab="Distance (nt)", ylab="",
    cuts=1000, colorkey=list(col=bcols), col.regions=bcols,
    useRaster=T, aspect="fill", at=seq(0, 0.5, length.out=1000),
    scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
);
grid.arrange(amap, cmap, widths=c(0.2, 0.8), padding=unit(0.2, 'line'));


edata = cbind(  Tot = log(csum[pOrder[rows]]),
                Esc = mTSB$ElongR[pOrder[rows]],
                Loc = mTSB$Location[pOrder[rows]],
                Cls =   pClass[pOrder[rows]]
            );
plotA = t(edata) / apply(edata, 2, max);
amap = levelplot( x=plotA,
    main="Reads", ylab="TSS", xlab="",
    cuts=1000, colorkey=F, col.regions=acols,
    aspect="fill", at=seq(0, 1, length.out=1000), useRaster=T,
    scales=list(
        x=list(labels=colnames(edata), rot=90, tck=c(1,0), cex=0.7),
        y=list(draw=F)
    )
);
cmap = levelplot( x=plotC,
    main="Capped", xlab="Distance (nt)", ylab="",
    cuts=1000, colorkey=list(col=bcols), col.regions=bcols,
    useRaster=T, aspect="fill", at=seq(0, 0.5, length.out=1000),
    scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
);
grid.arrange(amap, cmap, widths=c(0.2, 0.8), padding=unit(0.2, 'line'));



edata = cbind(     Tot = log(usum[pOrder[rows]]),
                Esc = mTSB$ElongR[pOrder[rows]],
                Loc = mTSB$Location[pOrder[rows]],
                Cls =   pClass[pOrder[rows]]
        );
plotA = t(edata) / apply(edata, 2, max);
amap = levelplot( x=plotA,
    main="Reads", ylab="TSS", xlab="",
    cuts=1000, colorkey=F, col.regions=acols,
    aspect="fill", at=seq(0, 1, length.out=1000), useRaster=T,
    scales=list(
        x=list(labels=colnames(edata), rot=90, tck=c(1,0), cex=0.7),
        y=list(draw=F)
    )
);
cmap = levelplot( x=plotU,
    main="Uncapped", xlab="Distance (nt)", ylab="",
    cuts=1000, colorkey=list(col=bcols), col.regions=bcols,
    useRaster=T, aspect="fill", at=seq(0, 0.5, length.out=1000),
    scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
);
grid.arrange(amap, cmap, widths=c(0.2, 0.8), padding=unit(0.2, 'line'));


dev.off();
