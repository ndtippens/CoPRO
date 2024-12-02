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


TSS$normC = ifelse(strand(TSS) == "+", 1, -1);
pTRE = TREs[TREs$nTSS>0 & TREs$pTSS>0];
pTSS = subsetByOverlaps(TSS, pTRE);
pTSS = pTSS[ pTSS$Pause >= 1 ];
mTSB = mTSB[ mTSB$TSS %in% pTSS$ID5 ];
pTSB = mTSB;
unst = pTRE$nS == 0;
summary(factor(ifelse(unst, "Unstable", "Stable")));

bcols = colorRampPalette(c("dodgerblue4", "white", "firebrick"))(501);
ncols = colorRampPalette(c("white", "dodgerblue4"))(501);
pcols = colorRampPalette(c("white", "firebrick"))(501);
qcols = c('#B2222280', '#22222280', '#AAAAAA80', '#228B2280');
mkcols = c('#888888', '#B22222', 'steelblue4');


binSize = 5;
nBins = 201; # 1 kbp window
xat = 101 + c(-100, 0, 100);
xlabels = c(-1, 0, 1);

sTSB$normC = ifelse(strand(sTSB) == "+", sTSB$C, -sTSB$C);
sTSS = pTSS[ pTSS$Stability == 'Stable' ];

#parallel
nucpos = c(200,400,600);
nucw = c(200,200,200);


#strand(TSSpl) = "-";
allTSB = mTSB[ order(seqnames(mTSB), start(mTSB)) ];

coreTSS = promoters( allTSB, upstream=300, downstream=100);
sure = import.bw("~/Rmain/STARR/data/SuRE23.plasmid.norm.combined.45.55.plus.160504.bw", which=coreTSS);
suremn = import.bw("~/Rmain/STARR/data/SuRE23.plasmid.norm.combined.45.55.minus.160504.bw", which=coreTSS);
strand(sure) = '+';
strand(suremn) = '-';
sure = append(sure, suremn);
rm(suremn);
invisible(gc());
hits = findOverlaps( allTSB, sure );
allTSB$SuRE = 0;
allTSB$SuRE[hits@from] = sure$score[hits@to];

tdist = diff(start(allTSB), lag=1);
rootT = allTSB[ which(abs(tdist) <= 2000) ];
leafT = allTSB[ which(abs(tdist) <= 2000)+1 ];
maxT = c( rootT[rootT$C >= leafT$C], leafT[leafT$C > rootT$C] );
minT = c( leafT[rootT$C >= leafT$C], rootT[leafT$C > rootT$C] );
neg = strand(maxT) == "-";
tclass = ifelse( allTSB$ID5 %in% maxT$ID5, tblue, tblack );


maxTPause = GenomicRanges::shift( maxT, ifelse(maxT$ID5>0, maxT$Peak, -maxT$Peak) );
Tgap = pgap( maxTPause, minT, ignore.strand=T );
Tgap$neg = strand(maxT) == "-";
#Tgap$neg = rowMaxs( cbind( DP1pl$C, DP2pl$C, DP1mn$C, DP2mn$C ) ) > 2;
Tgap$Stable = maxT$Stability == "Stable";
Tgap$Promoter = maxT$ID5 %in% TX$maxTSS;
Tgap$type = ifelse( strand(maxT) == strand(minT), 'samestrand', 'convergent' );
Tgap$type[ which(strand(maxT) != strand(minT) & (maxT$ID5 > abs(minT$ID5) | abs(maxT$ID5) < minT$ID5)) ] = 'divergent';
#  intTSS = as.matrix(findOverlaps(Tgap, Tgap, type='within'));
#  intTSS = unique(intTSS[ intTSS[,1] != intTSS[,2], 1 ]);
pdf("../plots/TSS_distributions.pdf", height=4, width=3);

hist(width(Tgap)[isdiv], xlab="Distance (bp)", breaks=seq(0, 2050, 10), xlim=c(0, 1000), col='dodgerblue', border=F, main="Divergent");
hist(width(Tgap)[isconv], xlab="Distance (bp)", breaks=seq(0, 2050, 10), xlim=c(0, 1000), col='dodgerblue', border=F, main="Convergent");
hist(width(Tgap)[isss], xlab="Distance (bp)", breaks=seq(0, 2050, 10), xlim=c(0, 1000), col='dodgerblue', border=F, main="Same Strand");
hist(width(Tgap), xlab="Distance (bp)", breaks=seq(0, 2050, 10), xlim=c(0, 1000), col='dodgerblue', border=F, main="All TSSes");

dev.off();

pdf("../plots/CoPRO_vs_SuRE.pdf", width=4, height=5)
isdiv = Tgap$type == "divergent";
isconv = Tgap$type == "convergent";
isss = Tgap$type == "samestrand";

plot( maxT$SuRE[isdiv]+1, maxT$C[isdiv], log='xy', col=tblue, pch=18, main="maxTSS of Div Pair", xlim=c(1,500), ylim=c(1,500) );
plot( minT$SuRE[isdiv]+1, minT$C[isdiv], log='xy', col=tblack, pch=18, main="minTSS of Div Pair", xlim=c(1,500), ylim=c(1,500) );
plot( maxT$SuRE[isconv]+1, maxT$C[isconv], log='xy', col=tblue, pch=18, main="maxTSS of Conv Pair", xlim=c(1,500), ylim=c(1,500) );
plot( minT$SuRE[isconv]+1, minT$C[isconv], log='xy', col=tblack, pch=18, main="minTSS of Conv Pair", xlim=c(1,500), ylim=c(1,500) );
plot( maxT$SuRE[isss]+1, maxT$C[isss], log='xy', col=tblue, pch=18, main="maxTSS of SS Pair", xlim=c(1,500), ylim=c(1,500) );
plot( minT$SuRE[isss]+1, minT$C[isss], log='xy', col=tblack, pch=18, main="minTSS of SS Pair", xlim=c(1,500), ylim=c(1,500) );

dev.off();

for( i in 1:3 ) {
  bstart = nucpos[i]-nucw[i];
  bstop = nucpos[i]+nucw[i];
  bounds = width(Tgap) >= bstart & width(Tgap) <= bstop;
  msgout( bstart, '-', bstop, ':   ', sum(bounds) );
  Output[[i]] = Tgap[bounds];

  pclass = factor(ifelse(bounds, paste(i-1, 'Nucs'), 'Other'));
  #summary(pclass);
  PsVstab = chisq.test(pclass, Tgap$Stable);
  print(PsVstab);
  print(PsVstab$observed / PsVstab$expected);
}

Factors = c(
  'MNase-seq', 'MNase_H3K4me1', 'MNase_H3K122ac', 'MNase_H3K27ac',
  'RoadMap-H2A.Z', 'RoadMap-H3K4me2',  'RoadMap-H3K9ac',
  'RoadMap-H3K4me1', 'RoadMap-H3K4me3',  'RoadMap-H3K27ac',
  'PhastCons_100way'
);
Signal = list();
Stable = list();
Unstable = list();
for( f in Factors ) {
  Signal[[f]] = matrix(0, nrow=nBins, ncol=3);
  Stable[[f]] = matrix(0, nrow=nBins, ncol=3);
  Unstable[[f]] = matrix(0, nrow=nBins, ncol=3);
}

for( i in 1:3 ) {
  bstart = nucpos[i]-nucw[i];
  bstop = nucpos[i]+nucw[i];
  bounds = width(Tgap) >= bstart & width(Tgap) <= bstop;
  window = resize(Tgap[bounds], width=binSize*nBins, fix='center');

  for( f in Factors ) {
    bwf = BigWigFile( paste0('/media/nate/Diatom/Sequencing/ChIP/Encode/bw/', f, '.bw') );
    rawm = summary( bwf, window, size=nBins, type='mean', defaultValue=0, as='matrix' );
    rawm[ which(window$neg), ] = rawm[ which(window$neg), nBins:1 ];
    Signal[[f]][,i] = colMeans(rawm);
    Stable[[f]][,i] = colMeans(rawm[window$Stable,]);
    Unstable[[f]][,i] = colMeans(rawm[!window$Stable,]);
    #Signal[[f]][,i] = Signal[[f]][,i] / max(Signal[[f]][,i]);
  }
}

Capped = list();
hits = findOverlaps( maxT, Output[[1]], maxgap=300, ignore.strand=T );
Capped[[1]] = aggregate( maxT$C[hits@from] ~ hits@to, FUN=mean )[,2];
hits = findOverlaps( maxT, Output[[2]], maxgap=300, ignore.strand=T );
Capped[[2]] = aggregate( maxT$C[hits@from] ~ hits@to, FUN=mean )[,2];
hits = findOverlaps( maxT, Output[[3]], maxgap=300, ignore.strand=T );
Capped[[3]] = aggregate( maxT$C[hits@from] ~ hits@to, FUN=mean )[,2];
names(Capped) = c('None', 'One', 'Two');

pdf("../plots/Spacing_RedPairs.pdf", width=4, height=4);

matplot( ghist$breaks[-1], ghist$counts, type='l', lwd=2, lty=1,
  xlab='Distance to Convergent TSS (bp)', ylab='Count',
  ylim=c(10, 1000), xlim=c(0, 1100) );
#abline(v=nucpos, col=mkcols, lty=c(2,1,1), lwd=2);
#legend( 'topright', title='Nucleosomes', legend=c('None', 'One', 'Two'), lty=c(3,1,1), col=mkcols, lwd=2, bty='n');

boxplot( Capped, col=mkcols, main='', xlab='Internal Nucleosomes',
  ylab='maxTSN Occupancy', log='y', notch=T, outline=F );

for( f in Factors ) {
  matplot( x=seq(-500, 500, binSize), y=Signal[[f]], type='l', lty=c(2,1,1), lwd=2, xlab='Distance from Midpoint', ylab='', col=mkcols, main=f);
}
for( f in Factors ) {
  matplot( x=seq(-500, 500, binSize), y=Stable[[f]], type='l', lty=c(2,1,1), lwd=2, xlab='Distance from Midpoint', ylab='', col=mkcols, main=paste('Stable', f));
}
for( f in Factors ) {
  matplot( x=seq(-500, 500, binSize), y=Unstable[[f]], type='l', lty=c(2,1,1), lwd=2, xlab='Distance from Midpoint', ylab='', col=mkcols, main=paste('Unstable', f));
}

dev.off();
