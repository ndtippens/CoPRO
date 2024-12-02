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
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(fpc))
suppressPackageStartupMessages(library(scatterplot3d))
registerDoParallel(cores=7);
source("./CoPRO_Functions.r");

load("../data/newK562_CoPRO_maxTSB.Rdata");
load("../data/newK562_CoPROcalls.Rdata");

# require 21mer mappability (20 nt + 1 run-on)
load("../data/Mappable21mers.Rdata");
mTSB = mTSB[ mTSB$ID5 %in% MapID5 ];

mTSB = mTSB[ mTSB$C >= 2 & mTSB$U >= 0.5 & mTSB$EscOcc > 0 ];
TSSid = mTSB$ID5;
nTSS = length(mTSB);
print(nTSS);

windowSize=500;
pReads = TSB[ width(TSB) >= 20 & width(TSB) <= windowSize & TSB$ID5 %in% TSSid ];

cmat = matrix(data=0, nrow=nTSS, ncol=windowSize);
rownames(cmat) = TSSid;
colnames(cmat) = as.character(1:windowSize);
tmat=cmat;
cmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$C + 0.12;
tmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$C + pReads$U + 0.12;

  capt = apply(cmat, 1, cumsum);
  allt = apply(tmat, 1, cumsum);
cratio = 100*t( capt / allt );
  LCap = rowSums(cmat[,60:windowSize]);
fincap = LCap/(rowSums(tmat[,60:windowSize]));

pdf("../plots/CappingKinetics.pdf", width=3, height=3);
pcols=c('dodgerblue4', 'forestgreen');

pidx = cbind(as.character(mTSB$ID5), mTSB$Peak);
fcap = cmat[ pidx ] / tmat[ pidx ];
boxplot( 100*fcap ~ mTSB$Peak, xlim=c(0,40), ylim=c(0,100),
	col='gray', lty=1, lwd=1, outline=F, notch=F,
	main='Capping Yield', xlab='Pause Position', ylab="Percent Capped"
);


boxplot( cratio[,20:60], xlim=c(0,40), ylim=c(0,100),
	col='gray', lty=1, lwd=1, outline=F, notch=T,
	main='Capping Kinetics', xlab='Distance from TSB (nt)', ylab="Percent Capped"
);

Efrac = colSums(cmat[mTSB$PsClass=='Early',])/colSums(tmat[mTSB$PsClass=='Early',]);
Lfrac = colSums(cmat[mTSB$PsClass=='Late',])/colSums(tmat[mTSB$PsClass=='Late',]);
matplot( 100*cbind(Efrac, Lfrac), ylim=c(60,100), xlim=c(20,60), type='l',
	col=pcols, lty=1, lwd=1, main='Capping vs Length', xlab='Distance from maxTSB', ylab="Percent Capped"
);
legend("bottomright", legend=c('Early', 'Late'), lty=1, col=pcols, lwd=2, bty='n')

Efrac = colSums(cmat[mTSB$Stability=='Unstable',])/colSums(tmat[mTSB$Stability=='Unstable',]);
Lfrac = colSums(cmat[mTSB$Stability=='Stable',])/colSums(tmat[mTSB$Stability=='Stable',]);
matplot( 100*cbind(Efrac, Lfrac), ylim=c(60,100), xlim=c(20,60), type='l',
	col=pcols, lty=1, lwd=1, main='Capping vs Length', xlab='Distance from maxTSB', ylab="Percent Capped"
);
legend("bottomright", legend=c('Unstable', 'Stable'), lty=1, col=pcols, lwd=2, bty='n')


boxplot( rowMaxs(cratio[,61:windowSize], na.rm=T) ~ mTSB$PsClass, ylim=c(60,100), notch=T, outline=F,
	col=c('gray', pcols), main='Capping Yield', xlab='', ylab="% RNA Capped by 60 nt"
);

boxplot( rowMaxs(cratio[,61:windowSize], na.rm=T) ~ mTSB$Stability, ylim=c(60,100), notch=T, outline=F,
	col=c('gray', pcols), main='Capping Yield', xlab='', ylab="% RNA Capped by 60 nt"
);



#PrePause = list();
#InPause = list();
#PostPause = list();
#for( i in 22:50 ) {
#	nonPaused = (i < mTSB$Peak-5) & (capt[,i] < allt[,i]);
#	PrePause[[i]] = unname(cratio[nonPaused,i]);
#}
#boxplot( CapRate, xlim=c(22,50), ylim=c(0,100),
#	col='gray', lty=1, lwd=1, outline=F, notch=T,
#	main='Capping Kinetics', xlab='RNA Length', ylab="Percent Capped"
#);


dev.off();
