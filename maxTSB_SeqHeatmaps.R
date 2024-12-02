#!/usr/bin/Rscript
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFiles))
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

setwd("~/Rmain/CoPRO/CoPRO/")
source("./CoPRO_Functions.r");
genome = "~/bin/genomes/hg19/hg19.2bit"
registerDoParallel(cores=7);


load("../data/newK562_CoPRO_maxTSB.Rdata");
load("../data/newK562_CoPROcalls.Rdata");


eTRE = TREs[ TREs$nTSS & TREs$pTSS & TREs$Cap >= 4 ];
unst = eTRE$nS == 0;

aTSB = mTSB[ order(-mTSB$C) ];
hits = findOverlaps( aTSB, eTRE );
# find first & second max TSBs
maxt = hits@from[!duplicated(hits@to)];
TSB1 = mTSB[  maxt ];
aTSB  = aTSB[ -maxt ];
hits = findOverlaps( aTSB, eTRE );
maxt = hits@from[!duplicated(hits@to)];
TSB2 = aTSB[  maxt ];

hits = findOverlaps( TSB1, eTRE );
TSB1 = TSB1[ order(hits@to) ];
hits = findOverlaps( TSB2, eTRE );
TSB2 = TSB2[ order(hits@to) ];

# plot relative to the initiation base (and sorted by distance to the max pause)
pTSB = c(TSB1, TSB2);
pTSB = pTSB[ order(-pTSB$Peak) ];
TSBprom = promoters( pTSB, upstream=100, downstream=101 );
Slen = width(TSBprom)[1];
seqs = import.2bit(genome, which=TSBprom);
seqs[ strand(pTSB) == '-' ] = reverseComplement( seqs[ strand(pTSB) == '-' ] );
gmat = foreach( s=iter(seqs), .combine='cbind' ) %dopar% {
	unname( letterFrequencyInSlidingView( s, view.width=1, letters='GC' ) );
}
CpGm = foreach( s=iter(seqs), .combine='cbind' ) %dopar% {
	matchC = unname( letterFrequencyInSlidingView( s, view.width=1, letters='C' ) );
	matchG = unname( letterFrequencyInSlidingView( s, view.width=1, letters='G' ) );
	return( matchC[ 1 : (Slen-1) ] + matchG[ 2:Slen ] );
}
CpAm = foreach( s=iter(seqs), .combine='cbind' ) %dopar% {
	matchC = unname( letterFrequencyInSlidingView( s, view.width=1, letters='C' ) );
	matchA = unname( letterFrequencyInSlidingView( s, view.width=1, letters='A' ) );
	return( matchC[ 1 : (Slen-1) ] + matchA[ 2:Slen ] );
}
GpAm = foreach( s=iter(seqs), .combine='cbind' ) %dopar% {
	matchG = unname( letterFrequencyInSlidingView( s, view.width=1, letters='G' ) );
	matchA = unname( letterFrequencyInSlidingView( s, view.width=1, letters='A' ) );
	return( matchG[ 1 : (Slen-1) ] + matchA[ 2:Slen ] );
}

gmat = 100*gmat;
CpGm = ifelse(CpGm < 2, 0, 100);
CpGm = rbind( mean(CpGm), CpGm, mean(CpGm) );
CpAm = ifelse(CpAm < 2, 0, 100);
CpAm = rbind( mean(CpAm), CpAm, mean(CpAm) );
GpAm = ifelse(GpAm < 2, 0, 100);
GpAm = rbind( mean(GpAm), GpAm, mean(GpAm) );


#export(seqs[pTSB$Peak < 32], "../out/K562_mTSB_Early_seq5.fasta.gz")
#export(seqs[pTSB$Peak > 36], "../out/K562_mTSB_Late_seq5.fasta.gz")


fTRE = resize( pTSB, width=210, fix='center' );
strand(fTRE) = "*";


# permute by max pause dist
byPauseD = matrix(0, ncol=55, nrow=width(TSBprom)[1]);
for( pausep in 20:55) {
	byPauseD[,pausep] = rowMeans(gmat[,pTSB$Peak==pausep]);
}

pdf(file="../plots/maxTSB_SeqHeatmaps.pdf", width=2, height=3, useDingbats=F);

xats  = 101 + c(-100, 0, 100);
xlabs = c(-100, 'maxTSB', '+100');
plotN = 1000;

gcols = colorRampPalette(c("white", "steelblue3"))(501);
bcols = colorRampPalette(c("firebrick", "white", "steelblue3"))(501);
#bcols = c(bcols, 'black');

cbreaks = seq(0, 100, length.out=501);
levelplot( x=byPauseD[,55:20],
	main="GC Content", xlab="Position (nt)", ylab="",
	colorkey=list(col=bcols, at=cbreaks), col.regions=bcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

# linear walk, 100 steps
linMat = foreach( chunk = iter(gmat, by='column', chunksize=floor(ncol(gmat)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GC Content", xlab="Position (nt)", ylab="",
	colorkey=list(col=bcols, at=cbreaks), col.regions=bcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);


linMat = foreach( chunk = iter(gmat[,order(-pTSB$LateOcc)], by='column', chunksize=floor(ncol(gmat)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GC by Late", xlab="Position (nt)", ylab="",
	colorkey=list(col=bcols, at=cbreaks), col.regions=bcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

linMat = foreach( chunk = iter(gmat[,order(-pTSB$PsOcc)], by='column', chunksize=floor(ncol(gmat)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GC by Pause", xlab="Position (nt)", ylab="",
	colorkey=list(col=bcols, at=cbreaks), col.regions=bcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

linMat = foreach( chunk = iter(gmat[,order(-pTSB$Stability)], by='column', chunksize=floor(ncol(gmat)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GC by Stability", xlab="Position (nt)", ylab="",
	colorkey=list(col=bcols, at=cbreaks), col.regions=bcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);


cbreaks = seq(0, 50, length.out=501);
linMat = foreach( chunk = iter(CpGm, by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="CpG Content", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

cbreaks = seq(0, 30, length.out=501);
linMat = foreach( chunk = iter(CpGm[,order(-pTSB$LateOcc)], by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="CpG by Late", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

linMat = foreach( chunk = iter(CpGm[,order(-pTSB$PsOcc)], by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="CpG by Pause", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);


cbreaks = seq(0, 80, length.out=501);
linMat = foreach( chunk = iter(CpAm, by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="CpA Content", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

cbreaks = seq(0, 30, length.out=501);
linMat = foreach( chunk = iter(CpAm[,order(-pTSB$LateOcc)], by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="CpA by Late", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

linMat = foreach( chunk = iter(CpAm[,order(-pTSB$PsOcc)], by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="CpA by Pause", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);



cbreaks = seq(0, 50, length.out=501);
linMat = foreach( chunk = iter(GpAm, by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GpA Content", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

cbreaks = seq(0, 30, length.out=501);
linMat = foreach( chunk = iter(GpAm[,order(-pTSB$LateOcc)], by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GpA by Late", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

linMat = foreach( chunk = iter(GpAm[,order(-pTSB$PsOcc)], by='column', chunksize=floor(ncol(CpGm)/100)), .combine='cbind' ) %dopar% {
	rowMeans(chunk);
}
levelplot( x=linMat,
	main="GpA by Pause", xlab="Position (nt)", ylab="",
	colorkey=list(col=gcols, at=cbreaks), col.regions=gcols,
	useRaster=T, aspect="fill", at=cbreaks,
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);



dev.off();
