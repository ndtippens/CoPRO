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

binSize = 10;
nBins = 201; # 1 kbp window

eTRE = TREs[ width(TREs) <= 2010 ];
unst = eTRE$nS == 0;
eTRE$bwidth = floor(width(eTRE)/binSize);
eTRE = eTRE[ order( -unst, -eTRE$bwidth, -eTRE$Cap ) ];

bcols = colorRampPalette(c("white", "#27003B"))(500);
bcols=c(bcols, '#999999');

fnames = list.files('/media/nate/Diatom/Sequencing/ChIP/Encode/bw');
pnames = list.files('../plots/ChIP_inTIDs/');

fTRE = resize( eTRE, width=nBins*binSize, fix='center' );
strand(fTRE) = "*";
xat=c(1, 101, 201);
xlabels=c(-1, 0, 1);
cksz = floor(length(eTRE)/200);
TIDw = floor(width(eTRE) / 20);
TIDst = cbind( 101-TIDw, 1:length(TIDw) );
TIDed = cbind( 101+TIDw, 1:length(TIDw) );
sTIDx = round( (1:200 - 0.5) * cksz );
sTIDst = cbind( TIDst[ sTIDx, 1 ], 1:200 );
sTIDed = cbind( TIDed[ sTIDx, 1 ], 1:200 );

mTSB = mTSB[ order(-mTSB$PsOcc) ];
hits = findOverlaps( mTSB, eTRE );
maxT = mTSB[  hits@from[!duplicated(hits@to)] ];
hits = findOverlaps( maxT, eTRE );
maxT = maxT[ order(hits@to) ];
revTID = which(strand(maxT) == "-");
bordern = floor(199*mean(unst))+1;

for( fn in fnames ) {
	fn = gsub( '.bw', '', fn );
	if( paste0(fn, '.pdf') %in% pnames )
		next;

	print( fn );
	bwf = BigWigFile( paste0('/media/nate/Diatom/Sequencing/ChIP/Encode/bw/', fn, '.bw') );
	smat = abs(summary( bwf, fTRE, size=nBins, type='mean', defaultValue=0, as='matrix' ));
	smat[ revTID, ] = smat[ revTID, nBins:1 ];

	pdf(paste0("../plots/ChIP_inTIDs/", fn, ".pdf"), width=3, height=6, useDingbats=F);

	meta = LinearHeatmap( smat, 200, 1:nrow(smat), FUN='mean' )
	meta = t(meta);
	lmeta = log2(meta+1);
	cscale = seq(min(meta), max(meta), length.out=500);
	lscale = seq(min(lmeta), max(lmeta), length.out=500);
	cscale = c(cscale, max(meta)+1);
	lscale = c(lscale, max(lmeta)+1);
	meta[sTIDst] = max(meta)+1;
	meta[sTIDed] = max(meta);
	lmeta[sTIDst] = max(lmeta)+1;
	lmeta[sTIDed] = max(lmeta);
	meta[,bordern] = max(meta);
	lmeta[,bordern] = max(lmeta);

	print( levelplot( x=meta,
		main=paste(fn, "at TIDs"), xlab="Distance from Center (kb)", ylab="Unstable                                      Stable",
		colorkey=F, col.regions=bcols,
		useRaster=T, aspect="fill", at=cscale,
		scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
	) );

	print( levelplot( x=lmeta,
		main=paste(fn, "at TIDs"), xlab="Distance from Center (kb)", ylab="Unstable                                      Stable",
		colorkey=F, col.regions=bcols,
		useRaster=T, aspect="fill", at=lscale,
		scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
	) );

	dev.off();
}
