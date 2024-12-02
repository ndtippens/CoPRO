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


pTSB = mTSB[ match(mTSB$TSS, TSS$ID5) ];
eTRE = TREs[ width(TREs) <= 2010 & TREs$pTSS & TREs$nTSS ];
unst = eTRE$nS == 0;
summary(factor(ifelse(unst, "Unst", "Stable")));

TSS = TSS[ order(-TSS$Pause) ];
hits = findOverlaps( TSS, eTRE );
# find first & second max TSS
maxt = hits@from[!duplicated(hits@to)];
TSS1 = TSS[  maxt ];
TSS  = TSS[ -maxt ];
hits = findOverlaps( TSS1, eTRE );
TSS1 = TSS1[ order(hits@to) ]; # now 1:1 correspondence with eTRE ordering

# look for TSS2 on opposite strand only
strand(eTRE) = ifelse(strand(TSS1) == "+", "-", "+");
hits = findOverlaps( TSS, eTRE );
maxt = hits@from[!duplicated(hits@to)];
TSS2 = TSS[  maxt ];
strand(eTRE) = "*";

hits = findOverlaps( mTSB, TSS1 );
TSB1 = mTSB[ hits@from ];
hits = findOverlaps( mTSB, TSS2 );
TSB2 = mTSB[ hits@from ];

hits = findOverlaps( TSS2, eTRE );
TSS2 = TSS2[ order(hits@to) ];
hits = findOverlaps( TSB1, eTRE );
TSB1 = TSB1[ order(hits@to) ];
hits = findOverlaps( TSB2, eTRE );
TSB2 = TSB2[ order(hits@to) ];

bcols = colorRampPalette(c("white", "#27003B"))(500);
bcols=c(bcols, '#999999');

binSize = 10;
nBins = 201; # 1 kbp window
fTRE = resize( TSB1, width=nBins*binSize, fix='center' );
strand(fTRE) = "*";
xat=c(1, 101, 201);
xlabels=c(-1, 0, 1);
bordern = floor(199*mean(unst))+1;

tOff = ifelse(strand(TSB1) == '+',   end(TSB2)-start(TSB1), end(TSB1)-start(TSB2) );

fnames = list.files('/media/nate/Diatom/Sequencing/ChIP/Encode/bw');
pnames = list.files('../plots/ChIP_maxTSS/');

for( fn in fnames ) {
	fn = gsub( '.bw', '', fn );
	if( paste0(fn, '.pdf') %in% pnames )
		next;


	bwf = BigWigFile( paste0('/media/nate/Diatom/Sequencing/ChIP/Encode/bw/', fn, '.bw') );
	smat = abs(summary( bwf, fTRE, size=nBins, type='mean', defaultValue=0, as='matrix' ));
	# flip orientation on minus strand
	smat[ which(strand(TSS1) == '-'), ] = smat[ which(strand(TSS1) == '-'), nBins:1 ];

	meta = LinearHeatmap( smat, 200, order(-unst, tOff), FUN='mean' );
	meta = t(meta);
	lmeta = log2(1+meta);
	cscale = seq(min(meta), max(meta), length.out=500);
	lscale = seq(min(lmeta), max(lmeta), length.out=500);
	cscale = c(cscale, max(meta)+1);
	lscale = c(lscale, max(lmeta)+1);

	meta[,bordern] = max(meta)+1;
	lmeta[,bordern] = max(lmeta)+1;


	pdf(paste0("../plots/ChIP_maxTSS/", fn, ".pdf"), width=3, height=6, useDingbats=F);

	print( levelplot( x=meta,
		main=paste(fn, "at maxTSS pairs"), xlab="Distance from Max TSS (kb)", ylab="Unstable                                           Stable",
		cuts=500, colorkey=F, col.regions=bcols,
		useRaster=T, aspect="fill", at=cscale,
		scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
	) );

	print( levelplot( x=lmeta,
		main=paste(fn, "at maxTSS pairs"), xlab="Distance from Max TSS (kb)", ylab="Unstable                                           Stable",
		cuts=500, colorkey=F, col.regions=bcols,
		useRaster=T, aspect="fill", at=lscale,
		scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
	) );

	dev.off();
}
