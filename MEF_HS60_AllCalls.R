#!/usr/bin/Rscript

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFiles))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
source("./CoPRO_Functions.r");
genomef = "~/bin/genomes/mm10/mm10.2bit";

mm10chr = paste0( 'chr', c(1:19, 'X', 'Y') );
mm10len   = c(  195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459,
				129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244,
				104043685,  98207768,  94987271,  90702639,  61431566, 171031299,  91744698
			 );


compute_mmIDs = function(reads) {
	chrOffsets = cumsum(as.numeric(c(0, mm10len)));
	chrOffsets = chrOffsets[1:length(chrOffsets)-1];
	names(chrOffsets) = mm10chr;
	offset = chrOffsets[as.character(seqnames(reads))];
	names(offset) = NULL;
	strdir = ifelse(strand(reads) == '+', 1, -1);
	RNA5 = ifelse( strdir>0, start(reads),   end(reads) );
	RNA3 = ifelse( strdir>0,   end(reads), start(reads) );
	reads$ID5 = strdir * (offset + RNA5);
	reads$ID3 = strdir * (offset + RNA3);
	return(reads);
}


load("../data/Quake_ReadWeights.Rdata");
load("~/Rmain/DualMap/data/mm10/RppHp-NHS-PROcap.Rdata");
NHSc = short;
load("~/Rmain/DualMap/data/mm10/RppHm-NHS-PROcap.Rdata");
NHSu = short;
load("~/Rmain/DualMap/data/mm10/RppHp-HS60-PROcap.Rdata");
 HSc = short;
load("~/Rmain/DualMap/data/mm10/RppHm-HS60-PROcap.Rdata");
 HSu = short;

NCnorm = 1;#sum(NHSc$count)*1E-6;
NUnorm = 1;#sum(NHSu$count)*1E-6;
HSnorm = 1;#sum(HSc$count)*1E-6;
HUnorm = 1;#sum(HSu$count)*1E-6;

NHSc$NC = round(NHSc$count/NCnorm, 3);
NHSc$NU = 0L;
NHSc$HC = 0L;
NHSc$HU = 0L;

NHSu$NC = 0L;
NHSu$NU = round(NHSu$count/NUnorm, 3);
NHSu$HC = 0L;
NHSu$HU = 0L;

HSc$NC = 0L;
HSc$NU = 0L;
HSc$HC = round(HSc$count/HSnorm, 3);
HSc$HU = 0L;

HSu$NC = 0L;
HSu$NU = 0L;
HSu$HC = 0L;
HSu$HU = round(HSu$count/HUnorm, 3);

NHSc$count = NULL;
NHSu$count = NULL;
HSc$count = NULL;
HSu$count = NULL;
Pooled = NHSc;

dupes = findOverlaps(NHSu, Pooled, type='equal');
Pooled[dupes@to]$NU = NHSu[dupes@from]$NU;
Pooled = append(Pooled, NHSu[-dupes@from]);

dupes = findOverlaps(HSc, Pooled, type='equal');
Pooled[dupes@to]$HC = HSc[dupes@from]$HC;
Pooled = append(Pooled, HSc[-dupes@from]);

dupes = findOverlaps(HSu, Pooled, type='equal');
Pooled[dupes@to]$HU = HSu[dupes@from]$HU;
Pooled = append(Pooled, HSu[-dupes@from]);

rm(NHSc, NHSu, HSc, HSu);
invisible(gc());

Pooled = compute_mmIDs(Pooled);
Pooled=Pooled[order(Pooled$ID5)];
save(Pooled, file="../data/MEF_HS60_Pooled.Rdata");

load("../data/MEF_HS60_Pooled.Rdata");
path = paste0("../data/mm10/CoPRO_", 'NHScap');
#CoPROtoBigWig( Pooled[ mcols(Pooled)[,'NC'] > 0, ], path, 'NC' );
path = paste0("../data/mm10/CoPRO_",  'HScap');
#CoPROtoBigWig( Pooled[ mcols(Pooled)[,'HC'] > 0, ], path, 'HC' );
path = paste0("../data/mm10/CoPRO_", 'NHSbg');
#CoPROtoBigWig( Pooled[ mcols(Pooled)[,'NU'] > 0, ], path, 'NU' );
path = paste0("../data/mm10/CoPRO_",  'HSbg');
#CoPROtoBigWig( Pooled[ mcols(Pooled)[,'HU'] > 0, ], path, 'HU' );

ReadTable = as.data.table(mcols(Pooled));
setkey(ReadTable, ID5);
ReadTable[, width := abs(ID5-ID3)+1];

starts = ReadTable[ width >= 18, .(PsPos=sum(HC>0|NC>0), HS=sum(HC), NHS=sum(NC), HU=sum(HU), NU=sum(NU)), by=ID5 ];
# TSS = at least 5 3' ends
starts = starts[ PsPos >= 5, ];
ReadTable = ReadTable[ ID5 %in% starts$ID5, ];
starts = round( ReadTable[, SigR:=max(HC/(HC+HU), NC/(NC+NU), na.rm=T), by=ID5 ], 2 );
TSB    = Pooled[ Pooled$ID5 %in% starts[ SigR>0.9, ID5 ] ];
rm(ReadTable, starts, Pooled);
invisible(gc());

# sum over 3' ends to simplify TSBs
sTSB = TSB[ !duplicated(TSB$ID5) ];
# collapse onto 5' base
sTSB = resize( sTSB, width=1 );
sTSB$ID3 = NULL;

sTSB$HC = aggregate( HC ~ ID5, data=TSB, FUN=sum )[,2];
sTSB$NC = aggregate( NC ~ ID5, data=TSB, FUN=sum )[,2];
sTSB$HU = aggregate( HU ~ ID5, data=TSB, FUN=sum )[,2];
sTSB$NU = aggregate( NU ~ ID5, data=TSB, FUN=sum )[,2];
HsigR = with(sTSB, HC/(HC+HU+1));
NsigR = with(sTSB, NC/(NC+NU+1));

sTSB = sTSB[ HsigR>0.9 | NsigR>0.9 ];
TSB  =  TSB[ TSB$ID5 %in% sTSB$ID5 ];

msgout(pn(length(sTSB)), "TSB calls (at least 5 different 3' ends, >90% Cap-specific reads)");


HS.TRE = compute_TREs( sTSB[sTSB$HC>0], 600 );
NHS.TRE = compute_TREs( sTSB[sTSB$NC>0], 600 );
HS.TRE = HS.TRE[width(HS.TRE)>=10];
NHS.TRE = NHS.TRE[width(NHS.TRE)>=10];

 HS.TSB = subsetByOverlaps(sTSB[sTSB$HC>0],  HS.TRE);
NHS.TSB = subsetByOverlaps(sTSB[sTSB$NC>0], NHS.TRE);
TSB  =  TSB[  TSB$ID5 %in% c(HS.TSB$ID5, NHS.TSB$ID5) ];
sTSB = sTSB[ sTSB$ID5 %in% c(HS.TSB$ID5, NHS.TSB$ID5) ];

#PsBs = GenomicRanges::resize(TSB[1:(10^5)], width=1, fix='end');
#flank3 = GenomicRanges::promoters(PsBs, upstream=6, downstream=5); # skip run-on base
#seq3 = import.2bit(genomef, which=flank3);
#seq3[ strand(flank3) == '-' ] = reverseComplement(seq3[ strand(flank3) == '-' ]);
#export(seq3, "../out/MEFs_All3p_seqs.fasta.gz");


 HS.TSS = compute_TSS( HS.TSB, 60 );
NHS.TSS = compute_TSS( NHS.TSB, 60 );
hits=findOverlaps( HS.TSS, NHS.TSS );
msgout(pn(length(HS.TSS)), "HS60 TSS calls");
msgout(pn(length(NHS.TSS)), "NHS TSS calls");
msgout(pn(length(hits)), "overlapping TSS");
msgout(pn(length(HS.TSS[-unique(hits@from)])), "HS-specific TSS");
msgout(pn(length(NHS.TSS[-unique(hits@to)])), "NHS-specific TSS");

hits=findOverlaps( HS.TRE, NHS.TRE );
msgout(pn(length(HS.TRE)), "HS TREs");
msgout(pn(length(NHS.TRE)), "NHS TREs");
msgout(pn(length(hits)), "overlapping TREs");
msgout(pn(length(HS.TRE[-unique(hits@from)])), "HS-specific TRE");
msgout(pn(length(NHS.TRE[-unique(hits@to)])), "NHS-specific TRE");

# assign TSB to TSS
hits = findOverlaps(  HS.TSB,  HS.TSS, type="any" );
HS.TSB$TSS = HS.TSS$ID5[ hits@to ];
HS.TSS$C[unique(hits@to)]   = aggregate( HC ~ TSS, data=HS.TSB, FUN=sum )[,2];
hits = findOverlaps( NHS.TSB, NHS.TSS, type="any" );
NHS.TSB$TSS = NHS.TSS$ID5[ hits@to ];
NHS.TSS$C[unique(hits@to)]   = aggregate( NC ~ TSS, data=NHS.TSB, FUN=sum )[,2];

# assign TSS to TRE
hits=findOverlaps( HS.TSS, HS.TRE );
HS.TSS$TID[hits@from] = HS.TRE$ID5[hits@to];

hits=findOverlaps( NHS.TSS, NHS.TRE );
NHS.TSS$TID[hits@from] = NHS.TRE$ID5[hits@to];

maxTSS = as.data.table( mcols( HS.TSB) );
# determine which site has max capped reads at each TID
maxTSS[, isMax := (which.max(HC) == seq_len(.N)), by=TSS ];
HS.mSB = HS.TSB[ maxTSS[,isMax] ];
HS.mSB = HS.mSB[ order(HS.mSB$ID5) ];

maxTSS = as.data.table( mcols(NHS.TSB) );
# determine which TSS has max capped reads at each TID
maxTSS[, isMax := (which.max(NC) == seq_len(.N)), by=TSS ];
NHS.mSB = NHS.TSB[ maxTSS[,isMax] ];
NHS.mSB = NHS.mSB[ order(NHS.mSB$ID5) ];

hits = findOverlaps( HS.mSB, HS.TSS );
HS.TSS$mTSB[hits@to] = HS.mSB[hits@from]$ID5;
hits = findOverlaps( NHS.mSB, NHS.TSS );
NHS.TSS$mTSB[hits@to] = NHS.mSB[hits@from]$ID5;

save(TSB, sTSB, HS.TSB, NHS.TSB, HS.TSS, NHS.TSS, HS.TRE, NHS.TRE, file='../data/mm10/MEF_AllCalls.Rdata');

qHeatmap = function( path, regions, rindex ) {
	xat=c(1, 101, 201);
	xlabels=c(-1000, 0, 1000);
	bcols = colorRampPalette(c("white", "deepskyblue4"))(500);

	bwf = BigWigFile( path );
	smat = abs(summary( bwf, regions, size=nBins, type='mean', defaultValue=0, as='matrix' ));
	meta = smat[ order( rindex ), ];
	meta = t(log2(meta+0.25));
	print( levelplot( x=meta,
		main=path, xlab="Distance to Center (nt)", ylab="",
		cuts=500, colorkey=list(col=bcols), col.regions=bcols,
		useRaster=T, aspect="fill", at=seq(min(meta), max(meta), length.out=500),
		scales=list(x=list(at=xat, tck=c(1,0), lab=xlabels), y=list(draw=F))
	) );
	return(smat);
}

qData = function( path, regions, rindex ) {
	bwf = BigWigFile( path );
	smat = abs(summary( bwf, regions, size=nBins, type='mean', defaultValue=0, as='matrix' ));
	return( smat );
}
