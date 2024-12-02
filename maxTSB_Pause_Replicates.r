#!/usr/bin/Rscript

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicAlignments)
	library(GenomicFiles)
	library(data.table)
	library(rtracklayer)
	library(GenomicRanges)
});
source("./CoPRO_Functions.r");



load("../data/newK562_CoPROcalls.Rdata");

maxTSS = as.data.table( mcols(sTSB) );
# determine which site has max capped reads at each TSS
maxTSS[, isMax := (which.max(C) == seq_len(.N)), by=TSS ];
mTSB = sTSB[ maxTSS[,isMax] ];
mTSB = mTSB[ order(mTSB$ID5) ];

mTSB = mTSB[ mTSB$C >= 2 & mTSB$R >= 2 ];
TSSid = mTSB$ID5;
nTSS = length(mTSB);
print(nTSS);

load("../data/Mappable19mers.Rdata");
load("../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Rep1.Rdata");
Rep1 = Pooled[Pooled$ID5 %in% MapID5];
load("../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Rep2.Rdata");
Rep2 = Pooled[Pooled$ID5 %in% MapID5];

windowSize = 100;
Rep1 = Rep1[ width(Rep1) >= 20 & width(Rep1) <= windowSize & Rep1$ID5 %in% TSSid ];
Rep2 = Rep2[ width(Rep2) >= 20 & width(Rep2) <= windowSize & Rep2$ID5 %in% Rep1$ID5 ];
RepIDs = TSSid[ TSSid %in% Rep1$ID5 ];
Nrep = length(RepIDs);
invisible(gc());

cmat1 = matrix(data=0, nrow=Nrep, ncol=windowSize);
rownames(cmat1) = RepIDs;
colnames(cmat1) = as.character(1:windowSize);
rmat1 = cmat1;
cmat2 = cmat1;
rmat2 = cmat1;
emat1 = cmat1;
emat2 = cmat1;
cmat1[ cbind(as.character(Rep1$ID5), width(Rep1)) ] = Rep1$C;
rmat1[ cbind(as.character(Rep1$ID5), width(Rep1)) ] = Rep1$R;
emat1[ cbind(as.character(Rep1$ID5), width(Rep1)) ] = Rep1$C + Rep1$U + Rep1$R;
cmat2[ cbind(as.character(Rep2$ID5), width(Rep2)) ] = Rep2$C;
rmat2[ cbind(as.character(Rep2$ID5), width(Rep2)) ] = Rep2$R;
emat2[ cbind(as.character(Rep2$ID5), width(Rep2)) ] = Rep2$C + Rep2$U + Rep2$R;

# only test sites with at least 12 reads in Rep1
cmat2 = cmat2[ rowSums(emat1) >= 12, ];
rmat2 = rmat2[ rowSums(emat1) >= 12, ];
emat2 = emat2[ rowSums(emat1) >= 12, ];
emat1 = emat1[ rowSums(emat1) >= 12, ];

cmat1 = cmat1/rowSums(cmat1);
rmat1 = rmat1/rowSums(rmat1);
emat1 = emat1/rowSums(emat1);
cmat2 = cmat2/rowSums(cmat2);
rmat2 = rmat2/rowSums(rmat2);
emat2 = emat2/rowSums(emat2);

# pause classes from each replicate
PsClasses=c('Early', 'Late', 'Both');
SCS = function(x) { sum(cumsum(x)) };
classifyPause = function( rdmat ) {
	totPs = rowSums(rdmat[,20:60]);
	totPs[!totPs] = 1;
	fearly = rowSums(rdmat[,20:32])/totPs;
	flate = rowSums(rdmat[,33:60])/totPs;
	#cls = ifelse(fearly > flate, 1, 2);
	#cls[ cls == 1 & fearly < 0.6 ] = 3;
	#cls[ cls == 2 &  flate < 0.6 ] = 3;
	#nrmat = rdmat[,20:60]/totPs;
	#pscore = apply(nrmat, 1, SCS);
	cutoffs = quantile(flate, na.rm=T);
	cls = ifelse( flate <= cutoffs[2], 2, 3 );
	cls[ flate >= cutoffs[4] ] = 1;
	return(cls);
}

pclass1 = classifyPause( emat1 );
pclass2 = classifyPause( emat2 );
#pclass3 = classifyPause( (cmat2+rmat2) );

for( c2 in 1:3 ) {
	for( c1 in 1:3 ) {
		ct1 = sum(pclass1==c1 & pclass2==c2, na.rm=T);
		tot1 = round(100.0*ct1/sum(pclass2==c2, na.rm=T), digits=2);
		#ct3 = sum(pclass3==c1 & pclass2==c2, na.rm=T);
		#tot3 = round(100*ct3/sum(pclass2==c2, na.rm=T));
		out = c(PsClasses[c2], PsClasses[c1], ct1, paste0(tot1, '%') ); #, ct3, paste0(tot3, '%'));
		writeLines( out, sep="\t" );
		msgout();
	}
}
