#!/usr/bin/Rscript

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicAlignments)
	library(GenomicFiles)
	library(data.table)
	library(rtracklayer)
	library(GenomicRanges)
	library(ggplot2)
	library(scales)
	library(viridis)
	library(lattice)
	library(doParallel)
});
register(MulticoreParam(8));
source("./CoPRO_Functions.r");

load("../data/newK562_CoPROcalls.Rdata");
load("../data/newK562_CoPRO_maxTSB.Rdata");
load("../data/Mappable19mers.Rdata");
mTSB = mTSB[ mTSB$ID5 %in% MapID5 ];

mTSB = mTSB[ mTSB$PsOcc >= 1 & mTSB$R >= 1 ];
TSSid = mTSB$ID5;
nTSS = length(mTSB);
print(nTSS);

windowSize = 100;
pReads = TSB[ width(TSB) >= 18 & width(TSB) <= windowSize & TSB$ID5 %in% TSSid ];
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
tmat[ cbind(as.character(pReads$ID5), width(pReads)) ] = pReads$C + pReads$U + pReads$R;



Stable = mTSB$Stability == "Stable";
Unstable = mTSB$Stability == "Unstable";
Early = mTSB$PsClass == "Early";
Late = mTSB$PsClass == "Late";
Both = mTSB$PsClass == "Both";

treatments = list(Capped=cmat, RppH=rmat, Uncapped=umat, Total=tmat);
subsets    = list(All=T, Stable=Stable, Unstable=Unstable, Early=Early, Late=Late, Mixed=Both);

for( s in names(subsets) ) {
	cond = subsets[[s]];
	pdf(paste0("../plots/maxTSB_", s, "_JointProbMatrices.pdf"), width=5, height=5);
	finished = c();

	for( t1 in names(treatments) ) {
		mat1 = treatments[[t1]][cond,];

		for( t2 in names(treatments) ) {
			if( paste0(t1, t2) %in% finished ) next;
			if( paste0(t2, t1) %in% finished ) next;
			finished = append( finished, paste0(t1, t2) );
			mat2 = treatments[[t2]][cond,];

			JPM = correlation_matrix( mat1, mat2 );
			JPM[cbind(1:100,1:100)] = 0;

			print(levelplot( x=JPM, xlim=c(18,75), ylim=c(18,75), at=seq(0, max(JPM), length.out=500),
				main=s, xlab=paste(t1, "RNA Length (nt)"), ylab=paste(t2, "RNA Length (nt)"),
				colorkey=list(col=plasma), col.regions=plasma, useRaster=T, aspect="fill"
			));
		}
	}

	dev.off();
}
