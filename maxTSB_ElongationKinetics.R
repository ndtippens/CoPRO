#!/usr/bin/Rscript
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFiles))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
registerDoParallel(cores=8);

load("../data/newK562_CoPROcalls.Rdata");
load("../data/newK562_CoPRO_maxTSB.Rdata");

setwd("~/Rmain/CoPRO/CoPRO/");
source("./CoPRO_Functions.r");
genome = "~/bin/genomes/hg19/hg19.2bit"

# Length of DNA bubble
BubLen = 12;
# Length of RNA-DNA hybrid (9 in yeast and human pol2 crystal structures)
HybLen = 9; # +1 after NTP incorporation
# Position of RNA-DNA hybrid relative to DNA bubble
HStart = 2;

# Pairing and stacking energies for DNA bases (kcal/mol).
# Values are for 0.1 M NaCl, 37 C
DpairWt  = c(AT=0.64, CG=0.12);
DstackWt = c(
		  	AA=-1.49, AT=-1.72, AG=-1.44, AC=-2.19,
		  	TA=-0.57, TT=-1.49, TG=-0.93, TC=-1.81,
		  	GA=-1.81, GT=-2.19, GG=-1.82, GC=-2.55,
		  	CA=-0.93, CT=-1.44, CG=-1.29, CC=-1.82
			);

# Propagation energies for RNA hybridization (kcal/mol).
# Dinucleotides describe the non-template DNA strand.
# Currently values are for 1 M NaCl, 37 C
RpropWt = c(
		  	AA=-1.0, AC=-2.1, AG=-1.8, AT=-0.9,
		  	CA=-0.9, CC=-2.1, CG=-1.7, CT=-0.9,
		  	GA=-1.3, GC=-2.7, GG=-2.9, GT=-1.1,
		  	TA=-0.6, TC=-1.5, TG=-1.6, TT=-0.2
			);

#TUNABLE
#kmax = 24.7;   # for T7 RNAP. NTP & PPi catalysis per second.
kmax = 8.33;    # Jonkers early class, mid = 8.33 nt/sec. Human RNA Pol2
Kdc  = 15.6E-6; # Kd of complementary base
Kdnc = 2E-2;    # Kd of noncomplementary base
kbT  = 9.83E-22;# cal

K0   = 1E9;     # rate pre-factor
Ebt  = 46.2*kbT;# energy barrier, back
E12  = 40  *kbT;# energy barrier from 1 -> 2
Sfwd = 3.1 *kbT;# barrier slope, fwd

# ribonucleotide concentrations [citation?]
rConc = c(  A=3152E-6,
            G= 468E-6,
            U= 567E-6,
            T= 567E-6, # (T=UTP for simplicity)
            C=  29E-6
);


# Gibb's Free Energy of dsDNA
dsDNA_GFE = function( seqs, w ) {
	N=width(seqs)[1];
	M=w-1;

	deltaG = foreach( s=iter(seqs), .combine='cbind', .maxcombine=10^3 ) %dopar% {
		# get AT and GC counts
		seqCt = letterFrequencyInSlidingView( s, view.width=w, letters=names(DpairWt) );
		# apply energy weights & sum
		pairing = seqCt %*% DpairWt;

		# get all digrams
		digrams = as.character(Views( s, start=1:(N-1), width=2 ));
		# compute scores in each sliding window
		stacking = foreach( i=0:(N-w), .combine='c', .maxcombine=10^3 ) %do% {
			# sum digram weights within this window
			swindow = digrams[ (1:M)+i ];
			sum( DstackWt[ swindow ] );
		}

		# convert to cals/mol
		cpm = 1000*(pairing + stacking);
		return(cpm);
	}

	return(deltaG);
};

#RNA-DNA Gibb's Free Energy
RNADNA_GFE = function( seqs, w ) {
	N=width(seqs)[1];
	M=w-1;

	deltaG = foreach( s=iter(seqs), .combine='cbind', .maxcombine=10^3 ) %dopar% {
		# get all digrams
		digrams = as.character(Views( s, start=1:(N-1), width=2 ));
		# compute scores in each sliding window
		propagation = foreach( i=0:(N-w), .combine='c', .maxcombine=10^3 ) %do% {
			# sum digram weights within this window
			swindow = digrams[ (1:M)+i ];
			sum( RpropWt[ swindow ] );
		}

		return(propagation);
	}

	# convert to cals/mol
	return(1000*deltaG);
};

# compute transcriptional elongation complex Gibbs Free Energy
# see doi:10.1016/j.jmb.2004.08.107 for details
# M = offset from active site (typically between -2 and +2)
TEC_GFE = function( M, seqs ) {
	bubG = -dsDNA_GFE( seqs, BubLen	);
	N = 1:nrow(bubG);
	U = max(M, 0);
	H = HybLen-M;

	hybG = RNADNA_GFE( seqs, H )[N+HStart,]; # Hybrid starts 2 bases after DNA bubble
	unpG = 0;
	if(U)
		unpG = RNADNA_GFE( seqs, U )[N+HStart+H,];

	return(bubG + hybG + unpG/2);
}


# compute "mainRate" term??
mainRate = function( seqs, ddG ) {
	cNTP = rConc[ simplify2array(strsplit(as.character(seqs), '')) ];
	num = kmax * cNTP;
	denom = Kdc*(1 + exp(ddG/kbT)) + cNTP;

	return(num/denom);
}

# require at least 100 reads (capped read weight is 0.22)
pTSB = mTSB[ mTSB$PsOcc >= 22 ];

# require mappability at 20 nt
load("../data/Mappable21mers.Rdata");
pTSB = pTSB[ pTSB$ID5 %in% MapID5 ];

# plot relative to the initiation base (and sorted by distance to the max pause)
pTSB = pTSB[ order(-pTSB$Peak) ];

TSBprom = promoters( pTSB, upstream=20+HStart+HybLen, downstream=81+BubLen );

seqs = import.2bit(genome, which=TSBprom);
seqs[ strand(pTSB) == '-' ] = reverseComplement( seqs[ strand(pTSB) == '-' ] );


#dGNn2 = TEC_GFE(-2, seqs );
#dGNn1 = TEC_GFE(-1, seqs );
dGN0  = TEC_GFE( 0, seqs );
dGN1  = TEC_GFE( 1, seqs );
dGN2  = TEC_GFE( 2, seqs );

Kmain = mainRate( subseq(seqs, start=HStart+HybLen, width=nrow(dGN0)), dGN1-dGN0 );
Kback = K0 * exp(          (dGN0 - Ebt)/kbT );
Kfwd1 = K0 * exp(          (dGN0 - E12)/kbT );
Kfwd2 = K0 * exp( (dGN1 - E12 -   Sfwd)/kbT );

# label pause bases
#idx = cbind( 80+TSBprom$Peak, 1:ncol(gmat) );
#gmat[idx] = 130;

xats  = 51 + c(-50, 0, 32, 80);
xlabs = c(-50, 'maxTSB', 32, 80);

pdf(file=paste0("../plots/maxTSB_ElongationKinetics.pdf"), width=4, height=6);

out = LinearHeatmap( t(dGN1-dGN0)/1000, 200, T, na.rm=T );
levelplot( x=t(out),
	main="dG of Incorporation", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

out = LinearHeatmap( t(dGN2-dGN1)/1000, 200, T, na.rm=T );
levelplot( x=t(out),
	main="dG of Fwd Tracking", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

out = LinearHeatmap( t(Kmain), 200, T, na.rm=T );
levelplot( x=t(out),
	main="Kmain", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

dev.off();


out = LinearHeatmap( t(Kback), 200, T, FUN='mean' );
levelplot( x=t(out),
	main="Kback", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

out = LinearHeatmap( t(Kfwd1), 200, order(pTSB$Stability), FUN='mean' );
levelplot( x=t(out),
	main="Kfwd1", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

out = LinearHeatmap( t(Kfwd2), 200, order(pTSB$Stability), FUN='mean' );
levelplot( x=t(out),
	main="Kfwd2", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
);

out = LinearHeatmap( t(Kmain-Kfwd1-Kback), 200, order(pTSB$Stability), FUN='mean' );
print(levelplot( x=t(out),
	main="Total", xlab="Position of Active Site (nt)", ylab="",
	colorkey=list(col=plasma), col.regions=plasma,
	useRaster=T, aspect="fill",
	scales=list(x=list(at=xats, tck=c(1,0), lab=xlabs), y=list(draw=F))
));
