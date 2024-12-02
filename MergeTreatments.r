#!/usr/bin/Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(data.table)
    library(rtracklayer)
    library(GenomicRanges)
    library(doParallel)
});
registerDoParallel(cores=7);
source("./CoPRO_Functions.r");

# use (monotonic) read weights from http://clinchem.aaccjnls.org/content/56/8/1279
len = 3:25 * 20;
nrd = c(2.35,2.62,2.3,2.31,1.3,1.48,1.65,1.2,.5,.8,.79,.81,.65,.54,.67,.32,.51,.5,.48,.47,.25,.55,.27);
quake = data.frame(len, nrd);
exp.m = lm(log(nrd) ~ len, data=quake);
quake$mod = exp( predict(exp.m, list(len=quake$len)) );

RdW = data.frame( len=15:1000, w=1 );
RdW$mod = exp( predict(exp.m, list(len=RdW$len)) );
RdW$w = RdW$mod[RdW$len==80] / RdW$mod;
RdW[RdW$len<80,'w'] = 1;
RdW[RdW$len>400,'w'] = RdW$w[RdW$len==400];

pdf("../plots/Quake_ReadWeights.pdf", height=4, width=4);
plot(x=quake$len, y=quake$nrd, ylim=c(0,3), xlim=c(0,400), xlab="Insert Length", ylab="Read Number (Millions)");
lines(quake$len, quake$mod, lwd=2);
plot(x=RdW$len, y=RdW$w, ylim=c(0,10), xlim=c(0,400), xlab="Insert Length", ylab="Weight", type='l', lwd=2);
dev.off();
save(RdW, file="../data/Quake_ReadWeights.Rdata");


# NEW: use updated size bias estimate for NextSeq, from Gohl et al, https://www.biorxiv.org/content/biorxiv/early/2018/08/09/388108.1
len = 1:10 * 150 - 128;
nrd = c(77.434, 11.790, 5.183, 1.701, 1.820, 1.009, 0.645, 0.270, 0.128, 0.020);
nextseq = data.frame(len, nrd);

RdW = data.frame( len=22:500, w=1 );
exp.m = lm(log(nrd) ~ len, data=nextseq[1:4,]);
RdW$mod = exp( predict(exp.m, list(len=RdW$len)) );
RdW$w = max(RdW$mod) / RdW$mod;

pdf("../plots/Gohl_NextSeq_ReadWeights.pdf", height=4, width=4);
plot(x=len, y=nrd, ylim=c(0.001,100), xlim=c(0,1400), xlab="Insert Length", ylab="Read Abundance", log="y");
plot(x=c(nextseq$len, RdW$len), y=c(max(nrd)/nrd, RdW$w), col=c(rep('black', length(nrd)), rep('gray', 501-22)), log='y')
plot(x=RdW$len, y=RdW$w, ylim=c(0,50), xlim=c(0,500), xlab="Insert Length", ylab="Weight", type='l', lwd=2);
dev.off();
save(RdW, file="../data/NextSeq_ReadWeights.Rdata");


for( repn in 1:2 ) {
    load(paste0('../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Uncapped_Rep', repn, '.Rdata'));
    uncap = aligned[[1]];
    load(paste0('../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Capped_Rep', repn, '.Rdata'));
    cap = aligned[[1]];
    load(paste0('../data/CoPRO_hg19_dm6/CoPRO_AllMerge_RppH_Rep', repn, '.Rdata'));
    rpph = aligned[[1]];
    rm(aligned);
    invisible(gc());

    # Pool all reads
    rpph$C = 0L;
    rpph$U = 0L;
    rpph$R = as.integer(rpph$count);
    #rpph$E = 0L;
    rpph$count = NULL;
    cap$C = as.integer(cap$count);
    cap$U = 0L;
    cap$R = 0L;
    #cap$E = 0L;
    cap$count = NULL;
    uncap$C = 0L;
    uncap$U = as.integer(uncap$count);
    uncap$R = 0L;
    #uncap$E = 0L;
    uncap$count = NULL;
    Pooled = rpph;

    dupes = findOverlaps(cap, Pooled, type='equal');
    Pooled[dupes@to]$C = cap[dupes@from]$C;
    Pooled = append(Pooled, cap[-dupes@from]);
    rm(rpph, cap);
    invisible(gc());

    dupes = findOverlaps(uncap, Pooled, type='equal');
    Pooled[dupes@to]$U = uncap[dupes@from]$U;
    uncap = uncap[-dupes@from];
    invisible(gc());
    Pooled = append(Pooled, uncap);
    rm(uncap);
    invisible(gc());

    Pooled = Pooled[ order(Pooled$ID5) ];
    invisible(gc());

    # For this, we need the efficiency of the data.table library
    ReadTable = as.data.table(mcols(Pooled));
    # index table by ID5
    setkey(ReadTable, ID5);
    
    ReadTable[, rnalen := abs(ID5-ID3)+1];
    # Normalize counts by read length & round
    ReadTable[, C := round(C * RdW[min(rnalen,400)+6,'w'], 1) ];
    ReadTable[, U := round(U * RdW[min(rnalen,400)+6,'w'], 1) ];
    ReadTable[, R := round(R * RdW[min(rnalen,400)+6,'w'], 1) ];

    rm(ReadTable);
    invisible(gc());

    save(Pooled, file=paste0("../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Rep", repn, ".Rdata"));
}

# now merge the two replicates (optional)
load("../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Rep1.Rdata");
Rep1=Pooled;
load("../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Rep2.Rdata");

# Pool reads from the two replicates
dupes = findOverlaps(Rep1, Pooled, type='equal');
Pooled[dupes@to]$C = Pooled[dupes@to]$C + Rep1[dupes@from]$C;
Pooled[dupes@to]$U = Pooled[dupes@to]$U + Rep1[dupes@from]$U;
Pooled[dupes@to]$R = Pooled[dupes@to]$R + Rep1[dupes@from]$R;

Rep1 = Rep1[-dupes@from];
invisible(gc());
Pooled = append(Pooled, Rep1);
rm(Rep1);
invisible(gc());

save(Pooled, file="../data/CoPRO_hg19_dm6/CoPRO_AllMerge_Pooled.Rdata");

Pooled$ID5 = NULL;
Pooled$ID3 = NULL;
# write individual bigwigs for each treatment
for( treat in c('RppH', 'Cap', 'Uncapped') ) {
    t = substr(treat, 1, 1);
    path = paste0("../data/CoPRO_hg19_dm6/CoPRO_AllMerge_", treat);
    CoPROtoBigWig( Pooled[ mcols(Pooled)[,t] > 0 ], path, t );
    mcols(Pooled)[,t] = NULL;
    invisible(gc());
}
