#!/usr/bin/Rscript

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicAlignments)
	library(GenomicFiles)
	library(GenomicRanges)
	library(rtracklayer)
	library(RColorBrewer)
	library(ggplot2)
	library(RWebLogo)
});
genomef = "~/bin/genomes/mm10/mm10.2bit";


load('../data/mm10/MEF_AllCalls.Rdata');


PsBs = TSB[ which(TSB$ID5 %in% NHS.TSS$mTSB & TSB$NC >= 10) ];
PsBs = PsBs[order(-PsBs$NC)];
PsBs = PsBs[!duplicated(PsBs$ID5)];
NHSd = density(width(PsBs));
plot(NHSd);

PsBs = GenomicRanges::resize(PsBs, width=1, fix='end');
flank3 = GenomicRanges::promoters(PsBs, upstream=6, downstream=5);
seq3 = import.2bit(genomef, which=flank3);
seq3[ which(strand(PsBs) == '-') ] = reverseComplement(seq3[ which(strand(PsBs) == '-') ]);
writeXStringSet(seq3, file="../out/MEF_maxPause.fa");
weblogo(file.in="../out/MEF_maxPause.fa", first.index=-5, title="MEF maxPause motif");

PsBs = TSB[ which(TSB$ID5 %in% HS.TSS$mTSB & TSB$HC >= 10) ];
PsBs = PsBs[order(-PsBs$HC)];
PsBs = PsBs[!duplicated(PsBs$ID5)];
HSd = density(width(PsBs));
plot(HSd);

PsBs = GenomicRanges::resize(PsBs, width=1, fix='end');
flank3 = GenomicRanges::promoters(PsBs, upstream=6, downstream=5);
seq3 = import.2bit(genomef, which=flank3);
seq3[ which(strand(PsBs) == '-') ] = reverseComplement(seq3[ which(strand(PsBs) == '-') ]);
writeXStringSet(seq3, file="../out/MEF_maxPause.fa");
weblogo(
	file.in="../out/MEF_maxPause.fa",
	file.out="../out/MEF_maxPause.pdf",
	first.index=-5,
	title="MEF Max Pause",
	verbose=F
);

PsBs = TSB[ which(TSB$NC >= 3 & TSB$NC <= 10) ];
NHSd = density(width(PsBs));
plot(NHSd);

PsBs = GenomicRanges::resize(PsBs[1:10^5], width=1, fix='end');
flank3 = GenomicRanges::promoters(PsBs, upstream=6, downstream=5);
seq3 = import.2bit(genomef, which=flank3);
seq3[ which(strand(PsBs) == '-') ] = reverseComplement(seq3[ which(strand(PsBs) == '-') ]);
writeXStringSet(seq3, file="../out/MEF_allPause.fa");
weblogo(
	file.in="../out/MEF_allPause.fa",
	file.out="../out/MEF_allPause.pdf",
	first.index=-5,
	title="MEF All Pause",
	verbose=F
);


PsBs = TSB[ TSB$ID5 %in% NHS.TSS$mTSB & TSB$NC >= 5 ];
PsBs = PsBs[order(-PsBs$NC)];
PsBs = PsBs[!duplicated(PsBs$ID5) & width(PsBs) <= 50];
psRNA = resize(PsBs, width=1, fix='start');
psRNA = GenomicRanges::promoters(psRNA, upstream=0, downstream=32);
seq3 = import.2bit(genomef, which=psRNA);
seq3[ which(strand(psRNA) == '-') ] = reverseComplement(seq3[ which(strand(psRNA) == '-') ]);
names(seq3) = width(PsBs);
writeXStringSet(seq3, file="../seqlearn/MEF_maxPause.fa.gz", compress=T);


PsBs = TSB[ TSB$ID5 %in% HS.TSS$mTSB & TSB$HC >= 20 ];
PsBs = PsBs[order(-PsBs$HC)];
PsBs = PsBs[!duplicated(PsBs$ID5)];
PsBs = resize(PsBs, width=1, fix='start');
psRNA = GenomicRanges::promoters(PsBs, upstream=0, downstream=32);
seq3 = import.2bit(genomef, which=psRNA);
seq3[ which(strand(psRNA) == '-') ] = reverseComplement(seq3[ which(strand(psRNA) == '-') ]);
names(seq3) = 1:length(seq3);
writeXStringSet(seq3, file="../seqlearn/MEF_goodPause.fa.gz", compress=T);
ngood = length(seq3);


PsBs = TSB[ !TSB$ID5 %in% HS.TSS$mTSB & !TSB$ID5 %in% NHS.TSS$mTSB & TSB$NC <= 15 & TSB$HC <= 15 ];
PsBs = PsBs[!duplicated(PsBs$ID5)];
PsBs = resize(PsBs[1:ngood], width=1, fix='start');
psRNA = GenomicRanges::promoters(PsBs, upstream=0, downstream=32);
seq3 = import.2bit(genomef, which=psRNA);
seq3[ which(strand(psRNA) == '-') ] = reverseComplement(seq3[ which(strand(psRNA) == '-') ]);
names(seq3) = 1:length(seq3);
writeXStringSet(seq3, file="../seqlearn/MEF_badPause.fa.gz", compress=T);

