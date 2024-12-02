#!/usr/bin/Rscript


hg19chr = paste0( 'chr', c(1:22, 'X', 'Y', 'M') );
hg19len   = c(  249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
                135006516, 133851895, 115169878, 107349540, 102531392,  90354753,  81195210,  78077248,  59128983,  63025520,
                 48129895,  51304566, 155270560,  59373566,  16571
             );

# Shorthand for "pretty number" formatting
pn = function(value) {
    prettyNum(value, big.mark=",")
}

# Shorthand to print the full argument list
msgout = function(...) {
    write(paste(...), stdout());
}

# Compute universal coordinate IDs for 5' and 3' ends
compute_IDs = function(reads) {
    chrOffsets = cumsum(as.numeric(c(0, hg19len)));
    chrOffsets = chrOffsets[1:length(chrOffsets)-1];
    names(chrOffsets) = hg19chr;
    offset = chrOffsets[as.character(seqnames(reads))];
    names(offset) = NULL;
    strdir = ifelse(strand(reads) == '+', 1, -1);
    RNA5 = ifelse( strdir>0, start(reads),   end(reads) );
    RNA3 = ifelse( strdir>0,   end(reads), start(reads) );
    reads$ID5 = strdir * (offset + RNA5);
    reads$ID3 = strdir * (offset + RNA3);
    return(reads);
}

# Function to export reads in bigWig format
CoPROtoBigWig = function( reads, path, wname ) {
    make_bigwig = function(reads, filename, wname) {
        pl = strand(reads) == '+';
        cov = coverage(reads[pl], weight=wname);
        export(cov, paste0(filename, '_pl.bw'));

        mcols(reads)[,wname] = -mcols(reads)[,wname];
        cov = coverage(reads[!pl], weight=wname);
        export(cov, paste0(filename, '_mn.bw'));
    }

    # Export 5' ends
    out = resize(reads, width=1, fix='start');
    make_bigwig( out, paste0(path, '_5p'), wname );

    # Export 3' ends
    out = resize(reads, width=1, fix='end');
    make_bigwig( out, paste0(path, '_3p'), wname );
}

# Find local maxima among the given peaks within
# a given distance. Strand-specific.
# pos = peak positions
# h   = peak heights
# d   = maximum distance allowed between peaks.
# returns: indices of local maxima peaks.
find_local_max = function( pos, h, d ) {
    if( length(pos) != length(h) ) {
        print("len(peaks) != len(heights)");
        return(NULL);
    }
    N = 1:length(pos);
    isMin = T;

    while( any(isMin) & length(N) > 1 ) {
        O = N[1:length(N)-1];
        N = N[2:length(N)  ];
        # find peaks within the same window
        isClust = pos[N] - pos[O] <= d;

        # remove minima from each cluster
        NMin = c( F, (h[N] <  h[O]) & isClust    );
        OMin = c(    (h[O] <= h[N]) & isClust, F );
        isMin = NMin | OMin;
        N = c(O[1], N);
        N = N[ !isMin ];
    }
    return(N);
}


# Cluster TSBs within a given distance. Strand-specific.
# tsb = Transcription Start Bases
#   d = distance to cluster
# tsb must be sorted by ID5 prior to using this function!
compute_TSS = function( tsb, d ) {
    # find left- or right-most TSB within each cluster
    # left = find_local_max( tsb$ID5, -abs(tsb$ID5), d );
    #right = find_local_max( tsb$ID5,  abs(tsb$ID5), d );
    clusters = resize(tsb, width=2*d+1, fix="center");
    clusters = reduce(clusters);
    start(clusters) = start(clusters) + d;
    end(clusters) = end(clusters) - d;
    
    #cstr = strand(tsb)[left];
    #cname = seqnames(tsb)[left];
    #clusters = paste0(cname, ':', start(tsb)[left], '-', start(tsb)[right]+1, ':', cstr);
    #head(clusters);
    #clusters = as( clusters, "GRanges" );
    #clusters$ID5 = tsb[left]$ID5;

    return(clusters);
}

# Cluster TSBs within a given distance. Strand is ignored.
# tsb = Transcription Start Bases
#   d = distance to cluster
compute_TREs = function( tsb, d ) {
    # reorder on unstranded position
    temp = tsb[ order(abs(tsb$ID5)) ];
    temp$ID5 = abs(temp$ID5);
    temp = temp[ !duplicated(temp$ID5) ];
    # find left- or right-most TSB within each cluster
     left = find_local_max( temp$ID5, -temp$ID5, d );
    right = find_local_max( temp$ID5,  temp$ID5, d );

    # create new GRanges object spanning left to right
    cname = seqnames(temp)[left];
    clusters = as( paste(cname, paste0(start(temp)[left], '-', start(temp)[right]), sep=':'), "GRanges");
    clusters$ID5 = tsb[left]$ID5;

    return(clusters);
}


# compute chromosome strand and position from a given ID
findChr = function( pos ) {
    chrOffsets = cumsum(as.numeric(c(0, hg19len)));
    chrOffsets = chrOffsets[1:length(chrOffsets)-1];
    names(chrOffsets) = hg19chr;

    m = min(which( chrOffsets > abs(pos) ))-1;
    ptitle = paste( names(chrOffsets)[m], format(abs(pos)-chrOffsets[m], scientific=F) )

    return(ptitle);
}


# draw a 2D heatmap of read counts
pairplot = function( reads, target, windowSize ) {
    sameStrand = abs(reads$ID5 - target) <= windowSize;
    diffStrand = abs(reads$ID5 + target) <= windowSize;
    sr = reads[ (sameStrand | diffStrand) & width(reads) <= windowSize ];
    sameStrand = abs(sr$ID5 - target) <= windowSize;
    sr$pos5 = abs(sr$ID5) - abs(target);
    sr$pos3 = abs(sr$ID3) - abs(target);
    sr = sr[ abs(sr$pos5) <= windowSize & abs(sr$pos3) <= windowSize, ];

    ptitle = findChr( target );
    #dev.new(width=8, height=8);
    ggplot( as.data.frame(sr), aes(x=pos5, y=pos3) ) + ggtitle(ptitle) +
        stat_bin2d( aes(fill=C), binwidth=c(1,1)  ) +
        scale_x_continuous(name = "5' Position (nt)", lim=c(-1*windowSize, windowSize)) +
        scale_y_continuous(name = "3' Position (nt)", lim=c(-1*windowSize, windowSize)) +
        scale_fill_gradientn(name = "Reads", colours=c('azure2', 'steelblue3', 'red3') ) +
        theme( panel.background = element_blank(), legend.position=c(0.85, 0.2),
            axis.line = element_line(colour = "black") );
}

# draw a plot of capped and uncapped pause profiles at a given site
# pooled = pooled reads containing R, C, U counts
# target = ID5 to be plotted
# windowSize = max distance to be plotted
pauseplot = function( pooled, target, windowSize = 100 ) {
    hits = pooled[ pooled$ID5 == target & width(pooled) <= windowSize & pooled$MapQ >= 30 ];
    hits = hits[ order(width(hits)) ];
    wtA = 0.17;
    wtB = 2.6;

    pcols = c('steelblue3', 'red3', 'green3');
    ptitle = findChr( target );
    sumC   = sum(hits$C);
    sumU   = sum(hits$U);
    dev.new(width=6, height=4);
    par(mar=c(5.1, 4.1, 3.1, 2.1));
    matplot( x=width(hits), y=hits[, c('C', 'U', 'R')], col=pcols, main=ptitle,
        type='l', lty=1, lwd=1.5, xlim=c(18, windowSize), xlab='Distance from TSS (nt)', ylab='Normalized Reads'
    );
    legend('topright', c('Capped', paste('Uncapped', sumU/wtB), 'RppH'), lty=1, col=pcols, lwd=1.5, bty='n' );
}


# compute a pause matrix from the given reads
pause_matrix = function( reads, fname, windowSize, as.prob=F ) {
    reads = reads[ width(reads) <= windowSize ];
    TSSid = unique(reads$ID5);
    nTSS = length(TSSid);
    windowSize = 100;
    mat3p = matrix(data=0, nrow=nTSS, ncol=windowSize);
    rownames(mat3p) = TSSid;
    colnames(mat3p) = as.character(1:windowSize);
    mat3p[ cbind(as.character(reads$ID5), width(reads)) ] = mcols(reads)[,fname];

    if( as.prob ) mat3p = mat3p / rowSums(mat3p);
    return(mat3p);
}


# compute the joint probability matrix between two matrices
# correlations are computed row-by-row
# A and B must have equal dimensions
# rows of A and B should be probabilities that sum to 1
correlation_matrix = function( A, B ) {
    library(foreach);
    library(parallel);
    library(iterators)

    if( any(dim(A) != dim(B)) ) {
        print("correlation_matrix: the supplied matrices have different dimensions");
        return(NA);
    }
    getProbs = function(x) {
        N=length(x)/2;
        A = x[1:N];
        B = x[(N+1):(2*N)];
        p = A %*% t(B);
        norm = sum(p, na.rm=T);
        if(norm > 0) {
            return( p / norm );
        } else
            return( p );
    }

    joined = cbind(A, B);
    # iterate over the joined matrix by rows (00 at a time)
    # combine results by addition
    cormat = foreach( j=iter(joined, by='row', chunksize=500), .combine='+' ) %dopar% {
        rowSums(
            apply(j, 1, getProbs)
        );
    }
    cormat = matrix(cormat, nrow=ncol(A)) / sum(cormat);
    return( cormat );
}


# Input  : two GRanges (windows first, counts second)
# Returns: a matrix of counts within windows
# windowSize = size of all 'windows' in the matrix
# fix = how to position windows. [start, center, end]
OverlapMatrix = function( windows, counts, windowSize, field='score', fix='center' ) {
    fixed = resize( windows, width=windowSize, fix=fix );
    cmat = matrix(data=0, nrow=length(fixed), ncol=windowSize);

    hits = findOverlaps( fixed, counts );
    offset = start(counts)[hits@to] - start(fixed)[hits@from];
    for( x in 1:max(width(counts)[hits@to]) ) {
        valid = width(counts)[hits@to] >= x & offset > -x & offset+x < windowSize;
        cmat[ cbind( hits@from, offset+x )[valid,] ] = mcols(counts)[hits@to[valid],field];
    }
    return(matrix(data=cmat, nrow=length(fixed), ncol=windowSize));
}


MetaHeatmap = function( CountMatrix, nRows, nPermut, SortVector, FUN='mean' ) {
    CountMatrix = CountMatrix[SortVector,];
    MetaVecs = matrix(0, nrow=nRows*ncol(CountMatrix), ncol=nPermut);
    for(i in 1:nPermut) {
        Indices = sample(nrow(CountMatrix), nRows, replace = F);
        Indices = sort(Indices);
        MetaVecs[,i] = as.vector( CountMatrix[Indices,] );
    }
    MetaVecs = apply(MetaVecs, 1, FUN);
    MetaMatrix = matrix( MetaVecs, nrow=nRows, ncol=ncol(CountMatrix) );
    return( MetaMatrix );
}


LinearHeatmap = function( CountMatrix, nRows, SortVector, FUN='mean', ... ) {
    CountMatrix = CountMatrix[SortVector,];
    cksz = floor(nrow(CountMatrix)/nRows);
    myit = iter(CountMatrix, by='row', chunksize=cksz);
    linMat = foreach( chunk = myit, .combine='rbind' ) %dopar% {
        apply(chunk, 2, FUN, ...);
    }
    return( linMat[1:nRows,] );
}




heatcols = c('white', "#A7A7A7", "dodgerblue", "forestgreen","gold", "firebrick");
tblue = adjustcolor('deepskyblue4', alpha.f=0.1);
tred  = adjustcolor('firebrick', alpha.f=0.1);
tgreen= adjustcolor('green3', alpha.f=0.1);
tblack= adjustcolor('black', alpha.f=0.2);
