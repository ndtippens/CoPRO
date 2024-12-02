# CoPRO analysis scripts
Scripts to facilitate analysis of CoPRO datasets.

## Change Log:
v1.0: Initial commit of Preprocessing and Mappability files. There are slight differences from the analysis reported in Nature Genetics:
 - R2 barcode trimming was conditioned on the R1 barcode to minimize mismatched R1/R2 trims.
 - Insert size bias normalization uses recently reported measurements for the NextSeq (https://www.biorxiv.org/content/early/2018/08/09/388108.1).
 - These changes do not significantly affect any conclusions reported in the manuscript.

## Terminology
TSN = Transcription start nucleotide
TSS = Transcription start site
TID = Transcription initiation domain
maxTSN = Max start nucleotide within a TSS
maxTSS = Max start site within a TID
MMRL = Minimum mappable RNA length for a TSN

## Preprocessing
1. Trim barcodes & perform sequence alignment (CoPRO_align.sh)
2. Pool multiple sequencing lanes (samtools merge)
3. Read, filter, and adjust alignments (ConvertBams.r)
4. Normalize & create an integrated dataset (MergeTreatments.r)
5. Global analysis to identify TSNs, TSSes, and TIDs (CoPRO_IdentifyAll.ipynb)

## Mappability lengths
6. Extract all TSN sequences as 19-24mers (Mapping_ExtractSeqs.r)
7. Align extracted sequences (Mapping_AlignSeqs.sh)
8. Record MMRL for each TSN (Mapping_RecordMMRL.r)

## Initiation, capping & pausing analyses
- TSN_PauseSharing.ipynb: analysis of pause sites shared by TSNs
- maxTSN_Pausing.ipynb: pause classification and analysis at maxTSNs
- maxTSN_PauseReplicates.r: comparison of maxTSN pausing in each replicate
- maxTSN_JointProbAnalysis.r: comparison of pause distributions using joint probabilities
- maxTSN_Capping.ipynb: capping analysis at maxTSNs
- maxTSN_SeqHeatmaps.r: analysis of GC content vs pause distance at maxTSNs
- maxTSN_PauseClassMotifs.r: analysis of TF motifs around maxTSNs
- maxTSN_Spacing.r: analysis of spacing between all maxTSNs (convergent, divergent, same-strand)

## TID analyses
- TID_LinearHeatmaps.ipynb: Generate linear row-collapsed (and raw) heatmaps of ChIP datasets around TIDs
