library(QDNAseq)
library(QDNAseq.hg38)
library(ggplot2)


cutoffDEL <- 0.5 + 0.5 - 0.5
cutoffLOSS <- 1.5 + 0.5 - 0.5
cutoffGAIN <-2.5 + 0.5 - 0.5
#copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff', cutoffs=log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10)/2))



bins <- getBinAnnotations(binSize=1000, genome="hg38")

readCounts <- binReadCounts(bins, bamfiles = ("W:/Clinical Working/SVI Project/Oxford Nanopore Technologies/Nanopore data from Maple/GM23366/megaruptor_15/GM23366_megaruptor15.sorted.bam"),minMapq =20 )
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
#readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, chromosomes =  "MT")
#readCountsFiltered
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)  
isobarPlot(readCountsFiltered)
readCountsFiltered <- estimateCorrection(readCountsFiltered)  
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes=NA)
isobarPlot(readCountsFiltered)
#readCountsFiltered <- estimateCorrection(readCountsFiltered)
noisePlot(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff', cutoffs=log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10)/2))
plot(copyNumbersCalled)


exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_1000bin_filter_cutoff.igv", format="igv",filter = TRUE)
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_1000bin_calls_cutoff.igv", format="igv",type = "calls")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_1000bin_segments_cutoff.igv", format="igv",type = "segments")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_1000bin_copynumber_cutoff.igv", format="igv",type = "copynumber")

#Bin500
bins <- getBinAnnotations(binSize=1000, genome="hg38")

readCounts <- binReadCounts(bins, bamfiles = ("W:/Clinical Working/SVI Project/Oxford Nanopore Technologies/Nanopore data from Maple/GM23366/megaruptor_15/GM23366_megaruptor15.sorted.bam"),minMapq =20 )
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
#readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, chromosomes =  "MT")
#readCountsFiltered
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)  
isobarPlot(readCountsFiltered)
readCountsFiltered <- estimateCorrection(readCountsFiltered)  
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes=NA)
isobarPlot(readCountsFiltered)
#readCountsFiltered <- estimateCorrection(readCountsFiltered)
noisePlot(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff', cutoffs=log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10)/2))
plot(copyNumbersCalled)


exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_500bin_filter_cutoff.igv", format="igv",filter = TRUE)
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_500bin_calls_cutoff.igv", format="igv",type = "calls")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_500bin_segments_cutoff.igv", format="igv",type = "segments")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_500bin_copynumber_cutoff.igv", format="igv",type = "copynumber")

#Bin100

bins <- getBinAnnotations(binSize=100, genome="hg38")

readCounts <- binReadCounts(bins, bamfiles = ("W:/Clinical Working/SVI Project/Oxford Nanopore Technologies/Nanopore data from Maple/GM23366/megaruptor_15/GM23366_megaruptor15.sorted.bam"),minMapq =20 )
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
#readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, chromosomes =  "MT")
#readCountsFiltered
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)  
isobarPlot(readCountsFiltered)
readCountsFiltered <- estimateCorrection(readCountsFiltered)  
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes=NA)
isobarPlot(readCountsFiltered)
#readCountsFiltered <- estimateCorrection(readCountsFiltered)
noisePlot(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff', cutoffs=log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10)/2))
plot(copyNumbersCalled)


exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_100bin_filter_cutoff.igv", format="igv",filter = TRUE)
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_100bin_calls_cutoff.igv", format="igv",type = "calls")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_100bin_segments_cutoff.igv", format="igv",type = "segments")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_100bin_copynumber_cutoff.igv", format="igv",type = "copynumber")

#bin50

bins <- getBinAnnotations(binSize=50, genome="hg38")

readCounts <- binReadCounts(bins, bamfiles = ("W:/Clinical Working/SVI Project/Oxford Nanopore Technologies/Nanopore data from Maple/GM23366/megaruptor_15/GM23366_megaruptor15.sorted.bam"),minMapq =20 )
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
#readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, chromosomes =  "MT")
#readCountsFiltered
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)  
isobarPlot(readCountsFiltered)
readCountsFiltered <- estimateCorrection(readCountsFiltered)  
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes=NA)
isobarPlot(readCountsFiltered)
#readCountsFiltered <- estimateCorrection(readCountsFiltered)
noisePlot(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff', cutoffs=log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10)/2))
plot(copyNumbersCalled)


exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_50bin_filter_cutoff_noy.igv", format="igv",filter = TRUE)
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_50bin_calls_cutoff_noy.igv", format="igv",type = "calls")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_50bin_segments_cutoff_noy.igv", format="igv",type = "segments")
exportBins(copyNumbersCalled,file = "GM23366_QDNAseq_50bin_copynumber_cutoff_noy.igv", format="igv",type = "copynumber")

