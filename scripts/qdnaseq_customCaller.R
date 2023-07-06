
## QDNASEq Custom Script 07/06/2023

#############################
## Import BAM file ##########
#############################
## Create New
#bamfile="W:/Clinical Working/SVI Project/Oxford Nanopore Technologies/Nanopore data from Maple/GM23366/megaruptor_15/GM23366_megaruptor15.sorted.bam"
#bins=getBinAnnotations(binSize=50, genome="hg38")
#readCounts=binReadCounts(bins, bamfiles=(bamfile), minMapq=20)
#saveRDS(readCounts, file="Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2023-02-01 scripts for ONT/QDNASeq_test_object_gm23366_bin50k.RData")
## OR Read in one that was already made
readCounts=readRDS(file="Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2023-02-01 scripts for ONT/QDNASeq_test_object_gm23366_bin500k.RData")


############################
##### run base QDNAseq #####
############################
readCountsFiltered <- applyFilters(readCounts, chromosomes=c("X","Y"))  
readCountsFiltered <- estimateCorrection(readCountsFiltered)
#readCountsFiltered <- estimateCorrection(readCountsFiltered, control=loess.control(surface="direct"))
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes=NA)
#Biobase::fData(readCountsFiltered)$use=Biobase::fData(readCountsFiltered)$use & !is.na(rowMeans(assayDataElement(readCountsFiltered, "fit")))
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)


############################################################
#### Pull Info from QDNAseq object into dataframe 'cn' #####
############################################################
bins=rownames(copyNumbersSegmented@assayData$copynumber)
cn=data.frame(bin=bins, chr=NA, start=NA, end=NA, cn=copyNumbersSegmented@assayData$copynumber[,1], cn.sexAdjusted=NA, cn.sexAdjusted.log2=NA, cn.call=NA)
rownames(cn)=1:nrow(cn)
cn$chr=as.character(sapply(strsplit(bins, ":"), `[`, 1))
cn$start=as.numeric(sapply(strsplit(as.character(sapply(strsplit(bins, ":"), `[`, 2)), "-"), `[`,1))
cn$end=as.numeric(sapply(strsplit(as.character(sapply(strsplit(bins, ":"), `[`, 2)), "-"), `[`,2))


################################################
## Determine Sex; if Male multiply X/Y cn 2x ###
################################################
cn$cn.sexAdjusted=cn$cn
meanX=mean(cn$cn[which(cn$chr=="X")], na.rm=T)
meanY=mean(cn$cn[which(cn$chr=="Y")], na.rm=T)
if((meanX>0.25&meanX<0.75) & (meanY>0.25&meanY<0.75)) {sex="Male"}
if((meanX>0.75&meanX<1.25) & (meanY<0.25)) {sex="Female"}
cat(paste0("ChrX: ", round(meanX,2), "\nChrY: ", round(meanY,2), "\nSex Determination is ", sex))
if(sex=="Female") {cn$cn.sexAdjusted[which(cn$chr=="Y")]=NA}
if(sex=="Male") {cn$cn.sexAdjusted[which(cn$chr=="X"|cn$chr=="Y")]=(2*cn$cn[which(cn$chr=="X"|cn$chr=="Y")])}


#############################
## transform cn to log2 #####
#############################
cn$cn.sexAdjusted.log2=log2(cn$cn.sexAdjusted)


############################################
## Specify cutoffs and identify cn calls ###
#############################################
#specify cutOffs
cutoffLOSS=log2(1.5/2); cutoffGAIN=log2(2.5/2)
#loop through each bin and determine call
cn$cn.call=0
for (i in 1:nrow(cn)) {
  if(is.na(cn$cn.sexAdjusted.log2[i])) {next}
  if(cn$cn.sexAdjusted.log2[i]<=cutoffLOSS) {cn$cn.call[i]=(-1); next}
  if(cn$cn.sexAdjusted.log2[i]>=cutoffGAIN) {cn$cn.call[i]=(1); next}
}


###########################################################
#### Merge bins with adjacent calls of the same type #####
#########################################################
cn$callToWrite=NA; variant_count=0
for (i in 1:nrow(cn)) {
  if(cn$cn.call[i]==0) {next}
  if((abs(cn$cn.call[i])==1&abs(cn$cn.call[i-1])!=1)) {
    variant_count=variant_count+1
    cn$callToWrite[i]=paste0("CNV_",variant_count)
  }
  if((abs(cn$cn.call[i])==1&abs(cn$cn.call[i-1])==1) & cn$chr[i]==cn$chr[i-1]) {
    cn$callToWrite[i]=paste0("CNV_",variant_count)
  }
}
cn$cn.call[which(is.na(cn$cn))]=NA
cn$callToWrite=factor(cn$callToWrite, levels=unique(cn$callToWrite))



# ####################################################
# ###### merge calls seperated by nbases (NAs) #######
# ####################################################
# all.cnCalls=setdiff(unique(cn$callToWrite), NA)
# for (i in 1:length(all.cnCalls)) {
#   if(isFALSE(all.cnCalls[i]%in%cn$callToWrite)) {next}
#   idx.start=which(cn$callToWrite==all.cnCalls[i])[1]
#   newRange=c()
#   for (callExtend in idx.start:nrow(cn)) {
#     if(identical(cn$cn.call[callExtend], 0)) {break}
#     newRange=c(newRange, callExtend)
#   }
#   cn$callToWrite[newRange]=all.cnCalls[i]
# }
# cn$callToWrite=factor(cn$callToWrite, levels=unique(cn$callToWrite))




###################################################
## filter for merged calls and write to new df ####
###################################################
split.calls=split(cn, cn$callToWrite)
cn.merge=data.frame(chr=NA, start=NA, end=NA, type=NA, length=NA, meanCN=NA, bins=NA, id=names(split.calls))
for (i in 1:length(split.calls)) {
  cn.merge$chr[i]=split.calls[[i]]$chr[1]
  cn.merge$start[i]=min(split.calls[[i]]$start)
  cn.merge$end[i]=max(split.calls[[i]]$end)
  if(split.calls[[i]]$cn.call[1]==(-1)) {cn.merge$type[i]="DEL"}
  if(split.calls[[i]]$cn.call[1]==(1)) {cn.merge$type[i]="DUP"}
  cn.merge$length[i]=as.numeric(cn.merge$end[i])-as.numeric(cn.merge$start[i])
  cn.merge$meanCN[i]=mean(split.calls[[i]]$cn.sexAdjusted)
  cn.merge$bins[i]=nrow(split.calls[[i]])
  cn.merge$id[i]=as.character(split.calls[[i]]$callToWrite[1])
}
View(cn.merge)


##############################
## Plot a single chromosome ##
##############################
chrOI="21"
cn.ToPlot=cn[which(cn$chr==chrOI),]
cn.ToPlot$callToWrite_forPlot="NA"
cn.ToPlot$callToWrite_forPlot[which(!is.na(cn.ToPlot$callToWrite))]="CALL"
colors=c("NA"="black", "CALL"="red")
ggplot(cn.ToPlot, aes(x=start, y=cn.sexAdjusted, color=callToWrite_forPlot))+
  geom_point()+
  scale_color_manual(values=colors)+
  ylab("Copy Number")+xlab(paste0("chromosome ", chrOI, " position"))+
  theme_bw()+
  #ylim(-2,2)+
  theme(legend.position="none")



#############################
## Filter by cnv size #####
###########################
#minSize=100000
#cn.merge=cn.merge[which(cn.merge$length>minSize),]
View(cn.merge)







############################
#### compare to array ######
############################
cma=read.xlsx("W:/Clinical Working/SVI Project/Bionano/Bionano Working Group/Bionano output/Validation_comparison_output/", sheet=4)
View(cma)
