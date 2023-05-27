##**## currently too slow, need faster code ##**##

tic=Sys.time()


# ## svim-read
args=commandArgs(trailingOnly=T)

## read in vcf
#vcf=read.csv(args[1], sep="\t", comment.char="#", header=F, stringsAsFactors=F)
vcf=read.csv("Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2023-02-01 scripts for ONT/test_vcfs/svim.vcf", sep="\t", comment.char="#", header=F, stringsAsFactors=F)
colnames(vcf)=c("chr","start","id","n","type","n2","qual","info","gt","geno")
#vcf=head(vcf, 10000)
## split INFO into individual columns
allFields=paste0(vcf$info, collapse=";")
allFields=strsplit(allFields, ";")[[1]]
#allFields=unique(sapply(strsplit(allFields, "="), `[`, 1))
allFields=c("SVTYPE","END","SVLEN", "SUPPORT")
vcfInfo=strsplit(vcf$info, ";")
vcf.new=vcf
info=vcf.new$info
vcf.new$info=NA
for (j in 1:length(allFields)) {
  print(paste0("processing field ", j, " of ", length(allFields)))
  vcf.new[[allFields[j]]]=NA
  for (i in 1:nrow(vcf)) {
    field.n=paste0("\\b", allFields[j], "\\b")
    text.pull=vcfInfo[[i]][grep(field.n, vcfInfo[[i]])]
    if(length(text.pull)==0) {next}
    vcf.new[[allFields[j]]][i]=text.pull
    vcf.new[i, grep("=",vcf.new[i,])]=sapply(strsplit(as.character(vcf.new[i,grep("=", vcf.new[i,])]), "="),`[`, 2)
  }
}
## fix insertion end coords
for (i in 1:nrow(vcf.new)) {
  if (vcf.new$SVTYPE[i]=="INS") {vcf.new$END[i]=vcf.new$start[i]+as.numeric(vcf.new$SVLEN[i])}
}

## absolute value of sv length
vcf.new$SVLEN=abs(as.numeric(vcf.new$SVLEN))

## fix inversion size calculation
idx.inv=which(vcf.new$SVTYPE=="INV")
vcf.new$SVLEN[idx.inv]=as.numeric(vcf.new$END[idx.inv])-as.numeric(vcf.new$start[idx.inv])

## fix translocation breakends
vcf.new$bnd.pair=NA
vcf.new$bnd.pair[which(vcf.new$SVTYPE=="BND")]=vcf.new$type[which(vcf.new$SVTYPE=="BND")]
for (i in 1:nrow(vcf.new)) {
  if(vcf.new$SVTYPE[i]!="BND") {next}
  sub.idx=grep("[[]|[]]", unlist(strsplit(vcf.new$bnd.pair[i], "")))
  vcf.new$bnd.pair[i]=substr(vcf.new$bnd.pair[i], sub.idx[1]+1, sub.idx[2]-1)
}

## write vcf
vcf.new_toWrite=vcf.new[,c(1,2,12,11,13,14,10,15,3)]
colnames(vcf.new_toWrite)=c("chr","start","end","type","size","support","geno","bnd2","id")
write.table(vcf.new_toWrite, "../Downloads/gm23366_svs/svim.bed", sep="\t", quote=F, col.names = T, row.names=F)

## write vcf
#newfname=paste0(strsplit(args[1], "[.]")[[1]][1], ".bed")
#write.table(vcf.new_toWrite, newfname, sep="\t", quote=F, col.names = F, row.names=F)



Sys.time()-tic