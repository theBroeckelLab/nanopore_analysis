## sniffles
#args=commandArgs(trailingOnly = T)

## file i/o
fpath="/data/svi/prom/training/flowcell_10.4.1/GM23366/megaruptor_speed15/read_SVs/sniffles2/"


## read in vcf
#vcf=read.csv(args[1], sep="\t", comment.char="#", header=F, stringsAsFactors=F)
vcf=read.csv(paste0(fpath, "output.vcf"), sep="\t", comment.char="#", header=F, stringsAsFactors=F)
colnames(vcf)=c("chr","start","id","n","type","n2","qual","info","gt","geno")
#vcf=head(vcf, 1000)
## split INFO into individual columns
allFields=paste0(vcf$info, collapse=";")
allFields=strsplit(allFields, ";")[[1]]
#allFields=unique(sapply(strsplit(allFields, "="), `[`, 1))
allFields=c("SVTYPE","SVLEN","END","SUPPORT","CHR2")
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
vcf.new$SVLEN=abs(as.numeric(vcf.new$SVLEN))

## write vcf
vcf.new_toWrite=vcf.new[,c(1,2,13,11,12,14,10,15,3)]
colnames(vcf.new_toWrite)=c("chr","start","end","type","size","support","geno","bnd2","id")
write.table(vcf.new_toWrite, paste0(fpath, "sniffles2.bed"), sep="\t", quote=F, col.names = T, row.names=F)

## write vcf
#newfname=paste0(strsplit(args[1], "[.]")[[1]][1], ".bed")
#write.table(vcf.new_toWrite, newfname, sep="\t", quote=F, col.names = T, row.names=F)

