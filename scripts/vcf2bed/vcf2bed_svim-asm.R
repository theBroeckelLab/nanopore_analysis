## SVIM-ASM
#args=commandArgs(trailingOnly=T)

## file i/o
fpath="/data/svi/prom/training/flowcell_10.4.1/GM23366/megaruptor_speed15/assembly_SVs/hapdup/hapdup/svim_asm/"


## read in vcf
#vcf=read.csv(args[1], sep="\t", comment.char="#", header=F)
vcf=read.csv(paste0(fpath, "variants.vcf"), sep="\t", comment.char="#", header=F, stringsAsFactors = F)
colnames(vcf)=c("chr","start","id","n","type","n2","qual","info","gt","geno")
#vcf=head(vcf, 1000)
## split INFO into individual columns
allFields=paste0(vcf$info, collapse=";")
allFields=strsplit(allFields, ";")[[1]]
allFields=unique(sapply(strsplit(allFields, "="), `[`, 1))
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

## fix DUP encoding
vcf.new$SVTYPE[which(vcf.new$SVTYPE=="DUP:TANDEM")]="DUP"

## absolute value of lengths
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
vcf.new_toWrite=vcf.new[,c(1,2,12,11,13,8,10,14,3)]
colnames(vcf.new_toWrite)=c("chr","start","end","type","size","support","geno","bnd2","id")
write.table(vcf.new_toWrite, paste0(fpath, "svim-asm.bed"), sep="\t", quote=F, col.names = T, row.names=F)

## write vcf
#newfname=paste0(strsplit(args[1], "[.]")[[1]][1], ".bed")
#write.table(vcf.new_toWrite, newfname, sep="\t", quote=F, col.names = F, row.names=F)