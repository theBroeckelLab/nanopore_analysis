## SMAP to BED
#"chr","start","end","type","length", "support","geno","bnd2","id"

## file i/o
fpath="/data/svi/prom/training/flowcell_10.4.1/GM23366/bionano_results/"

#Read in SMAP
smap <-read.csv(paste0(fpath, "variants_combine_filters_inMoleRefine1.smap"), skip = 6, sep="\t", stringsAsFactors = F)
colnames(smap)[1]="ID"
headers<-colnames(smap)
smap <- read.csv(paste0(fpath, "variants_combine_filters_inMoleRefine1.smap"), skip=7, col.names=headers, sep="\t", stringsAsFactors = F)

#Read in CNV file
cnv<-read.csv(paste0(fpath, "cnv_calls_exp_annotation_results.txt"), skip = 3, sep="\t", stringsAsFactors = F)
colnames(cnv)[1]="Id"


#Manipulate SMAP file
smap$ID=paste0("SMAP_", smap$ID)
smap$OverlapGenes[which(smap$OverlapGenes=="-")]=FALSE
smap$OverlapGenes[which(smap$OverlapGenes!=FALSE)]=TRUE
smap$bnd2=paste0("chr",smap$RefcontigID2, ":", round(smap$RefEndPos))
#colnames(smap)
#plot(smap$SVsize, smap$Size) ##Size vs SVsize
smap=smap[,c(1,3,7,8,9,10,30,25,32,33,44,46,42,48,49,53,28)]    ##SMAP format for 'newer' samples
#smap=smap[,c(1,3,7,8,9,10,30,25,32,33,38,40,36,42,43)]     ##SMAP format for training/coriell samples


#Manipulate CNV file
cnv$Id=paste0("CNV_", cnv$Id)
cnv$OverlapGenes[which(cnv$OverlapGenes=="-")]=FALSE
cnv$OverlapGenes[which(cnv$OverlapGenes!=FALSE)]=TRUE
cnv$Present_in_BNG_control_samples=NA
cnv$Present_in_BNG_control_samples_w_same_enzyme=NA
cnv$Fail_chimeric_score=NA
cnv$Self_molecules=NA
cnv$Self_molecules_count=NA
cnv$bnd2=NA
colnames(cnv)
cnv=cnv[,c(1:4,9,6,12,5,19,20,13,15,21:23,24,10)]

#Combine SMAP and CNVs
colnames(smap)=colnames(cnv)
comb=as.data.frame(rbind(smap, cnv))
comb=comb[order(comb$Chromosome, comb$Start),]

#Add column for masked variants
##**## masking only tracked for SVs in CNV file, not SMAP file ##**##
comb$Masked=FALSE
comb$Masked[grep("masked", comb$Type)]=TRUE

#Recode labels for SV types
comb$Type_general=comb$Type
comb$Type_general[c(grep("deletion", comb$Type), grep("loss", comb$Type))]="DEL"
comb$Type_general[c(grep("duplication", comb$Type), grep("gain", comb$Type))]="DUP"
comb$Type_general[grep("insertion", comb$Type)]="INS"
comb$Type_general[grep("inversion", comb$Type)]="INV"
comb$Type_general[grep("trans", comb$Type)]="BND"

##fix SV sizes
idx.inv=which(comb$Type_general=="INV")
comb$End[idx.inv]=comb$Start[idx.inv]+comb$Size[idx.inv]
comb$Size=round(comb$End-comb$Start)
comb$Size[which(comb$Type_general=="BND")]=NA
comb$bnd2[which(comb$Type_general!="BND")]=NA
colnames(comb)
bng_full_toWrite=comb[,c(2,3,4,19,8,15,17,16,1)]
colnames(bng_full_toWrite)=c("chr","start","end","type","size","support","geno","bnd2","id")
write.table(bng_full_toWrite, paste0(fpath, "bionano_full.bed"), quote=F, sep="\t", row.names = F, col.names = T)

#Apply filters
#general filters - not type specific
# idx.general=which((comb$Present_in_BNG_control_samples<=1 | 
#                      is.na(comb$Present_in_BNG_control_samples)) &
#                     (comb$Present_in_BNG_control_samples_w_same_enzyme<=1 | 
#                        is.na(comb$Present_in_BNG_control_samples_w_same_enzyme)) &
#                     (comb$Self_molecules_count>5 | is.na(comb$Self_molecules)) &
#                     (comb$OverlapGenes==TRUE | 
#                        comb$OverlapGenes==FALSE & comb$NearestNonOverlapGeneDistance<=3000) &
#                     (comb$Fail_chimeric_score!="fail" | is.na(comb$Fail_chimeric_score)) &
#                     (comb$Self_molecules!="no" | is.na(comb$Self_molecules)))
#idx.general

#type-specific filters - length and qual
# idx.spec=which((comb$Type_general=="<INS>" & comb$Size>1500 & comb$Confidence>0) |
#                  (comb$Type_general=="<INV>" & comb$Size>0 & comb$Confidence>0.70) |
#                  (comb$Type_general=="<DEL>" & comb$Algorithm=="assembly_comparison" & 
#                     comb$Size>1500 & comb$Confidence>0) |
#                  (comb$Type_general=="<DEL>" & comb$Algorithm=="Region-based" & 
#                     comb$Size>500000 & comb$Confidence>0.99) | 
#                  (comb$Type_general=="<DUP>" & comb$Algorithm=="assembly_comparison" & 
#                     comb$Size>0 & comb$Confidence==(-1)) |
#                  (comb$Type_general=="<DUP>" & comb$Algorithm=="Region-based" & 
#                     comb$Size>500000 & comb$Confidence>0.99) |
#                  (comb$Type_general=="<BND>" & comb$Confidence>=0.05))
#idx.spec


#**modified general filters**#
idx.general=which((comb$Present_in_BNG_control_samples<=1 | 
                     is.na(comb$Present_in_BNG_control_samples)) &
                    (comb$Present_in_BNG_control_samples_w_same_enzyme<=1 | 
                       is.na(comb$Present_in_BNG_control_samples_w_same_enzyme)) &
                    (comb$Self_molecules_count>5 | is.na(comb$Self_molecules)) &
                    (comb$Fail_chimeric_score!="fail" | is.na(comb$Fail_chimeric_score)) &
                    (comb$Self_molecules!="no" | is.na(comb$Self_molecules)))
idx.spec=which((comb$Type_general=="INS" & comb$Size>1500 & comb$Confidence>0) |
                 (comb$Type_general=="INV" & comb$Size>0 & comb$Confidence>0.70) |
                 (comb$Type_general=="DEL" & comb$Algorithm=="assembly_comparison" & 
                    comb$Size>1500 & comb$Confidence>0) |
                 (comb$Type_general=="DEL" & comb$Algorithm=="Region-based" & 
                    comb$Size>200000 & comb$Confidence>0.99) | 
                 (comb$Type_general=="DUP" & comb$Algorithm=="assembly_comparison" & 
                    comb$Size>0 & comb$Confidence==(-1)) |
                 (comb$Type_general=="DUP" & comb$Algorithm=="Region-based" & 
                    comb$Size>200000 & comb$Confidence>0.99) |
                 (comb$Type_general=="BND" & comb$Confidence>=0.05))



## Filter for final set of SVs
filtered.svs=comb[intersect(idx.general, idx.spec),]
#filtered.svs=filtered.svs[,c(2:4,17,1)]
#filtered.svs$Chromosome=paste0("chr", filtered.svs$Chromosome)
#filtered.svs$Start=round(filtered.svs$Start); filtered.svs$End=round(filtered.svs$End) 
bng_filtered_toWrite=filtered.svs[,c(2,3,4,19,8,15,17,16,1)]
colnames(bng_filtered_toWrite)=c("chr","start","end","type","size","support","geno","bnd2","id")
write.table(bng_filtered_toWrite, paste0(fpath, "bionano_highQC_controlBng01.bed"), quote=F, sep="\t", row.names = F, col.names = T)
#write.table(filtered.svs, "../Downloads/gm23366_svs/bionano_results/bionano_highQC.bed", quote=F, sep="\t", row.names = F, col.names = T)

