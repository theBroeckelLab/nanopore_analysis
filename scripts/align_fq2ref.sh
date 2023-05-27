
#################################
## EXPORT FILE NAME AND REF #####
#################################
export PREFIX="mcw_svi_0042"
export REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
export FASTQ="/data/svi/prom/mcw_svi/flowcell_10.4.1/mcw_svi_0042/raw_data/fastq_pass/mcw_svi_0042.fastq.gz"

###################################
######################
### READ ALIGNMENT ###
######################
minimap2 -L -t 96 -ax map-ont $REF $FASTQ | samtools view -bS -@ 96 > sample.bam 
samtools sort -@ 96 sample.bam -o ${PREFIX}.sorted.bam 
samtools index -@ 16 ${PREFIX}.sorted.bam 
rm sample.bam


###########################
### QC - POST ALIGNMENT ###
###########################
#pycoQC -f sequencing_summary_PAK63063_a1d9c21d_a1480da3.txt -a ${PREFIX}.sorted.bam -o pycoQC_alignment.html   
  

