
export ASM="/data/svi/prom/training/flowcell_10.4.1/GM23366/megaruptor_speed15/assembly_SVs/assembly.fasta"
export FASTQ="/data/svi/prom/training/flowcell_10.4.1/GM23366/megaruptor_speed15/raw_data/fastq_pass/GM23366_megaruptor15.fastq.gz"


########################################
## ALIGN FQ TO ASSEMBLY ################
#######################################
minimap2 -L -t 72 -ax map-ont $ASM $FASTQ | samtools view -hbS > fq2asm.bam 
samtools sort -@ 72 fq2asm.bam -o fq2asm.sorted.bam 
samtools index -@ 24 fq2asm.sorted.bam 
rm fq2asm.bam

