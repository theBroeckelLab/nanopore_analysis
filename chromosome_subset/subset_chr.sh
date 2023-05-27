##### Subset a chromosome from ONT data ######
REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
BAM="/data/svi/prom/mcw_svi/flowcell_10.4.1/mcw_svi_0034/alignment/mcw_svi_0034.sorted.bam"

## extract chr of interest from GRCh38 reference
samtools faidx $REF chrX > GRCh38_chrX.fasta
samtools faidx GRCh38_chrX.fasta

## pull chr of interest from sorted bam file
samtools view -hb $BAM -@24 chrX > sample_chrX.sorted.bam
samtools index sample_chrX.sorted.bam

## convert bam to fastq
samtools fastq -@24 sample_chrX.sorted.bam > sample_chrX.fastq



