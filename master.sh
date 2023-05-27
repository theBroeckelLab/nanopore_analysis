
#################################
## EXPORT FILE NAME AND REF #####
#################################
export PREFIX="sample"
export REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
export FASTQ="/data/svi/prom/mcw_svi/flowcell_10.4.1/mcw_svi_0042/raw_data/fastq_pass/sample.fastq.gz"



######################
### READ ALIGNMENT ###
######################
minimap2 -L -t 96 -ax map-ont $REF $FASTQ | samtools view -bS -@ 96 > sample.bam 
samtools sort -@ 96 sample.bam -o ${PREFIX}.sorted.bam 
samtools index -@ 16 ${PREFIX}.sorted.bam 
rm sample.bam


###########################
### ALIGNMENT QC ##########
###########################
pycoQC -f sequencing_summary_PAK63063_a1d9c21d_a1480da3.txt -a ${PREFIX}.sorted.bam -o pycoQC_alignment.html   
  


#######################################
### SV CALLING FROM READS WITH SVIM ###
#######################################
svim alignment \
	$OUTDIR \
	$BAM \
	$REF \
	--min_mapq 20 \
	--min_sv_size 1500 \
	--max_sv_size 1000000 \
	--segment_gap_tolerance 10 \
 	--segment_overlap_tolerance 5 \
	--partition_max_distance 1000 \
	--position_distance_normalizer 900 \
 	--edit_distance_normalizer 1.0 \
 	--cluster_max_distance 0.5 \
 	--del_ins_dup_max_distance 1.0 \
 	--trans_sv_max_distance 500 \
 	--max_consensus_length 10000 \
  	--minimum_score 10 \
 	--minimum_depth 10 \
  	--homozygous_threshold 0.8 \
  	--heterozygous_threshold 0.2 \
	--symbolic_alleles \
  	--sample Sample 



###################################
## ASSEMBLY WITH FLYE ############
##################################
/usr/bin/flye --nano-hq $FASTQ -i 1 --out-dir ./ --threads 96


#############################
## ASSEMBLY QC ##############
#############################
python3 /usr/local/src/quast-5.2.0/quast.py assembly.fasta -r $REF --large -o assembly_qc
assess_assembly -i assembly.fasta -r $REF -t 16 -p pomoxis_out


########################################
## ALIGN FQ TO NEW ASSEMBLY ############
#######################################
##align fq to new assembly
minimap2 -L -t 72 -ax map-ont $ASM $FASTQ | samtools view -hbS > fq2asm.bam 
samtools sort -@ 72 fq2asm.bam -o fq2asm.sorted.bam 
samtools index -@ 24 fq2asm.sorted.bam 
rm fq2asm.bam


###############################################
## SEGREGATE HAPLOID ASSEMBLY TO DIPLOID ######
###############################################
export INPUT_DIR="$PWD"
time docker run -v $INPUT_DIR:$INPUT_DIR -u `id -u`:`id -g` mkolmogo/hapdup:0.11 hapdup --assembly $INPUT_DIR/assembly.fasta --bam $INPUT_DIR/fq2asm.sorted.bam --out-dir $INPUT_DIR/hapdup -t 48 --rtype ont



###############################################
### SV CALLING FROM ASSEMBLY WITH SVIM-ASM ####
###############################################
##hap1 alignment
minimap2 -a -x asm5 --cs -r2k -t 48 $REF ./hapdup_dual_1.fasta > ./hapdup_1_asm2ref.sam
samtools sort -m4G -@ 24 -o ./hapdup_1_asm2ref.sorted.bam ./hapdup_1_asm2ref.sam
samtools index ./hapdup_1_asm2ref.sorted.bam
rm ./hapdup_1_asm2ref.sam 
##hap2 alignment
minimap2 -a -x asm5 --cs -r2k -t 48 $REF ./hapdup_dual_2.fasta > ./hapdup_2_asm2ref.sam
samtools sort -m4G -@ 24 -o ./hapdup_2_asm2ref.sorted.bam ./hapdup_2_asm2ref.sam
samtools index ./hapdup_2_asm2ref.sorted.bam
rm ./hapdup_2_asm2ref.sam 
#run svim-asm
svim-asm diploid ./hapdup/svim_asm \
	./hapdup_1_asm2ref.sorted.bam \
	./hapdup_2_asm2ref.sorted.bam \
	$REF \
	--min_mapq 20 \
	--min_sv_size 40 \
	--max_sv_size 100000 \
	--query_gap_tolerance 50 \
	--query_overlap_tolerance 50 \
	--reference_gap_tolerance 50 \
	--reference_overlap_tolerance 50

