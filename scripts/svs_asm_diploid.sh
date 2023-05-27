export REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"


###############################################
### SV CALLING - POST ASSEMBLY (SVIM-ASM) ####
##############################################
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


#run sv calling
##default parameters#
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
##publcation parameters##
#svim-asm diploid ./svim_asm \
#	./hapdup/hapdup_1_asm2ref.sorted.bam \
#	./hapdup/hapdup_2_asm2ref.sorted.bam \
#	$REF \
#	--min_sv_size 20 \
#	--tandem_duplications_as_insertions \
#	--interspersed_duplications_as_insertions \
#	--reference_gap_tolerance 1000 \
#	--reference_overlap_tolerance 1000 \
#	--query_gap_tolerance 2000 \
#	--query_overlap_tolerance 2000 \
#	--max_edit_distance 200 \
#	--query_names 
##modified parameters##
#svim-asm diploid ./svim_asm \
#	./hapdup/hapdup_1_asm2ref.sorted.bam \
#	./hapdup/hapdup_2_asm2ref.sorted.bam \
#	$REF \
#	--min_sv_size 1500 \
#	--max_sv_size 1000000 \
#	--min_mapq 20 \
#	--reference_gap_tolerance 1000 \
#	--reference_overlap_tolerance 1000 \
#	--query_gap_tolerance 2000 \
#	--query_overlap_tolerance 2000 \
#	--max_edit_distance 200 \
#	--query_names \
#	--symbolic_alleles

