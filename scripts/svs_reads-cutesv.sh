
REF="/data/svi/prom/GM04927_081622_fc6/chr11_subset/GRCh38_chr11.fasta"
BAM="/data/svi/prom/GM04927_081622_fc6/chr11_subset/alignment/GM04927_chr11.sorted.bam"
PREFIX="GM04927_chr11"

####################
### cuteSV #########
#####################
mkdir cuteSV

##default parameters
cuteSV \
	$BAM \
	$REF \
	cuteSV/out.vcf \
	cuteSV \
	--max_split_parts 7 \
	--min_mapq 20 \
	--min_read_len 500 \
	--merge_del_threshold 0 \
	--merge_ins_threshold 100 \
	--min_support 10 \
	--min_size 30 \
	--max_size -1 \
	--min_siglength 10 \
	--max_cluster_bias_INS 100 \
	--diff_ratio_merging_INS 0.3  \
	--max_cluster_bias_DEL 100 \
	--diff_ratio_merging_DEL 0.3  \
	--max_cluster_bias_DUP 500 \
	--max_cluster_bias_INV 500 \
	--max_cluster_bias_TRA 50 \
	--diff_ratio_filtering_TRA 0.6 \
	--genotype \
	--threads 48



