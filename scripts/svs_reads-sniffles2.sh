

export REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
export FASTQ="/data/svi/prom/mcw_svi/flowcell_10.4.1/mcw_svi_0042/raw_data/fastq_pass/mcw_svi_0042.fastq.gz"
export BAM="/data/svi/prom/mcw_svi/flowcell_10.4.1/mcw_svi_0042/alignment/mcw_svi_0042.sorted.bam"

###############
## sniffles2 ##
###############
#conda run -n sniffles
#conda activate sniffles
mkdir sniffles2

##default parameters
sniffles \
	--input $BAM \
	--reference $REF \
	--vcf sniffles2/output.vcf \
	--threads 16 \
 	--minsupport auto \
 	--minsvlen 35 \
  	--minsvlen-screen-ratio 0.9 \
  	--mapq 25 \
  	--qc-stdev True \
  	--qc-stdev-abs-max 500 \
  	--qc-coverage 1 \
  	--long-ins-length 2500 \
  	--long-del-length 50000 \
  	--long-del-coverage 0.66 \
  	--long-dup-length 50000 \
  	--long-dup-coverage 1.33 \
	--max-splits-kb 0.1 \
 	--max-splits-base 3 \
  	--min-alignment-length 1000 \
  	--phase-conflict-threshold 0.1 \
  	--detect-large-ins True \
  	--cluster-binsize 100 \
  	--cluster-r 2.5 \
  	--cluster-repeat-h 1.5 \
  	--cluster-repeat-h-max 1000 \
  	--cluster-merge-pos 150 \
  	--cluster-merge-len 0.33 \
  	--cluster-merge-bnd 1500 \
  	--genotype-ploidy 2 \
  	--genotype-error 0.05 \
  	--sample-id SAMPLE_ID \
	--allow-overwrite \
	--symbolic \
  	--max-del-seq-len 50000 
## default False parameters: --non-germline, --phase, --no-qc, --qc-strand, --output-rnames, --no-consensus, --no-sort, --no-progress, --quiet, --symbolic, --allow-overwrite



##modified parameters
#sniffles \
#	--input $BAM \
#	--reference $REF \
#	--vcf sniffles2/output.vcf \
#	--threads 48 \
#  	--minsupport 10 \
#  	--minsvlen 1500 \
#  	--minsvlen-screen-ratio 0.9 \
#  	--mapq 20 \
#  	--qc-stdev True \
#  	--qc-stdev-abs-max 500 \
#  	--qc-coverage 1 \
#  	--long-ins-length 2500 \
#  	--long-del-length 50000 \
#  	--long-del-coverage 0.66 \
#  	--long-dup-length 50000 \
#  	--long-dup-coverage 1.33 \
# 	--max-splits-kb 0.1 \
# 	--max-splits-base 3 \
#  	--min-alignment-length 1000 \
#  	--phase-conflict-threshold 0.1 \
#  	--detect-large-ins True \
#  	--cluster-binsize 100 \
#  	--cluster-r 2.5 \
#  	--cluster-repeat-h 1.5 \
#  	--cluster-repeat-h-max 1000 \
#  	--cluster-merge-pos 150 \
#  	--cluster-merge-len 0.33 \
#  	--cluster-merge-bnd 1500 \
#  	--genotype-ploidy 2 \
#  	--genotype-error 0.05 \
#  	--sample-id SAMPLE_ID \
#	--allow-overwrite \
#	--symbolic \
#  	--max-del-seq-len 50000 
## default False parameters: --non-germline, --phase, --no-qc, --qc-strand, --output-rnames, --no-consensus, --no-sort, --no-progress, --quiet, --symbolic, --allow-overwrite

