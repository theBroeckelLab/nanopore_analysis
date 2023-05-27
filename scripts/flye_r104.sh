
#################################
## EXPORT FILE NAME AND REF #####
#################################
export REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
export FASTQ="/data/svi/prom/mcw_svi/flowcell_10.4.1/mcw_svi_0034/raw_data/fastq_pass/mcw_svi_0034.fastq.gz"

###################################
## ASSEMBLY WITH FLYE ############
##################################
/usr/bin/flye --nano-hq $FASTQ -i 1 --out-dir ./ --threads 96


################################
## POST-ASSEMBLY QC ############
################################
python3 /usr/local/src/quast-5.2.0/quast.py assembly.fasta -r $REF --large -o assembly_qc
python3 /data/svi/prom/scripts/fastmer.py --reference $REF --assembly assembly.fasta --min-mapping-quality 10 > ./assembly_qc/fastmerOut.txt
#assess_assembly -i assembly.fasta -r $REF -t 16 -p pomoxis_out

