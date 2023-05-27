export ASM="/data/svi/prom/training/flowcell_10.4.1/GM23366/megaruptor_speed15/assembly_SVs/assembly.fasta"
export FASTQ="/data/svi/prom/training/flowcell_10.4.1/GM23366/megaruptor_speed15/raw_data/fastq_pass/GM23366_megaruptor15.fastq.gz"


########################################
## SEGREGATE ASSEMBLY TO DIPLOID ######
#######################################
##run hapdup
export INPUT_DIR="$PWD"
time docker run -v $INPUT_DIR:$INPUT_DIR -u `id -u`:`id -g` mkolmogo/hapdup:0.11 hapdup --assembly $INPUT_DIR/assembly.fasta --bam $INPUT_DIR/fq2asm.sorted.bam --out-dir $INPUT_DIR/hapdup -t 48 --rtype ont


