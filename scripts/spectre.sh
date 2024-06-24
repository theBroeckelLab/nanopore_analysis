#!/bin/bash
##spectre needs mosdepth so run that 1st
##from general sample level
mkdir cnv
cd cnv/
mkdir mosdepth
cd mosdepth/
##need to activate the mosdepth conda enviroment
conda activate mosdepth
mosdepth -n -t 4 -x --by 1000 samplename_mosdepth sample.sorted.bam
conda deactivate
cd ..
##Run spectre tool
##If rerunning sample makesure to delete spectre output file or it will crash
##match bin sizes between mosdepth and spectre but that number can be changed
conda activate spectre
spectre CNVCaller  --bin-size 1000  --coverage mosdepth/  --sample-id sampleid --output-dir Spectre_results/  --reference /data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
conda deactivate spectre
