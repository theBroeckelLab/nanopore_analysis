
guppy_basecaller -i /data/prom_data/Covaris_GM23366_Opt_030823/GM23366_4800rpm/20230308_2031_1E_PAK63203_e821227d/basecall_testing/fast5_subset/ -s /data/prom_data/Covaris_GM23366_Opt_030823/GM23366_4800rpm/20230308_2031_1E_PAK63203_e821227d/basecall_testing/fast5_subset/basecalling_commandLine -c dna_r10.4.1_e8.2_400bps_sup.cfg -q 0 --min_qscore 9 -x "auto" --recursive --compress_fastq --verbose_logs --trim_adapters --do_read_splitting --chunks_per_runner 256


