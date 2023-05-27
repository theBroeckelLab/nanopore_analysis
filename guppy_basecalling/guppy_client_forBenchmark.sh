export FAST5="/data/prom_data/guppy_benchmarking/fast5/GM23366_rpm4800"
export OUTDIR="/data/prom_data/guppy_benchmarking/fastq_output"

####400bps-FAST####
time guppy_basecall_client \
	-i $FAST5 \
	-s $OUTDIR/fastq_fast_400bps/ \
	-c dna_r10.4.1_e8.2_400bps_fast.cfg \
	-q 10000000 \
	--min_qscore 9 \
	--recursive \
	--compress_fastq \
	--trim_adapters \
	--do_read_splitting \
	--chunks_per_runner 208 \
	--port "ipc:///tmp/.guppy/5555"
####400bps-HAC####
time guppy_basecall_client \
	-i $FAST5 \
	-s $OUTDIR/fastq_hac_400bps/ \
	-c dna_r10.4.1_e8.2_400bps_hac.cfg \
	-q 10000000 \
	--min_qscore 9 \
	--recursive \
	--compress_fastq \
	--trim_adapters \
	--do_read_splitting \
	--chunks_per_runner 208 \
	--port "ipc:///tmp/.guppy/5555"
####400bps-SUP####
time guppy_basecall_client \
	-i $FAST5 \
	-s $OUTDIR/fastq_sup_400bps/ \
	-c dna_r10.4.1_e8.2_400bps_sup.cfg \
	-q 10000000 \
	--min_qscore 9 \
	--recursive \
	--compress_fastq \
	--trim_adapters \
	--do_read_splitting \
	--chunks_per_runner 208 \
	--port "ipc:///tmp/.guppy/5555"
####260bps-FAST####
time guppy_basecall_client \
	-i $FAST5 \
	-s $OUTDIR/fastq_fast_260bps/ \
	-c dna_r10.4.1_e8.2_260bps_fast.cfg \
	-q 10000000 \
	--min_qscore 9 \
	--recursive \
	--compress_fastq \
	--trim_adapters \
	--do_read_splitting \
	--chunks_per_runner 208 \
	--port "ipc:///tmp/.guppy/5555"
####260bps-HAC####
time guppy_basecall_client \
	-i $FAST5 \
	-s $OUTDIR/fastq_hac_260bps/ \
	-c dna_r10.4.1_e8.2_260bps_hac.cfg \
	-q 10000000 \
	--min_qscore 9 \
	--recursive \
	--compress_fastq \
	--trim_adapters \
	--do_read_splitting \
	--chunks_per_runner 208 \
	--port "ipc:///tmp/.guppy/5555"
####260bps-SUP####
time guppy_basecall_client \
	-i $FAST5 \
	-s $OUTDIR/fastq_sup_260bps/ \
	-c dna_r10.4.1_e8.2_260bps_sup.cfg \
	-q 10000000 \
	--min_qscore 9 \
	--recursive \
	--compress_fastq \
	--trim_adapters \
	--do_read_splitting \
	--chunks_per_runner 208 \
	--port "ipc:///tmp/.guppy/5555"
