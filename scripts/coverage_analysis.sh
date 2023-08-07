#!/bin/bash
#Take genome file, and use bedtools to make into desired bins
#Run these commands on the  bins .bed
#sed -i -r 's/\s+/:/' "binsize".bed
#sed -i -r 's/\s+/-/' "binsize".bed

#Use this new bed file as input 1

inbed=$1
#Input 2 is desired bam to get coverage of bins
bam=$2
#Name for output file
output=$3

echo 'rname	startpos	endpos	numreads	covbases	 coverage	meandepth	meanbaseq	meanmapq' >> ${output}.txt

#Produces the coverage of the given bin regions by looping the read lines into samtools coverage and appending to output file
while read -r line;
do
	samtools coverage -H -r $line $bam >> ${output}.txt 
done < $inbed
