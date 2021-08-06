#!/bin/bash

inbam=$1
outfq=${inbam%.bam}.unmapped.fastq

samtools view -f 4 $inbam | awk '{OFS="\n"} {print "@"$1,$10,"+",$11}' > $outfq

gzip $outfq


