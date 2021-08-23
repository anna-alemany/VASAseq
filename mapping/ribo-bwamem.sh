#!/bin/bash
# Part 2 of the mapping pipeline: 
# BWA mapping to Ribosomal reference reads. 
# Produces  an output fastq.gz file containing non-ribosomal reads;
# and a bam file containing ribosomal reads. This last file contains also
# reads that map with very low qualities. 

if [ $# -ne 7 ]
then
    echo "Please, give:"
    echo "(1) reference fasta file (mouse, human or full path)"
    echo "(2) fastq file to map"
    echo "(3) prefix output"
    echo "(4) Path to bwa software"
    echo "(5) Path to samtools"
    echo "(6) Stranded protocol? (y/n)"
    echo "(7) Path to riboread-selection.py script"
    exit
fi

ref=$1
fq=$2
out=$3
p2bwa=$4
p2samtools=$5
stranded=$6
p2s=$7

folder=${fq%/*}
fq=${fq#*/}

# mapping short reads
${p2bwa}/bwa aln ${ref} ${folder}/${fq} > ${folder}/aln_${fq%.f*q}.sai 
${p2bwa}/bwa samse ${ref} ${folder}/aln_${fq%.f*q}.sai ${folder}/${fq} | ${p2samtools}/samtools view -Sb > ${out}.aln-ribo.bam &

# mapping normal reads
${p2bwa}/bwa mem -t 8 -h 15 ${ref} ${folder}/${fq} | ${p2samtools}/samtools view -Sb > ${out}.mem-ribo.bam & 

wait

${p2samtools}/samtools merge -n -r -h ${out}.aln-ribo.bam --threads 8 ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam 
rm ${out}.aln-ribo.bam ${out}.mem-ribo.bam ${folder}/aln_${fq%.f*q}.sai
${p2samtools}/samtools sort -n --threads 8 ${out}.all-ribo.bam -O BAM -o ${out}.nsorted.all-ribo.bam
rm ${out}.all-ribo.bam

${p2s}/riboread-selection.py ${out}.nsorted.all-ribo.bam $stranded ${out}

exit

#if [ $stranded == "n" ]
#then
#    # select reads that don't map and create a new fastq file
#    ${p2samtools}/samtools view -f 4 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq & 
#    # select reads that map (even low qualities, and create a new bam file
#    ${p2samtools}/samtools view -F 4 -Sb ${out}.ribo.sam > ${out}.ribo.bam & 
#elif [ $stranded == "y" ]
#then
#    ${p2samtools}/samtools view -f 16 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq
#    ${p2samtools}/samtools view -f 4 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq &
#    ${p2samtools}/samtools view -f 0 -F 4 -Sb ${out}.ribo.sam > ${out}.ribo.bam &
#fi

#wait

# zip fastq file and delete sam file from first mapping
#gzip ${out}.nonRibo.fastq &
#rm ${out}.ribo.sam

#wait

