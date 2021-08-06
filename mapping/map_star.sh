#!/bin/bash

if [ $# -ne 5 ]
then
    echo "Please, give as input:"
    echo "1) path to STAR"
    echo "2) path to samtools"
    echo "3) genome folder"
    echo "4) input fastq file (only one)"
    echo "5) prefix for output file names"
    exit
fi

p2star=$1
p2samtools=$2
genome=$3
inputfq=$4
outprefix=$5

#${p2star}/STAR --runThreadN 8 --genomeDir ${genome} --readFilesIn ${inputfq} --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted SortedByCoordinate --quantMode TranscriptomeSAM --outSAMattributes All  --outFileNamePrefix ${outprefix}
#${p2star}/STAR --runThreadN 8 --genomeDir ${genome} --readFilesIn ${inputfq} --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted SortedByCoordinate --outSAMattributes All  --outFileNamePrefix ${outprefix}

${p2star}/STAR --runThreadN 8 --genomeDir ${genome} --readFilesIn ${inputfq} --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${outprefix}


rm -r ${outprefix}_STARtmp
rm ${outprefix}Log.progress.out

mv ${outprefix}Log.out ${outprefix}Log.txt
mv ${outprefix}Log.final.out ${outprefix}Log.final.txt

#${p2samtools}/samtools view -q 255 ${outprefix}Aligned.sortedByCoord.out.bam -b -o ${outprefix}Aligned.sortedByCoord.filtered.bam
#rm ${outprefix}Aligned.sortedByCoord.out.bam


