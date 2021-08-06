#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give (1) input fastq file; (2) output directory; (3) path2trimgalore; (4) path2cutadapt;"
  exit
fi

file2trim=$1
outdir=$2
path2trimgalore=$3
path2cutadapt=$4

# trim adaptors
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt ${file2trim} -o ${outdir}
mv ${file2trim}_trimming_report.txt ${file2trim%.fastq.gz}_trimming_report.txt

# trim homopolymers
${path2cutadapt}/cutadapt -m 15 --trim-n -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" -o ${file2trim%.fastq.gz}_trimmed_homoATCG.fq.gz  ${file2trim%.fastq.gz}_trimmed.fq.gz

exit





#java -jar ${path2trimmomatic}/trimmomatic-0.36.jar SE -phred33 ${file2trim%.fastq.gz}_trimmed.fq.gz ${file2trim%.fastq.gz}_trimHomo.fq ILLUMINACLIP:/hpc/hub_oudenaarden/fsalmen/barcodes/homopolymers.fa:0::0:0 MINLEN:20

#gzip ${file2trim%.fastq.gz}_trimHomo.fq


